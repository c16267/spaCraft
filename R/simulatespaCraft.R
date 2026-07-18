#' @title Simulate Spatial Transcriptomics Data (Hybrid Mode, Optimized)
#' @description
#' Optimized rewrite of \code{simulatespaCraft}. Statistically identical to the
#' previous version under the same \code{seed_base}, but substantially faster.
#' The main algorithmic gains include:
#'   \itemize{
#'     \item The conditional-texture (\code{eta_cond}) KNN structure is built
#'       once per domain per sample and reused across all genes.
#'     \item Sample/group invariants are hoisted out of the per-sample and per-gene loops.
#'     \item The per-spot domain-mean lookup is fully vectorized.
#'   }
#'
#' @param object A \code{spaCraft} object.
#' @param n_sample_per_group Integer. Number of synthetic samples to generate per group.
#' @param scenario_settings List. Scenario configurations including \code{DE_lfc},
#'   \code{target_prop_case}, \code{kappa_m}, \code{kappa_bio}, \code{delta_rho},
#'   and \code{lambda_cond}.
#' @param target_domain Character. Target domain name for differential effects.
#' @param genes_to_simulate Character vector. Genes to simulate.
#' @param de_genes Character vector. Subset of genes that will have DE effects applied.
#' @param verbose Logical. Print progress.
#' @param n_cores Integer. Number of cores for parallel execution.
#' @param gp_rank Integer. Rank for Gaussian Process (legacy).
#' @param gp_jitter Numeric. Jitter to ensure positive definiteness.
#' @param spatial_mode Character. "original" or "basis".
#' @param graph_k Integer. KNN for the whole-tissue spatial graph.
#' @param cond_k_nn Integer. KNN for the pilot-guided conditional texture.
#' @param basis_rank Integer. Basis rank for spatial mode "basis".
#' @param basis_kernel Character. "gaussian" or "exponential".
#' @param basis_knot_method Character. "sample" or "grid".
#' @param basis_seed Integer. Seed for basis knot sampling.
#' @param lambda_cloud Optional numeric in \code{[0, 1]}. Used for FGKMM simulation.
#' @param seed_base Optional integer for reproducible simulations.
#'
#' @return A list of simulated synthetic samples.
#' @importFrom stats plogis rnorm rbinom setNames median cov
#' @importFrom utils modifyList
#' @export
#'
simulatespaCraft <- function(object,
                               n_sample_per_group = 10,
                               scenario_settings = list(
                                 DE_lfc            = 0,
                                 target_prop_case  = NULL,
                                 delta_pp          = NULL,
                                 kappa_m           = 1,
                                 kappa_bio         = 1,
                                 delta_rho         = 1,
                                 lambda_cond       = 0.1
                               ),
                               target_domain     = NULL,
                               genes_to_simulate = NULL,
                               de_genes          = NULL,
                               verbose           = TRUE,
                               n_cores           = 4,
                               gp_rank           = 20L,
                               gp_jitter         = 1e-8,
                               spatial_mode      = c("original", "basis"),
                               graph_k           = 20L,
                               cond_k_nn         = 20L,
                               basis_rank        = 25L,
                               basis_kernel      = c("gaussian", "exponential"),
                               basis_knot_method = c("sample", "grid"),
                               basis_seed        = 1L,
                               lambda_cloud      = NULL,
                               seed_base         = NULL) {

  spatial_mode      <- match.arg(spatial_mode)
  basis_kernel      <- match.arg(basis_kernel)
  basis_knot_method <- match.arg(basis_knot_method)

  if (!requireNamespace("FNN", quietly = TRUE)) stop("Package 'FNN' required.")

  if (.Platform$OS.type == "windows" && n_cores > 1L) {
    if (verbose) message(">>> Note: parallelisation limited to 1 core on Windows.")
    n_cores <- 1L
  }
  n_cores <- max(1L, as.integer(n_cores))
  cond_k_nn <- max(1L, as.integer(cond_k_nn))

  # CRAN-safe RNG state management for 1-core scenarios
  if (!is.null(seed_base) && is.finite(seed_base)) {
    if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
      .old_seed <- get(".Random.seed", envir = globalenv(), inherits = FALSE)
      on.exit(assign(".Random.seed", .old_seed, envir = globalenv()), add = TRUE)
    }
  }

  defaults <- list(
    DE_lfc           = 0,
    target_prop_case = NULL,
    delta_pp         = NULL,
    kappa_m          = 1,
    kappa_bio        = 1,
    delta_rho        = 1,
    lambda_cond      = 0.5,
    sd_logit_re      = 0
  )
  scenario <- utils::modifyList(defaults, scenario_settings)

  params_comp <- object@params_composition
  expr_slot   <- object@params_expression

  if (is.null(target_domain)) target_domain <- params_comp$beta_binomial$target
  ref_domain <- params_comp$beta_binomial$reference

  # ===========================================================================
  # Regime detection
  # ===========================================================================
  K_regime     <- expr_slot$K_regime %||% "multi"
  is_singleton <- identical(K_regime, "singleton")
  if (verbose) message(sprintf(">>> Detected K_regime = '%s'", K_regime))

  # ===========================================================================
  # Geometry structure detection + per-group accessor
  # ===========================================================================
  .detect_geometry_layout <- function(params_geometry) {
    if (is.null(params_geometry) || length(params_geometry) == 0L) return("empty")
    top_names <- names(params_geometry)
    if (is.null(top_names)) return("empty")
    if (all(top_names %in% c("0", "1"))) return("grouped")
    first <- params_geometry[[1]]
    if (is.list(first) && !is.null(first$placement)) return("flat")
    if (is.list(first) && length(first) > 0L) {
      sub_first <- first[[1]]
      if (is.list(sub_first) && !is.null(sub_first$placement)) return("grouped")
    }
    return("grouped")
  }

  geom_layout <- .detect_geometry_layout(object@params_geometry)
  if (verbose) message(sprintf(">>> Geometry layout: '%s'", geom_layout))

  .get_geometry_for_group <- function(grp_chr) {
    pg <- object@params_geometry
    if (geom_layout == "flat") return(pg)
    .extract_geometry_group(pg, grp_chr)
  }

  # ===========================================================================
  # Expression adapter for singleton regime
  # ===========================================================================
  .adapt_singleton_expr <- function(grp_chr) {
    grp_data <- expr_slot[[grp_chr]]
    if (is.null(grp_data) || length(grp_data) == 0L) return(list())

    adapted <- list()
    for (g in names(grp_data)) {
      p_ex <- grp_data[[g]]
      if (is.null(p_ex)) next

      mu_grand <- p_ex$mu_grand
      if (is.null(mu_grand) && !is.null(p_ex$mu_domain)) mu_grand <- p_ex$mu_domain

      th <- p_ex$theta
      alpha_val <- if (!is.null(th) && "alpha" %in% names(th)) unname(th["alpha"])
      else if (!is.null(th)) unname(th[1]) else 0.2
      rho_val   <- if (!is.null(th) && "rho" %in% names(th)) unname(th["rho"])
      else if (!is.null(th)) unname(th[2]) else 10
      nugget_val <- if (!is.null(th) && "sigma.sq" %in% names(th)) unname(th["sigma.sq"])
      else if (!is.null(th)) unname(th[3]) else 0.05

      if (!is.finite(alpha_val))  alpha_val  <- 0.2
      if (!is.finite(rho_val))    rho_val    <- 10
      if (!is.finite(nugget_val)) nugget_val <- 0.05

      theta <- c(alpha = alpha_val, rho = rho_val, sigma.sq = nugget_val)
      sigma_bio <- if (!is.null(p_ex$sigma_bio) && is.finite(p_ex$sigma_bio)) p_ex$sigma_bio else 0

      adapted[[g]] <- list(mu_grand = mu_grand, theta = theta, sigma_bio = sigma_bio)
    }
    adapted
  }

  # ===========================================================================
  # Geometry bounds for singleton regime
  # ===========================================================================
  .bounds_from_global_theta <- function(global_theta) {
    M <- length(global_theta)
    if (M == 0L) return(list(r_min = 0.01, r_max = 1.0, sigma_min = 1e-6,
                             sigma_max = 1.0, tau_min = 0.1, tau_max = 100))

    get_r_mean <- function(th) {
      if (!is.null(th$mu_r) && is.finite(th$mu_r)) return(as.numeric(th$mu_r))
      NA_real_
    }
    get_r_sd <- function(th) {
      if (!is.null(th$var_r) && is.finite(th$var_r) && th$var_r > 0) return(sqrt(th$var_r))
      mu_r <- if (!is.null(th$mu_r)) as.numeric(th$mu_r) else 0.1
      0.15 * abs(mu_r)
    }
    get_tau_mean <- function(th) {
      if (!is.null(th$mu_tau) && is.finite(th$mu_tau)) return(as.numeric(th$mu_tau))
      if (!is.null(th$b_tau) && is.finite(th$b_tau) && th$b_tau > 0 &&
          !is.null(th$a_tau) && is.finite(th$a_tau)) return(as.numeric(th$a_tau / th$b_tau))
      1
    }
    get_sigma_mean <- function(th) {
      if (!is.null(th$mu_sigma2) && is.finite(th$mu_sigma2) && th$mu_sigma2 > 0)
        return(as.numeric(th$mu_sigma2))
      if (!is.null(th$a_sigma) && is.finite(th$a_sigma) && th$a_sigma > 1 &&
          !is.null(th$b_sigma) && is.finite(th$b_sigma))
        return(as.numeric(th$b_sigma / (th$a_sigma - 1)))
      if (!is.null(th$b_sigma) && is.finite(th$b_sigma)) return(as.numeric(th$b_sigma))
      1e-4
    }

    r_vals      <- vapply(global_theta, get_r_mean,     numeric(1))
    r_sds       <- vapply(global_theta, get_r_sd,       numeric(1))
    tau_means   <- vapply(global_theta, get_tau_mean,   numeric(1))
    sigma_means <- vapply(global_theta, get_sigma_mean, numeric(1))

    r_vals_clean <- r_vals[is.finite(r_vals)]
    if (length(r_vals_clean) == 0L) r_vals_clean <- 0.1

    list(
      r_min     = max(min(r_vals_clean - 3 * r_sds, na.rm = TRUE), 0.001),
      r_max     = max(max(r_vals_clean + 3 * r_sds, na.rm = TRUE), 0.01),
      sigma_min = max(min(sigma_means, na.rm = TRUE) * 0.1, 1e-8),
      sigma_max = max(max(sigma_means, na.rm = TRUE) * 10, 1e-6),
      tau_min   = max(min(tau_means, na.rm = TRUE) * 0.1, 0.01),
      tau_max   = max(max(tau_means, na.rm = TRUE) * 10, 1)
    )
  }

  # ===========================================================================
  # Composition parameters
  # ===========================================================================
  beta0 <- params_comp$beta_binomial$beta0
  beta1 <- params_comp$beta_binomial$beta1 %||% 0
  p0    <- stats::plogis(beta0)
  OR_pilot <- exp(beta1)

  has_tpc   <- !is.null(scenario$target_prop_case)
  has_delta <- !is.null(scenario$delta_pp)
  if (has_tpc && has_delta)
    stop("Specify either 'target_prop_case' or 'delta_pp' in scenario_settings, not both.")

  if (has_tpc) {
    target_prop_case <- as.numeric(scenario$target_prop_case)
    if (target_prop_case <= 0 || target_prop_case >= 1)
      stop("'target_prop_case' must be in the open interval (0, 1).")
    OR_target <- (target_prop_case / (1 - target_prop_case)) / (p0 / (1 - p0))
    theta_DA  <- OR_target / OR_pilot
  } else if (has_delta) {
    delta_pp         <- as.numeric(scenario$delta_pp)
    target_prop_case <- p0 + delta_pp / 100
    if (target_prop_case <= 0 || target_prop_case >= 1)
      stop(sprintf("'delta_pp' = %.2f implies target_prop_case = %.4f, outside (0, 1).",
                   delta_pp, target_prop_case))
    OR_target <- (target_prop_case / (1 - target_prop_case)) / (p0 / (1 - p0))
    theta_DA  <- OR_target / OR_pilot
  } else {
    theta_DA         <- 1
    target_prop_case <- stats::plogis(beta0 + beta1)
  }

  OR_syn <- OR_pilot * theta_DA
  p1     <- stats::plogis(beta0 + log(OR_syn))

  if (verbose) {
    message(sprintf(
      paste0(">>> DA effect conversion:\n",
             "    pilot p0 = %.4f, pilot OR = %.4f (beta1 = %.4f)\n",
             "    user input  -> target_prop_case = %.4f (delta_pp = %+.2f pp)\n",
             "    backend     -> theta_DA = %.4f, OR_syn = %.4f\n",
             "    final p0 = %.4f, p1 = %.4f"),
      p0, OR_pilot, beta1, target_prop_case, 100 * (target_prop_case - p0),
      theta_DA, OR_syn, p0, p1))
  }

  # ===========================================================================
  # Gene list
  # ===========================================================================
  reserved_keys <- c("top_genes", "stable_genes", "stable_gene_stats",
                     "K_regime", "calibration", "settings", "by_group", "shared")
  group_keys <- intersect(setdiff(names(expr_slot), reserved_keys), c("0", "1"))
  available_genes <- sort(unique(unlist(lapply(group_keys, function(k) names(expr_slot[[k]])))))

  if (verbose)
    message(sprintf(">>> Available genes: %d (from groups: %s)",
                    length(available_genes), paste(group_keys, collapse = ", ")))

  if (is.null(genes_to_simulate)) genes_to_simulate <- available_genes
  genes_to_simulate <- intersect(genes_to_simulate, available_genes)
  if (is.null(de_genes)) de_genes <- genes_to_simulate
  de_genes <- intersect(de_genes, genes_to_simulate)
  de_set   <- as.character(de_genes)

  groups          <- c(rep(0, n_sample_per_group), rep(1, n_sample_per_group))
  n_total_samples <- length(groups)

  pilot_group_idx <- list(
    "0" = which(vapply(object@pilot_data, function(x) as.character(x$group) == "0", logical(1))),
    "1" = which(vapply(object@pilot_data, function(x) as.character(x$group) == "1", logical(1)))
  )
  if (length(pilot_group_idx[["0"]]) == 0L) pilot_group_idx[["0"]] <- seq_along(object@pilot_data)
  if (length(pilot_group_idx[["1"]]) == 0L) pilot_group_idx[["1"]] <- seq_along(object@pilot_data)

  # ===========================================================================
  # HOISTED INVARIANTS
  # ===========================================================================
  pilot_total_mean <- mean(vapply(object@pilot_data, function(x) nrow(x$coords), numeric(1)), na.rm = TRUE)
  Fixed_N <- round(pilot_total_mean)

  geom_by_group <- list("0" = .get_geometry_for_group("0"),
                        "1" = .get_geometry_for_group("1"))
  expr_by_group <- list(
    "0" = if (is_singleton) .adapt_singleton_expr("0") else expr_slot[["0"]],
    "1" = if (is_singleton) .adapt_singleton_expr("1") else expr_slot[["1"]]
  )

  pilot_coords_norm_all <- lapply(object@pilot_data, function(ps)
    as.matrix(.normalize_coords_unit(ps$coords)[, c("x", "y")]))
  pilot_dom_all <- lapply(object@pilot_data, function(ps)
    if (!is.null(ps$coords) && "domain" %in% names(ps$coords))
      as.character(ps$coords$domain) else character(0))

  scenario_de_shift  <- { v <- as.numeric(scenario$DE_lfc);      if (is.finite(v)) v else 0 }
  scenario_kappa_bio <- { v <- as.numeric(scenario$kappa_bio %||% 1); if (is.finite(v) && v >= 0) v else 1 }
  scenario_lambda    <- { v <- as.numeric(scenario$lambda_cond %||% 0); if (is.finite(v)) max(0, min(1, v)) else 0 }
  scenario_delta_rho <- as.numeric(scenario$delta_rho)

  # --------------------------------------------------------------------------
  # Per-sample worker
  # --------------------------------------------------------------------------
  .simulate_one_sample <- function(k) {
    if (!is.null(seed_base) && is.finite(seed_base)) set.seed(as.integer(seed_base) + k)

    grp_id  <- groups[k]
    grp_chr <- as.character(grp_id)
    is_case <- (grp_id == 1L)
    samp_id <- sprintf("Sim_Grp%d_%03d", grp_id, k)

    # 3.1: Domain composition
    m_k <- min(
      max(10L, round(params_comp$contrast_size$by_group[[grp_chr]]$mu *
                       as.numeric(scenario$kappa_m %||% 1))),
      Fixed_N
    )

    eta_base <- if (is_case) (beta0 + log(OR_syn)) else beta0
    sd_re    <- as.numeric(scenario$sd_logit_re %||% 0)
    eta_k    <- eta_base + if (sd_re > 0) stats::rnorm(1, 0, sd_re) else 0
    mu_c     <- max(min(stats::plogis(eta_k), 1 - 1e-5), 1e-5)
    n_target <- stats::rbinom(1, m_k, mu_c)
    n_ref    <- m_k - n_target

    bg_props   <- params_comp$background$proportions[[grp_chr]]
    n_bg_total <- Fixed_N - m_k

    if (n_bg_total > 0L && length(bg_props) > 0L) {
      bg_props <- bg_props / sum(bg_props)
      n_bg_vec <- floor(n_bg_total * bg_props)
      leftover <- n_bg_total - sum(n_bg_vec)
      if (leftover > 0L) {
        idx <- order(n_bg_total * bg_props - n_bg_vec, decreasing = TRUE)[seq_len(leftover)]
        n_bg_vec[idx] <- n_bg_vec[idx] + 1L
      }
    } else {
      n_bg_vec <- stats::setNames(rep(0L, length(bg_props)), names(bg_props))
    }

    dom_counts <- c(stats::setNames(n_target, target_domain),
                    stats::setNames(n_ref, ref_domain),
                    n_bg_vec)

    # 3.2: Spatial geometry
    params_geom_grp <- geom_by_group[[grp_chr]]
    sim_coords_list <- list()
    for (dom in names(dom_counts)) {
      n_d <- as.integer(dom_counts[[dom]])
      if (n_d <= 0L) next
      geom_p   <- params_geom_grp[[dom]]
      coords_d <- NULL

      if (!is.null(geom_p) && !is.null(geom_p$pattern$pooled)) {
        pattern <- geom_p$pattern$pooled
        if (!is.null(geom_p$fit_data) && length(geom_p$fit_data) > 0L) {
          bounds <- get_r_sigma_bounds(lapply(geom_p$fit_data, function(x) x$theta),
                                       length(pattern$global_theta))
        } else {
          bounds <- .bounds_from_global_theta(pattern$global_theta)
        }

        coords_d <- tryCatch({
          res <- suppressMessages(simulate_from_FG_post(
            n_d, pattern$global_theta, pattern$alpha_hat,
            bounds$r_min, bounds$r_max, bounds$sigma_min, bounds$sigma_max,
            bounds$tau_min, bounds$tau_max, lambda = lambda_cloud %||% 0.5))
          ctr <- .sample_domain_centroid(geom_p)
          res$x <- res$x + ctr[1]; res$y <- res$y + ctr[2]
          bad <- !is.finite(res$x) | !is.finite(res$y)
          if (any(bad)) {
            res$x[bad] <- ctr[1] + stats::rnorm(sum(bad), 0, 0.01)
            res$y[bad] <- ctr[2] + stats::rnorm(sum(bad), 0, 0.01)
          }
          res
        }, error = function(e) NULL)
      }

      if (is.null(coords_d)) {
        ctr <- if (!is.null(geom_p))
          tryCatch(.sample_domain_centroid(geom_p), error = function(e) c(0.5, 0.5))
        else c(0.5, 0.5)
        coords_d <- data.frame(x = ctr[1] + stats::rnorm(n_d, 0, 0.05),
                               y = ctr[2] + stats::rnorm(n_d, 0, 0.05))
      }

      coords_d$domain        <- dom
      sim_coords_list[[dom]] <- coords_d
    }

    sim_coords_df  <- do.call(rbind, sim_coords_list)
    sim_coords_mat <- as.matrix(sim_coords_df[, c("x", "y")])

    na_mask <- !is.finite(sim_coords_mat[, 1]) | !is.finite(sim_coords_mat[, 2])
    if (any(na_mask)) {
      sim_coords_mat[na_mask, 1] <- stats::runif(sum(na_mask), 0, 1)
      sim_coords_mat[na_mask, 2] <- stats::runif(sum(na_mask), 0, 1)
      sim_coords_df$x <- sim_coords_mat[, 1]
      sim_coords_df$y <- sim_coords_mat[, 2]
    }
    n_cells <- nrow(sim_coords_mat)

    # Spatial cache
    spatial_cache <- build_spatial_cache(
      sim_coords_df, spatial_mode = spatial_mode, graph_k = graph_k,
      basis_rank = basis_rank, basis_knot_method = basis_knot_method,
      basis_seed = basis_seed)
    sim_domains <- spatial_cache$sim_domains

    # Pilot bundle
    pick_pool <- pilot_group_idx[[grp_chr]]
    sel_idx   <- sample(pick_pool, min(3L, length(pick_pool)))
    n_pilot   <- length(sel_idx)
    pilot_samples_sel <- object@pilot_data[sel_idx]
    pilot_coords_sel  <- pilot_coords_norm_all[sel_idx]
    pilot_dom_sel     <- pilot_dom_all[sel_idx]
    pilot_mix_w       <- rep(1 / n_pilot, n_pilot)

    params_expr_grp <- expr_by_group[[grp_chr]]

    domain_vec <- if ("tilde_d" %in% names(sim_coords_df))
      as.character(sim_coords_df$tilde_d) else as.character(sim_coords_df$domain)
    uniq_doms      <- unique(domain_vec)
    syn_idx_by_dom <- split(seq_len(n_cells), domain_vec)

    # eta_cond KNN cache
    cond_cache <- vector("list", n_pilot)
    if (scenario_lambda > 1e-8) {
      for (p in seq_len(n_pilot)) {
        pcoords <- pilot_coords_sel[[p]]
        pdom    <- pilot_dom_sel[[p]]
        if (is.null(pcoords) || nrow(pcoords) == 0L) { cond_cache[[p]] <- NULL; next }

        per_dom <- vector("list", length(uniq_doms)); names(per_dom) <- uniq_doms
        for (d in uniq_doms) {
          syn_idx_d <- syn_idx_by_dom[[d]]
          pilot_idx <- which(pdom == d)
          if (length(pilot_idx) < 3L) pilot_idx <- seq_along(pdom)
          if (length(pilot_idx) == 0L) next

          k_use <- max(1L, min(cond_k_nn, length(pilot_idx)))
          nn <- FNN::get.knnx(data  = pcoords[pilot_idx, , drop = FALSE],
                              query = sim_coords_mat[syn_idx_d, , drop = FALSE],
                              k     = k_use)
          d2 <- nn$nn.dist^2
          h  <- stats::median(d2[is.finite(d2)], na.rm = TRUE)
          if (!is.finite(h) || h <= 0) h <- mean(d2[is.finite(d2)], na.rm = TRUE)
          if (!is.finite(h) || h <= 0) h <- 1
          W  <- exp(-d2 / (h + 1e-8))
          rs <- rowSums(W); rs[!is.finite(rs) | rs <= 0] <- 1
          W  <- W / rs

          per_dom[[d]] <- list(syn_idx   = syn_idx_d,
                               pilot_idx = pilot_idx,
                               nn_index  = nn$nn.index,
                               W         = W)
        }
        cond_cache[[p]] <- list(pdom = pdom, per_dom = per_dom)
      }
    }

    # Per-gene generator
    .one_gene <- function(g) {
      p_ex <- if (!is.null(params_expr_grp)) params_expr_grp[[g]] else NULL
      if (is.null(p_ex)) {
        mu_grand  <- stats::setNames(rep(0.1, length(sim_domains)), sim_domains)
        theta     <- c(alpha = 0.2, rho = 10, sigma.sq = 0.05)
        sigma_bio <- 0
      } else {
        mu_grand  <- as.numeric(p_ex$mu_grand)
        names(mu_grand) <- if (!is.null(names(p_ex$mu_grand))) names(p_ex$mu_grand) else sim_domains
        theta     <- p_ex$theta
        sigma_bio <- if (!is.null(p_ex$sigma_bio)) p_ex$sigma_bio else 0
      }

      th       <- extract_theta_triplet(theta)
      alpha_sp <- as.numeric(th["alpha"]); rho_sp <- as.numeric(th["rho"])
      nugget2  <- as.numeric(th["nugget"])
      if (!is.finite(alpha_sp) || alpha_sp < 0) alpha_sp <- 0
      if (!is.finite(rho_sp)   || rho_sp <= 0)  rho_sp   <- 10
      if (!is.finite(nugget2)  || nugget2 < 0)  nugget2  <- 0

      rho_eff <- rho_sp
      if (is_case && is.finite(scenario_delta_rho) && scenario_delta_rho > 0)
        rho_eff <- rho_sp * scenario_delta_rho

      w   <- pilot_mix_w
      sel <- sample(seq_len(n_pilot), size = 1L, prob = w)

      mu_fallback <- mean(mu_grand[is.finite(mu_grand)], na.rm = TRUE)
      if (!is.finite(mu_fallback)) mu_fallback <- 0.1

      g_is_de <- g %in% de_set
      mu_log  <- numeric(n_cells)
      for (dom in uniq_doms) {
        base <- if (!is.null(names(mu_grand)) && dom %in% names(mu_grand)) mu_grand[[dom]] else mu_fallback
        if (!is.finite(base)) base <- 0.1
        if (is_case && identical(as.character(dom), as.character(target_domain)) && g_is_de)
          base <- base + scenario_de_shift
        mu_log[syn_idx_by_dom[[dom]]] <- base
      }

      sigma_bio_new <- as.numeric(sigma_bio) * scenario_kappa_bio
      if (!is.finite(sigma_bio_new) || sigma_bio_new < 0) sigma_bio_new <- 0
      B_shift <- stats::rnorm(1, mean = 0, sd = sqrt(sigma_bio_new))

      eta_cond <- numeric(n_cells)
      if (scenario_lambda > 1e-8 && !is.null(cond_cache[[sel]])) {
        pilot_p <- pilot_samples_sel[[sel]]
        cp      <- cond_cache[[sel]]

        gene_in_log <- !is.null(pilot_p$logcounts) && (g %in% rownames(pilot_p$logcounts))
        gene_in_cnt <- !is.null(pilot_p$counts)    && (g %in% rownames(pilot_p$counts))

        if (gene_in_log || gene_in_cnt) {
          p_log <- if (gene_in_log) as.numeric(pilot_p$logcounts[g, ])
          else             log1p(as.numeric(pilot_p$counts[g, ]))

          p_mu <- unname(mu_grand[cp$pdom])
          p_mu[!is.finite(p_mu)] <- mu_fallback

          p_res <- p_log - p_mu
          p_res[!is.finite(p_res)] <- 0
          p_res <- p_res - mean(p_res)

          for (d in uniq_doms) {
            cd <- cp$per_dom[[d]]
            if (is.null(cd)) next
            p_res_dom <- p_res[cd$pilot_idx]
            gathered  <- matrix(p_res_dom[cd$nn_index], nrow = length(cd$syn_idx))
            eta_cond[cd$syn_idx] <- rowSums(cd$W * gathered)
          }
          eta_cond[!is.finite(eta_cond)] <- 0
          eta_cond <- eta_cond - mean(eta_cond)
        }
      }

      eta_para <- numeric(n_cells)
      if (alpha_sp > 0) {
        if (spatial_mode == "original") {
          eta_para <- draw_graph_gp_whole_tissue(spatial_cache, alpha_sp, rho_eff, include_self = TRUE)
        } else {
          eta_para <- draw_basis_gp_whole_tissue(spatial_cache, alpha_sp, rho_eff, basis_kernel)
        }
      }
      eta_para[!is.finite(eta_para)] <- 0

      eta <- scenario_lambda * eta_cond + (1 - scenario_lambda) * eta_para
      eta[!is.finite(eta)] <- 0

      eps_noise <- if (nugget2 > 0) stats::rnorm(n_cells, 0, sqrt(nugget2)) else numeric(n_cells)

      log1p_y <- mu_log + B_shift + eta + eps_noise
      log1p_y[!is.finite(log1p_y)] <- 0
      log1p_y[log1p_y < 0]         <- 0
      log1p_y
    }

    n_genes   <- length(genes_to_simulate)
    log1p_mat <- matrix(0.0, nrow = n_genes, ncol = n_cells)
    for (gi in seq_len(n_genes)) log1p_mat[gi, ] <- .one_gene(genes_to_simulate[gi])

    rownames(log1p_mat) <- genes_to_simulate
    colnames(log1p_mat) <- sprintf("spot_%05d", seq_len(n_cells))

    count_mat <- log1p_mat
    count_mat[] <- pmax(0, round(expm1(log1p_mat)))

    list(counts    = count_mat,
         logcounts = log1p_mat,
         coords    = sim_coords_df,
         meta      = data.frame(sample_id = samp_id, group = grp_id, stringsAsFactors = FALSE))
  }

  # --------------------------------------------------------------------------
  # Execute
  # --------------------------------------------------------------------------
  if (verbose) message(sprintf(">>> Generating %d samples (continuous FGKMM, n_cores=%d)...", n_total_samples, n_cores))

  if (n_cores == 1L) {
    simulated_data <- lapply(seq_along(groups), .simulate_one_sample)
  } else {
    simulated_data <- parallel::mclapply(seq_along(groups), .simulate_one_sample, mc.cores = n_cores, mc.set.seed = TRUE)
  }

  if (verbose) message(">>> Simulation completed successfully.")
  return(simulated_data)
}
