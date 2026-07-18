# =============================================================================
# SaLFC: spatially-aware LFC endpoint (expression) + pilot TREAT threshold.
# =============================================================================

#' @title Spatially-Aware LFC Test with Empirical-Bayes Moderation
#' @description
#' Tests for a differential Target-vs-Reference log-fold-change between two groups
#' of spatial samples, using ground-truth domain labels. Per sample, the statistic
#' \eqn{z = \bar{x}_T - \bar{x}_R} is formed together with its spatial sampling
#' variance \eqn{\tau^2} (derived from an exponential covariance with the per-gene
#' parameters in \code{object@params_expression}). Across samples, the biological
#' variance is estimated, optionally moderated via Smyth-style empirical Bayes
#' (\code{limma::squeezeVar}), and a TREAT-style threshold test against
#' \code{lfc_threshold} is performed with Satterthwaite degrees of freedom. A
#' rank-based (Wilcoxon) alternative is also available.
#'
#' @param processed_data List of samples (e.g. recovered synthetic data). Each
#'   element must contain \code{counts} (genes x spots) and \code{coords} with
#'   columns \code{x}, \code{y}, and the label column selected by
#'   \code{domain_col} (\code{hat_d} and/or \code{domain}); optionally
#'   \code{meta$group} and \code{meta$sample_id}.
#' @param object A spaCraft object, used for per-gene spatial parameters
#'   (\code{theta}) via \code{object@params_expression}.
#' @param genes Character vector of genes to test.
#' @param target_domain Character. Target domain label.
#' @param reference_domain Character or NULL. Reference domain; if NULL, all non-target
#'   spots are used as the reference.
#' @param domain_col Character. Which \code{coords} column defines the target and
#'   reference domains: \code{"hat_d"} (recovered clusters, default) or
#'   \code{"domain"} (ground truth). Falls back to \code{"domain"} when the
#'   requested column is absent (e.g. pilot data, which carries only
#'   ground-truth labels).
#' @param rho_TR Numeric. Assumed Target-Reference biological correlation (< 1).
#' @param adjust_method Multiple-testing correction: "BH" (FDR, default) or
#'   "bonferroni" (FWER).
#' @param test_method Group comparison on per-sample \code{z_hat}: "wilcoxon" or
#'   "t".
#' @param alternative Test direction: "greater" (default), "less", or
#'   "two.sided".
#' @param lfc_threshold Numeric scalar, named numeric vector keyed by gene, or
#'   the string \code{"pilot"}. If a named vector, an optional
#'   \code{"global_tau"} attribute supplies the fallback for missing/non-finite
#'   entries. If \code{"pilot"}, the threshold is estimated on the fly from
#'   \code{object@pilot_data} via \code{estimate_lfc_threshold_from_pilot()} and
#'   then used for the test.
#' @param pilot_threshold_args Named list of extra arguments forwarded to
#'   \code{estimate_lfc_threshold_from_pilot()} when \code{lfc_threshold = "pilot"}
#'   (e.g. \code{list(q = 0.75, min_tau = 0.05)}). \code{object}, \code{genes},
#'   \code{target_domain}, \code{reference_domain}, \code{rho_TR}, \code{adjust_method},
#'   \code{test_method}, and \code{alternative} are supplied automatically but may
#'   be overridden here.
#' @param max_spots_per_domain Integer. Per-domain spot cap; domains larger than
#'   this are randomly down-sampled (bounds the O(n^2) distance matrices). The
#'   global RNG state is restored on exit.
#' @param seed Integer base seed for reproducible down-sampling.
#'
#' @return A data.frame of per-gene moderated statistics and adjusted p-values,
#'   ordered by \code{padj} then \code{p_value}. An empty data.frame is returned
#'   when no gene has sufficient data.
#'
#' @seealso \code{\link{evaluatePowerSaLFC}}, \code{\link{estimate_lfc_threshold_from_pilot}}
#'
#' @importFrom stats var p.adjust pt wilcox.test median setNames
#' @export
#'
#' @examples
#' \dontrun{
#' data("mini_obj_fitted",   package = "spaCraft")
#' data("mini_custom_genes", package = "spaCraft")
#'
#' ## One-point SaLFC power calculation (single design point: K = 5, DE_lfc = 0.3):
#' ## simulate -> recover domains -> rearrange -> SaLFC -> rejection rate.
#' sim <- simulatespaCraft(mini_obj_fitted, n_sample_per_group = 5,
#'   scenario_settings = list(DE_lfc = 0.3), target_domain = "WM",
#'   genes_to_simulate = mini_custom_genes$G_svg,
#'   de_genes = mini_custom_genes$G_spike, seed_base = 1)
#'
#' syn <- pBANSKY(sim, mini_obj_fitted@pilot_data, lambda = 0.5, do_hungarian = TRUE)
#' syn <- setNames(lapply(syn, function(s) list(counts = s$counts, coords = s$coords,
#'          meta = data.frame(sample_id = s$sample_id, group = s$group))),
#'        vapply(syn, function(s) s$sample_id, character(1)))
#' syn <- rearrangeSyntheticToPilot(syn, mini_obj_fitted@pilot_data, match_by = "hat_d")
#'
#' res <- SaLFC(syn, mini_obj_fitted,
#'   genes = unique(c(mini_custom_genes$G_null, mini_custom_genes$G_spike)),
#'   target_domain = "WM", reference_domain = "Layer6", domain_col = "hat_d",
#'   lfc_threshold = "pilot")
#'
#' ## empirical power at this point = fraction of true spike-ins at FDR <= 0.05
#' mean(res$padj[res$gene %in% mini_custom_genes$G_spike] <= 0.05)
#' }
SaLFC <- function(processed_data,
                  object,
                  genes,
                  target_domain = "WM",
                  reference_domain    = NULL,
                  domain_col    = c("hat_d", "domain"),
                  rho_TR        = 0,
                  adjust_method = c("bonferroni", "BH"),
                  test_method   = c("wilcoxon", "t"),
                  alternative   = c("greater", "less", "two.sided"),
                  lfc_threshold = 0.0,
                  pilot_threshold_args = list(),
                  max_spots_per_domain = 5000L,
                  seed = 2026L) {

  adjust_method <- match.arg(adjust_method)
  test_method   <- match.arg(test_method)
  alternative   <- match.arg(alternative)
  domain_col    <- match.arg(domain_col)

  if (!requireNamespace("fields", quietly = TRUE)) stop("Package 'fields' required.")
  if (!is.list(processed_data) || length(processed_data) == 0)
    stop("processed_data must be a non-empty list.")
  if (length(genes) == 0) stop("genes is empty.")

  use_dt <- requireNamespace("data.table", quietly = TRUE)
  max_spots_per_domain <- as.integer(max_spots_per_domain)

  # ---------------------------------------------------------------------------
  # Optional: estimate the TREAT threshold from the pilot data on the fly.
  # Triggered by lfc_threshold = "pilot". The estimator calls this function back
  # with a numeric lfc_threshold (= 0), so there is no infinite recursion.
  # ---------------------------------------------------------------------------
  if (is.character(lfc_threshold) && length(lfc_threshold) == 1L &&
      tolower(lfc_threshold) == "pilot") {
    est_defaults <- list(
      object = object, genes = genes,
      target_domain = target_domain, reference_domain = reference_domain,
      rho_TR = rho_TR, adjust_method = adjust_method,
      test_method = test_method, alternative = alternative,
      max_spots_per_domain = max_spots_per_domain,
      test_fn = sys.function(),
      verbose = FALSE, return_table = FALSE
    )
    est_args <- utils::modifyList(est_defaults, pilot_threshold_args)
    lfc_threshold <- do.call(estimate_lfc_threshold_from_pilot, est_args)
  }

  # ---------------------------------------------------------------------------
  # Helpers
  # ---------------------------------------------------------------------------

  # Resolve a per-gene threshold vector plus a scalar fallback (global_tau).
  .resolve_lfc_threshold <- function(genes, lfc_threshold) {
    if (length(lfc_threshold) == 1L) {
      tau0 <- as.numeric(lfc_threshold)
      if (!is.finite(tau0) || tau0 < 0) tau0 <- 0
      return(list(tau_by_gene = stats::setNames(rep(tau0, length(genes)), genes),
                  global_tau  = tau0))
    }

    if (!is.numeric(lfc_threshold))
      stop("lfc_threshold must be numeric.")
    if (is.null(names(lfc_threshold)))
      stop("When lfc_threshold has length > 1, it must be a named numeric vector keyed by gene.")

    global_tau <- attr(lfc_threshold, "global_tau")
    if (is.null(global_tau) || !is.finite(global_tau) || global_tau < 0) {
      global_tau <- stats::median(as.numeric(lfc_threshold[is.finite(lfc_threshold)]),
                                  na.rm = TRUE)
      if (!is.finite(global_tau) || global_tau < 0) global_tau <- 0
    }

    tau_by_gene <- lfc_threshold[genes]
    tau_by_gene[!is.finite(tau_by_gene)] <- global_tau
    tau_by_gene <- pmax(as.numeric(tau_by_gene), 0)
    names(tau_by_gene) <- genes

    list(tau_by_gene = tau_by_gene, global_tau = global_tau)
  }

  # Per-gene spatial parameters (alpha, rho, nugget) from the spaCraft object.
  .build_theta_table <- function(object, genes) {
    params_list <- object@params_expression
    alpha_vec <- stats::setNames(rep(1.0,  length(genes)), genes)
    rho_vec   <- stats::setNames(rep(10.0, length(genes)), genes)
    sigma_vec <- stats::setNames(rep(0.1,  length(genes)), genes)
    found_vec <- stats::setNames(rep(FALSE, length(genes)), genes)

    if (is.list(params_list)) {
      for (grp in names(params_list)) {
        grp_list <- params_list[[grp]]
        if (is.null(grp_list) || !is.list(grp_list)) next

        g_hit <- intersect(genes[!found_vec], names(grp_list))
        if (length(g_hit) == 0) next

        for (g in g_hit) {
          theta <- grp_list[[g]]$theta
          if (is.null(theta) || !is.numeric(theta)) next

          if (exists("extract_theta_triplet", mode = "function")) {
            th <- extract_theta_triplet(theta)
            alpha_vec[g] <- as.numeric(th["alpha"])
            rho_vec[g]   <- as.numeric(th["rho"])
            sigma_vec[g] <- as.numeric(th["nugget"])
          } else {
            alpha_sq <- 1.0; rho_val <- 10.0; sigma_sq <- 0.1
            if ("alpha" %in% names(theta)) alpha_sq <- as.numeric(theta["alpha"])
            if ("rho"   %in% names(theta)) rho_val  <- as.numeric(theta["rho"])
            if ("phi"   %in% names(theta) && !("rho" %in% names(theta)))
              rho_val <- as.numeric(theta["phi"])
            if ("sigma.sq" %in% names(theta)) sigma_sq <- as.numeric(theta["sigma.sq"])
            alpha_vec[g] <- alpha_sq; rho_vec[g] <- rho_val; sigma_vec[g] <- sigma_sq
          }
          found_vec[g] <- TRUE
        }
        if (all(found_vec)) break
      }
    }

    alpha_vec[!is.finite(alpha_vec) | alpha_vec < 0] <- 1.0
    rho_vec[!is.finite(rho_vec)     | rho_vec  <= 0] <- 10.0
    sigma_vec[!is.finite(sigma_vec) | sigma_vec < 0] <- 0.1

    data.frame(
      gene     = genes,
      alpha    = as.numeric(alpha_vec[genes]),
      rho      = as.numeric(rho_vec[genes]),
      sigma.sq = as.numeric(sigma_vec[genes]),
      found    = as.logical(found_vec[genes]),
      stringsAsFactors = FALSE
    )
  }

  # Spatial sampling variance tau^2 of (mean_T - mean_R) under an exponential GP.
  .tau_sq_from_precomputed_dist <- function(D_TT, D_RR, D_TR, n_T, n_R,
                                            alpha_sq, rho_val, sigma_sq) {
    if (!is.finite(alpha_sq) || alpha_sq < 0) alpha_sq <- 1.0
    if (!is.finite(rho_val)  || rho_val <= 0) rho_val  <- 10.0
    if (!is.finite(sigma_sq) || sigma_sq < 0) sigma_sq <- 0.1

    term_TT <- alpha_sq * sum(exp(-D_TT / rho_val)) + n_T * sigma_sq
    term_RR <- alpha_sq * sum(exp(-D_RR / rho_val)) + n_R * sigma_sq
    term_TR <- alpha_sq * sum(exp(-D_TR / rho_val))

    (term_TT / n_T^2) + (term_RR / n_R^2) - (2 * term_TR / (n_T * n_R))
  }

  tau_info    <- .resolve_lfc_threshold(genes, lfc_threshold)
  tau_by_gene <- tau_info$tau_by_gene
  global_tau  <- tau_info$global_tau
  theta_tbl   <- .build_theta_table(object, genes)
  rownames(theta_tbl) <- theta_tbl$gene

  # CRAN-safe RNG: down-sampling sets seeds internally; restore global state on exit.
  if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
    .old_seed <- get(".Random.seed", envir = globalenv(), inherits = FALSE)
    on.exit(assign(".Random.seed", .old_seed, envir = globalenv()), add = TRUE)
  }

  # ---------------------------------------------------------------------------
  # 1) Per-sample z_hat and tau_sq
  # ---------------------------------------------------------------------------
  nm <- names(processed_data)
  if (is.null(nm) || any(nm == "")) nm <- paste0("s", seq_along(processed_data))

  sample_rows <- vector("list", length(processed_data))
  out_idx <- 0L

  for (i in seq_along(processed_data)) {
    samp   <- processed_data[[i]]
    coords <- samp$coords

    if (is.null(coords) || !all(c("x", "y") %in% names(coords))) next
    if (is.null(samp$counts)) next

    grp <- NA_real_; sid <- nm[i]
    if (!is.null(samp$meta) && "group" %in% names(samp$meta))
      grp <- as.numeric(samp$meta$group)
    if (!is.null(samp$meta) && "sample_id" %in% names(samp$meta))
      sid <- as.character(samp$meta$sample_id)

    # Resolve the label column: requested domain_col if present, else fall back
    # to "domain" (e.g. pilot data carries only ground-truth labels).
    lab_col <- if (domain_col %in% names(coords)) domain_col
    else if ("domain" %in% names(coords)) "domain" else NA_character_
    if (is.na(lab_col)) next
    labs <- as.character(coords[[lab_col]])

    idx_T <- which(labs == target_domain)
    idx_R <- if (is.null(reference_domain)) which(labs != target_domain)
    else which(labs == reference_domain)
    if (length(idx_T) < 2 || length(idx_R) < 2) next

    # Per-domain down-sampling bounds the O(n^2) distance matrices below.
    if (length(idx_T) > max_spots_per_domain) {
      set.seed(seed + i)
      idx_T <- sample(idx_T, max_spots_per_domain)
    }
    if (length(idx_R) > max_spots_per_domain) {
      set.seed(seed + 1000L + i)
      idx_R <- sample(idx_R, max_spots_per_domain)
    }

    g_use <- intersect(genes, rownames(samp$counts))
    if (length(g_use) == 0) next

    counts_mat <- samp$counts[g_use, , drop = FALSE]
    if (!is.matrix(counts_mat)) counts_mat <- as.matrix(counts_mat)

    mean_T <- rowMeans(counts_mat[, idx_T, drop = FALSE])
    mean_R <- rowMeans(counts_mat[, idx_R, drop = FALSE])
    z_hat_vec <- as.numeric(mean_T - mean_R)

    coords_T <- as.matrix(coords[idx_T, c("x", "y"), drop = FALSE])
    coords_R <- as.matrix(coords[idx_R, c("x", "y"), drop = FALSE])

    D_TT <- fields::rdist(coords_T, coords_T)
    D_RR <- fields::rdist(coords_R, coords_R)
    D_TR <- fields::rdist(coords_T, coords_R)

    n_T <- nrow(coords_T); n_R <- nrow(coords_R)
    th_sub <- theta_tbl[g_use, , drop = FALSE]

    tau_sq_vec <- mapply(
      FUN      = .tau_sq_from_precomputed_dist,
      alpha_sq = th_sub$alpha, rho_val = th_sub$rho, sigma_sq = th_sub$sigma.sq,
      MoreArgs = list(D_TT = D_TT, D_RR = D_RR, D_TR = D_TR, n_T = n_T, n_R = n_R),
      SIMPLIFY = TRUE, USE.NAMES = FALSE
    )

    out_idx <- out_idx + 1L
    sample_rows[[out_idx]] <- data.frame(
      gene      = g_use,
      sample_id = rep(sid, length(g_use)),
      group     = rep(grp, length(g_use)),
      z_hat     = z_hat_vec,
      tau_sq    = tau_sq_vec,
      stringsAsFactors = FALSE
    )
  }

  sample_rows <- sample_rows[seq_len(out_idx)]
  if (length(sample_rows) == 0) return(data.frame())

  stats_long <- if (use_dt) {
    as.data.frame(data.table::rbindlist(sample_rows, use.names = TRUE, fill = TRUE),
                  stringsAsFactors = FALSE)
  } else {
    do.call(rbind, sample_rows)
  }
  stats_long <- stats_long[is.finite(stats_long$group) &
                             stats_long$group %in% c(0, 1), , drop = FALSE]
  if (nrow(stats_long) == 0) return(data.frame())

  # ---------------------------------------------------------------------------
  # 2) Wilcoxon (rank-based) test
  # ---------------------------------------------------------------------------
  if (test_method == "wilcoxon") {
    by_gene <- split(stats_long, stats_long$gene)

    out_list <- lapply(by_gene, function(df) {
      if (length(unique(df$group)) < 2 ||
          sum(df$group == 0) < 2 || sum(df$group == 1) < 2) return(NULL)

      g     <- df$gene[1]
      tau_g <- as.numeric(tau_by_gene[[g]])
      if (!is.finite(tau_g) || tau_g < 0) tau_g <- global_tau

      z1 <- df$z_hat[df$group == 1]
      z0 <- df$z_hat[df$group == 0]
      delta_hat <- mean(z1) - mean(z0)
      p_value   <- NA_real_

      if (tau_g > 0) {
        if (alternative == "greater") {
          p_value <- tryCatch(
            stats::wilcox.test(z1, z0, mu = tau_g, alternative = "greater",
                               exact = FALSE)$p.value,
            error = function(e) NA_real_)
        } else if (alternative == "less") {
          p_value <- tryCatch(
            stats::wilcox.test(z1, z0, mu = -tau_g, alternative = "less",
                               exact = FALSE)$p.value,
            error = function(e) NA_real_)
        } else {  # two.sided
          p_one <- tryCatch(
            stats::wilcox.test(
              z1, z0,
              mu          = ifelse(delta_hat >= 0,  tau_g, -tau_g),
              alternative = ifelse(delta_hat >= 0, "greater", "less"),
              exact = FALSE)$p.value,
            error = function(e) NA_real_)
          p_value <- if (is.finite(p_one)) min(1, 2 * p_one) else NA_real_
        }
      } else {
        p_value <- tryCatch(
          stats::wilcox.test(z1, z0, mu = 0, alternative = alternative,
                             exact = FALSE)$p.value,
          error = function(e) NA_real_)
      }

      data.frame(
        gene = g, p_value = p_value, delta_hat = delta_hat,
        t_stat = NA_real_, treat_stat = NA_real_, df = NA_real_,
        v0 = NA_real_, v1 = NA_real_,
        sigma_bio_tilde_g0 = NA_real_, sigma_bio_tilde_g1 = NA_real_,
        df_eb_g0 = NA_real_, df_eb_g1 = NA_real_,
        K0 = sum(df$group == 0), K1 = sum(df$group == 1),
        lfc_threshold = tau_g, stringsAsFactors = FALSE
      )
    })

    out <- do.call(rbind, Filter(Negate(is.null), out_list))
    if (is.null(out) || nrow(out) == 0) return(data.frame())

    okp <- is.finite(out$p_value)
    out$padj <- NA_real_
    out$padj[okp] <- stats::p.adjust(out$p_value[okp], method = adjust_method)
    out$adjust_method <- adjust_method
    out$test_method   <- paste0("wilcoxon_", alternative)
    return(out[order(out$padj, out$p_value), , drop = FALSE])
  }

  # ---------------------------------------------------------------------------
  # 3) Theory-based, EB-moderated t test
  # ---------------------------------------------------------------------------
  denom <- 2 * (1 - rho_TR)
  if (!is.finite(denom) || denom <= 0) stop("rho_TR must be < 1.")

  agg <- by(stats_long, list(stats_long$gene, stats_long$group), function(df) {
    Kc <- nrow(df)
    if (Kc < 2) return(NULL)
    S2     <- stats::var(df$z_hat)
    barTau <- mean(df$tau_sq)
    data.frame(
      gene = df$gene[1], group = df$group[1], K = Kc,
      bar_Z = mean(df$z_hat), S2 = S2, bar_tau_sq = barTau,
      sigma_bio_hat = max(0, (S2 - barTau) / denom),
      stringsAsFactors = FALSE
    )
  })

  agg_df <- do.call(rbind, Filter(Negate(is.null), as.list(agg)))
  if (is.null(agg_df) || nrow(agg_df) == 0) return(data.frame())

  tab        <- table(agg_df$gene)
  keep_genes <- names(tab)[tab >= 2]
  agg_df     <- agg_df[agg_df$gene %in% keep_genes, , drop = FALSE]
  if (nrow(agg_df) == 0) return(data.frame())

  # Smyth-style empirical-Bayes moderation of the biological variance (per group).
  use_limma <- requireNamespace("limma", quietly = TRUE)
  agg_df$sigma_bio_tilde <- agg_df$sigma_bio_hat
  agg_df$df_eb           <- agg_df$K - 1

  if (use_limma) {
    for (grp in c(0, 1)) {
      sub <- agg_df[agg_df$group == grp, , drop = FALSE]
      if (nrow(sub) < 3) next
      sv <- sub$sigma_bio_hat
      df <- sub$K - 1
      if (all(!is.finite(sv)) || all(sv <= 0)) next
      sq <- limma::squeezeVar(var = sv, df = df)
      agg_df$sigma_bio_tilde[agg_df$group == grp] <- sq$var.post
      agg_df$df_eb[agg_df$group == grp]           <- sq$df.prior + (sub$K - 1)
    }
  }

  split0 <- agg_df[agg_df$group == 0, , drop = FALSE]
  split1 <- agg_df[agg_df$group == 1, , drop = FALSE]
  split0 <- split0[match(keep_genes, split0$gene), , drop = FALSE]
  split1 <- split1[match(keep_genes, split1$gene), , drop = FALSE]

  ok <- is.finite(split0$bar_Z) & is.finite(split1$bar_Z) &
    is.finite(split0$sigma_bio_tilde) & is.finite(split1$sigma_bio_tilde) &
    is.finite(split0$bar_tau_sq) & is.finite(split1$bar_tau_sq)
  split0 <- split0[ok, , drop = FALSE]
  split1 <- split1[ok, , drop = FALSE]
  if (nrow(split0) == 0) return(data.frame())

  v0        <- (denom * split0$sigma_bio_tilde + split0$K * split0$bar_tau_sq) / (split0$K^2)
  v1        <- (denom * split1$sigma_bio_tilde + split1$K * split1$bar_tau_sq) / (split1$K^2)

  v0        <- pmax(v0, 1e-5)
  v1        <- pmax(v1, 1e-5)

  delta_hat <- split1$bar_Z - split0$bar_Z
  se        <- sqrt(v0 + v1)

  nu0     <- split0$df_eb
  nu1     <- split1$df_eb
  df_satt <- (v0 + v1)^2 / (v0^2 / nu0 + v1^2 / nu1)

  genes_final <- split0$gene
  tau_final   <- tau_by_gene[genes_final]
  tau_final[!is.finite(tau_final)] <- global_tau
  tau_final <- pmax(as.numeric(tau_final), 0)

  ok_num <- is.finite(delta_hat) & is.finite(se) & se > 0 &
    is.finite(df_satt) & df_satt > 0
  t_stat     <- rep(NA_real_, length(delta_hat))
  treat_stat <- rep(NA_real_, length(delta_hat))
  p_value    <- rep(NA_real_, length(delta_hat))
  t_stat[ok_num] <- delta_hat[ok_num] / se[ok_num]

  if (alternative == "greater") {
    treat_stat[ok_num] <- (delta_hat[ok_num] - tau_final[ok_num]) / se[ok_num]
    p_value[ok_num]    <- stats::pt(treat_stat[ok_num], df = df_satt[ok_num], lower.tail = FALSE)
  } else if (alternative == "less") {
    treat_stat[ok_num] <- (-delta_hat[ok_num] - tau_final[ok_num]) / se[ok_num]
    p_value[ok_num]    <- stats::pt(treat_stat[ok_num], df = df_satt[ok_num], lower.tail = FALSE)
  } else {  # two.sided
    treat_stat[ok_num] <- (abs(delta_hat[ok_num]) - tau_final[ok_num]) / se[ok_num]
    p_value[ok_num]    <- pmin(1, 2 * stats::pt(treat_stat[ok_num], df = df_satt[ok_num],
                                                lower.tail = FALSE))
  }

  padj <- rep(NA_real_, length(p_value))
  okp  <- is.finite(p_value)
  padj[okp] <- stats::p.adjust(p_value[okp], method = adjust_method)

  out <- data.frame(
    gene = genes_final, p_value = p_value, padj = padj,
    adjust_method = adjust_method, test_method = paste0("t_", alternative),
    lfc_threshold = tau_final, delta_hat = delta_hat,
    t_stat = t_stat, treat_stat = treat_stat, df = df_satt,
    v0 = v0, v1 = v1,
    sigma_bio_tilde_g0 = split0$sigma_bio_tilde, sigma_bio_tilde_g1 = split1$sigma_bio_tilde,
    df_eb_g0 = nu0, df_eb_g1 = nu1, K0 = split0$K, K1 = split1$K,
    stringsAsFactors = FALSE
  )

  out[order(out$padj, out$p_value), , drop = FALSE]
}


#' @title Estimate Gene-Specific SaLFC TREAT Thresholds from Pilot Data
#' @description
#' Fast plug-in estimator (no bootstrap) of gene-specific LFC thresholds
#' \eqn{\tau_g} for the TREAT-style SaLFC test, using only pilot data. SaLFC is
#' run once on the pilot samples with \code{lfc_threshold = 0}; \eqn{\tau_g} is
#' then defined from the pilot effect size and its uncertainty in a
#' direction-aware way:
#' \deqn{\tau_g = \max(\mathrm{effect\_ref}_g,\; z_q \cdot \mathrm{se}_g),}
#' where \eqn{\mathrm{effect\_ref}} depends on \code{alternative} and
#' \eqn{z_q = \Phi^{-1}(q)}.
#'
#' This module is normally invoked automatically by
#' \code{SaLFC(..., lfc_threshold = "pilot")}, but it
#' can also be called directly and its result passed as \code{lfc_threshold}.
#'
#' @param object A spaCraft object with \code{@pilot_data}.
#' @param genes Character vector of (null) genes, e.g. the stable/null set.
#' @param target_domain Character. Target domain label.
#' @param reference_domain Character. Reference domain label.
#' @param rho_TR Numeric. Assumed Target-Reference biological correlation (< 1).
#' @param adjust_method Multiple-testing method passed through to SaLFC.
#' @param test_method SaLFC group comparison: "t" or "wilcoxon".
#' @param alternative Test direction: "greater" (default), "less", "two.sided".
#' @param q Numeric in (0,1). Standard-normal quantile level for \eqn{z_q}.
#' @param min_tau Numeric. Lower bound for \eqn{\tau_g}.
#' @param cap_effect Numeric scalar, named numeric vector, or NULL. Optional
#'   upper cap on \eqn{\tau_g} (per gene if named).
#' @param max_spots_per_domain Integer. Per-sample spot cap before running SaLFC
#'   on the pilot (bounds memory). The global RNG state is restored on exit.
#' @param verbose Logical.
#' @param return_table Logical. If TRUE return a data.frame, else a named vector.
#' @param test_fn Function. The SaLFC test to run on the pilot. If NULL, resolves
#'   \code{SaLFC} by name. When invoked via
#'   \code{SaLFC(..., lfc_threshold = "pilot")}, the
#'   running test function is passed automatically, so renaming it does not break
#'   this module.
#' @param B,shrink_to_global,seed Deprecated and ignored; retained for backward
#'   compatibility.
#'
#' @return A named numeric vector \eqn{\tau_g} carrying an \code{attr(., "global_tau")},
#'   or a data.frame when \code{return_table = TRUE}.
#'
#' @seealso \code{\link{SaLFC}}
#'
#' @importFrom stats qnorm median setNames
#' @export
estimate_lfc_threshold_from_pilot <- function(object,
                                              genes,
                                              target_domain = "WM",
                                              reference_domain    = "Layer6",
                                              rho_TR        = 0,
                                              adjust_method = "bonferroni",
                                              test_method   = "t",
                                              alternative   = c("greater", "less", "two.sided"),
                                              q             = 0.75,
                                              min_tau       = 0,
                                              cap_effect    = NULL,
                                              max_spots_per_domain = 5000L,
                                              verbose       = TRUE,
                                              return_table  = FALSE,
                                              test_fn       = NULL,
                                              B = 200, shrink_to_global = 0.25, seed = 1) {
  alternative <- match.arg(alternative)

  if (is.null(object@pilot_data) || length(object@pilot_data) == 0)
    stop("object@pilot_data is empty.")
  if (length(genes) == 0) stop("genes is empty.")
  if (!is.finite(q) || q <= 0 || q >= 1) stop("q must be in (0,1).")
  if (!is.finite(min_tau) || min_tau < 0) stop("min_tau must be >= 0.")

  # The SaLFC test function to run on the pilot. Decoupled from its bound name so
  # renaming the test function (e.g. to SaLFC) does not break this module.
  if (is.null(test_fn)) {
    if (exists("SaLFC", mode = "function")) {
      test_fn <- match.fun("SaLFC")
    } else if (exists("test_saLFC_theory_ground_truth_EB", mode = "function")) {
      test_fn <- match.fun("test_saLFC_theory_ground_truth_EB")
    } else {
      stop("No SaLFC test function supplied; pass it via `test_fn`.")
    }
  }
  if (!is.function(test_fn)) stop("`test_fn` must be a function.")

  max_spots_per_domain <- as.integer(max_spots_per_domain)
  if (verbose)
    message(">>> [pilot tau_g] fast plug-in mode: B / shrink_to_global / seed are ignored.")

  # CRAN-safe RNG: down-sampling sets seeds internally; restore global state on exit.
  if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
    .old_seed <- get(".Random.seed", envir = globalenv(), inherits = FALSE)
    on.exit(assign(".Random.seed", .old_seed, envir = globalenv()), add = TRUE)
  }

  # ---------------------------------------------------------------------------
  # 1) Build processed pilot data (with per-sample down-sampling)
  # ---------------------------------------------------------------------------
  pilot_pd <- lapply(seq_along(object@pilot_data), function(i) {
    s <- object@pilot_data[[i]]
    if (is.null(s$coords) || is.null(s$counts) || is.null(s$group)) return(NULL)
    coords <- s$coords
    if (!("domain" %in% names(coords))) return(NULL)

    N_spots <- nrow(coords)
    if (N_spots > max_spots_per_domain) {
      set.seed(2026 + i)
      keep_idx <- sample(seq_len(N_spots), max_spots_per_domain)
    } else {
      keep_idx <- seq_len(N_spots)
    }

    sid <- if (!is.null(names(object@pilot_data)) && names(object@pilot_data)[i] != "") {
      names(object@pilot_data)[i]
    } else {
      paste0("pilot_", i)
    }

    list(
      counts    = s$counts[, keep_idx, drop = FALSE],
      logcounts = s$logcounts[, keep_idx, drop = FALSE],
      coords    = coords[keep_idx, , drop = FALSE],
      meta      = data.frame(sample_id = sid, group = as.numeric(s$group),
                             stringsAsFactors = FALSE)
    )
  })

  pilot_pd <- Filter(Negate(is.null), pilot_pd)
  if (length(pilot_pd) == 0) stop("No valid pilot samples with counts/coords/group.")
  names(pilot_pd) <- vapply(pilot_pd, function(x) x$meta$sample_id, character(1))

  # ---------------------------------------------------------------------------
  # 2) Run SaLFC once on the pilot with threshold = 0 (base case: numeric tau)
  # ---------------------------------------------------------------------------
  # Pilot data is ground truth and carries only the "domain" column (no hat_d),
  # so force domain_col = "domain" here. Passed only if test_fn accepts it, so a
  # legacy test function without that argument still works.
  .pilot_args <- list(
    processed_data       = pilot_pd,
    object               = object,
    genes                = genes,
    target_domain        = target_domain,
    reference_domain           = reference_domain,
    rho_TR               = rho_TR,
    adjust_method        = adjust_method,
    test_method          = test_method,
    alternative          = alternative,
    lfc_threshold        = 0.0,
    max_spots_per_domain = max_spots_per_domain
  )
  if ("domain_col" %in% names(formals(test_fn))) .pilot_args$domain_col <- "domain"

  pilot_res <- tryCatch(
    do.call(test_fn, .pilot_args),
    error = function(e) {
      message("\n[ERROR in SaLFC] ", e$message)
      NULL
    }
  )

  if (is.null(pilot_res) || !is.data.frame(pilot_res) || nrow(pilot_res) == 0) {
    if (verbose) message(">>> [pilot tau_g] pilot result empty -> returning zeros")
    tau0 <- stats::setNames(rep(0, length(genes)), genes)
    attr(tau0, "global_tau") <- 0
    if (return_table) {
      out_df <- data.frame(
        gene = genes, tau_lfc = 0, effect_ref = NA_real_, se = NA_real_,
        z_q = stats::qnorm(q), alternative = alternative, stringsAsFactors = FALSE
      )
      attr(out_df, "global_tau") <- 0
      return(out_df)
    }
    return(tau0)
  }

  # ---------------------------------------------------------------------------
  # 3) Direction-aware fast plug-in tau_g
  # ---------------------------------------------------------------------------
  idx       <- match(genes, pilot_res$gene)
  delta_hat <- stats::setNames(rep(NA_real_, length(genes)), genes)
  v0        <- stats::setNames(rep(NA_real_, length(genes)), genes)
  v1        <- stats::setNames(rep(NA_real_, length(genes)), genes)

  ok <- which(!is.na(idx))
  delta_hat[ok] <- pilot_res$delta_hat[idx[ok]]
  v0[ok]        <- pilot_res$v0[idx[ok]]
  v1[ok]        <- pilot_res$v1[idx[ok]]

  z_q <- stats::qnorm(q)
  if (!is.finite(z_q)) z_q <- 0
  z_q <- max(z_q, 0)

  se <- sqrt(pmax(0, v0 + v1))

  effect_ref <- switch(
    alternative,
    greater   = pmax(delta_hat, 0),
    less      = pmax(-delta_hat, 0),
    two.sided = abs(delta_hat)
  )
  effect_ref[!is.finite(effect_ref)] <- 0

  tau_raw <- pmax(effect_ref, z_q * se)
  tau_raw[!is.finite(tau_raw)] <- effect_ref[!is.finite(tau_raw)]
  tau_raw[!is.finite(tau_raw)] <- (z_q * se)[!is.finite(tau_raw)]
  tau_raw[!is.finite(tau_raw)] <- 0
  tau_raw <- pmax(tau_raw, min_tau)

  if (!is.null(cap_effect)) {
    if (length(cap_effect) == 1L) {
      if (is.finite(cap_effect) && cap_effect > 0) tau_raw <- pmin(tau_raw, cap_effect)
    } else {
      if (is.null(names(cap_effect)))
        stop("If cap_effect has length > 1, it must be a named numeric vector.")
      cap_vec <- cap_effect[genes]
      ok_cap  <- is.finite(cap_vec) & cap_vec > 0
      tau_raw[ok_cap] <- pmin(tau_raw[ok_cap], cap_vec[ok_cap])
    }
  }

  tau_raw[!is.finite(tau_raw)] <- 0
  names(tau_raw) <- genes

  global_tau <- stats::median(tau_raw, na.rm = TRUE)
  if (!is.finite(global_tau) || global_tau < 0) global_tau <- 0

  if (verbose)
    message(sprintf(
      ">>> [pilot tau_g] done. alternative=%s, q=%.2f (z_q=%.3f), global_tau=%.4f, median(tau_g)=%.4f",
      alternative, q, z_q, global_tau, stats::median(tau_raw, na.rm = TRUE)))

  if (return_table) {
    out_df <- data.frame(
      gene = genes, tau_lfc = as.numeric(tau_raw[genes]),
      effect_ref = as.numeric(effect_ref[genes]), se = as.numeric(se[genes]),
      z_q = rep(z_q, length(genes)), alternative = rep(alternative, length(genes)),
      stringsAsFactors = FALSE
    )
    attr(out_df, "global_tau") <- global_tau
    return(out_df)
  }

  attr(tau_raw, "global_tau") <- global_tau
  tau_raw
}
