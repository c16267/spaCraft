#' @title Stability- and Validity-Driven Hyperparameter Tuning
#' @description
#' Tunes the spatial hyperparameters of a \code{spaCraft} object so that the
#' resulting simulation behaves as a \emph{faithful and stable power oracle}.
#' The recovery weight \code{lambda_p} (pBANKSY) is tuned once by maximizing
#' clustering reproducibility (ARI), exactly as before. The conditional-texture
#' weight \code{lambda_cond} is tuned with a regularized-noncentrality objective
#' that fixes two defects of a pure \eqn{Var(T)}-minimization:
#'
#' \enumerate{
#'   \item \strong{Objective.} Power is governed by the noncentrality
#'     \eqn{\delta = E[T]/sd[T]}, not by \eqn{sd[T]} alone. Minimizing variance
#'     can be achieved by attenuating the signal, which collapses power. The
#'     objective is therefore a regularized noncentrality (equivalently the
#'     reciprocal of a floored coefficient of variation),
#'     \deqn{\widehat\delta(\lambda) = \bar T / (sd(T) + \varepsilon\,s_0),}
#'     which is scale-free and cannot be gamed by uniform shrinkage. Both
#'     endpoints are placed on the \emph{same} standardized (t-statistic)
#'     footing: the SaLFC endpoint uses \code{mean(abs(t_stat))} and the LOR
#'     endpoint uses the standardized statistic \code{LOR()$stat} (not the raw
#'     \code{delta_hat}).
#'   \item \strong{Effect ordering.} The default
#'     \code{lambda_cond_selection = "per_effect"} selects the pointwise-optimal
#'     weight at every effect size,
#'     \eqn{\lambda^*(e)=\arg\max_\lambda \widehat\delta(e,\lambda)\,w(\lambda)},
#'     where \eqn{w(\lambda)} is the null-calibration weight below. Because the
#'     generator's noncentrality is monotone in effect size for each fixed
#'     \eqn{\lambda}, the pointwise maximum \eqn{\max_\lambda\delta(e,\lambda)} is
#'     itself monotone in \eqn{e} in expectation, so per-effect selection does not
#'     sacrifice the effect ordering; residual non-monotonicity is Monte-Carlo
#'     noise and is removed by isotonizing (PAVA) the \emph{reported}
#'     noncentrality (\code{SNR_iso}). The per-effect \code{lambda_cond} is stored
#'     in \code{by_effect} and consumed downstream (e.g.
#'     \code{\link{evaluatePowerSaLFC}}, \code{\link{evaluatePowerLOR}}), which
#'     already match \code{lambda_cond} to each effect size. Setting
#'     \code{lambda_cond_selection = "shared"} recovers the previous behaviour: a
#'     single \code{lambda_cond} is selected jointly across the effect grid by
#'     maximizing the worst-case (minimum-over-effects) calibrated noncentrality,
#'     trading a little per-effect optimality for lower selection variance.
#' }
#'
#' A null-calibration guard (\code{calibration_guard}) further down-weights any
#' \code{lambda_cond} whose empirical type-I error under the matched null
#' (SaLFC: \code{DE_lfc = 0}, no DE genes; LOR: \code{target_prop_case = p0}, the
#' pilot baseline proportion) is inflated, so noncentrality is never bought by
#' biasing the generator.
#'
#' @param object A \code{spaCraft} object.
#' @param target_domain Character. Target domain name.
#' @param reference_domain Character. Reference domain name.
#' @param G_DE Character vector. Genes with a true DE effect (H1 scenario).
#' @param scenario_H1 List. Joint scenario settings for H1. \code{DE_lfc} and
#'   \code{target_prop_case} may each be a scalar or a numeric vector; when a
#'   vector, \code{lambda_cond} is evaluated at every value and (depending on
#'   \code{lambda_cond_selection}) either a per-effect or a single shared value is
#'   selected across the grid.
#' @param lambda_p_grid Numeric vector. Grid of \code{lambda_p} to test.
#' @param lambda_c_grid Numeric vector. Grid of \code{lambda_cond} to test.
#' @param lambda_cond_selection Character. How the conditional-texture weight is
#'   chosen across the effect grid. \code{"per_effect"} (default) returns the
#'   pointwise-optimal (calibrated) \code{lambda_cond} at each effect size,
#'   stored in \code{by_effect$lambda_cond} and matched per effect by the
#'   downstream power functions; \code{"shared"} returns a single value maximizing
#'   the worst-case calibrated noncentrality across effects (the previous
#'   behaviour).
#' @param B Integer. Number of replicates per grid point (Default: 10).
#' @param K_syn Integer. Number of samples per group during tuning (Default: 10).
#' @param n_cores Integer. Number of cores for simulation.
#' @param seed Integer. Random seed for reproducible tuning. Default is 2026.
#' @param verbose Logical. Prints progress.
#' @param eps_floor Numeric. Denominator floor multiplier \eqn{\varepsilon} for
#'   the regularized noncentrality (Default: 1e-3). Stabilizes the \eqn{E[T]\to 0}
#'   regime (relevant for LOR near the null).
#' @param calibration_guard Logical. If \code{TRUE} (default) compute the
#'   null-calibration weight and down-weight inflated \code{lambda_cond}.
#' @param alpha Numeric. Nominal significance level for the calibration guard.
#' @param calib_mult Numeric > 1. Tolerance multiplier: the calibration weight is
#'   \eqn{\approx 1} while the empirical type-I error stays below
#'   \code{alpha * calib_mult} and decays smoothly beyond it.
#'
#' @return An object of class \code{spaCraftTuning}. \code{optimal_params} has
#'   one entry per endpoint (\code{SaLFC}, \code{LOR}); each contains the shared
#'   \code{lambda_p}, a representative scalar \code{lambda_cond} (the value at the
#'   median effect), a \code{by_effect} data frame (the selected \code{lambda_cond}
#'   \emph{per effect size}, realized noncentrality \code{SNR}, its isotonic
#'   projection \code{SNR_iso}, and the per-effect calibration weight
#'   \code{calib_weight}), and \code{shared_lambda_cond} (the single shared weight
#'   when \code{lambda_cond_selection = "shared"}, otherwise \code{NA}).
#'   \code{summary_stage2} holds the per-endpoint noncentrality surfaces (with
#'   \code{mean_T}, \code{sd_T}, \code{SNR}, and \code{Var_T} for reference).
#'   \code{diagnostics} stores the calibration table.
#' @importFrom mclust adjustedRandIndex
#' @importFrom stats median sd isoreg plogis
#' @export
#'
tuneSpaDesignLambdas <- function(object,
                                 target_domain, reference_domain, G_DE, scenario_H1,
                                 lambda_p_grid = seq(0.2, 0.8, 0.1),
                                 lambda_c_grid = seq(0.2, 0.8, 0.1),
                                 lambda_cond_selection = c("per_effect", "shared"),
                                 B = 10, K_syn = 10, n_cores = 1,
                                 seed = 2026,
                                 verbose = TRUE,
                                 eps_floor = 1e-3,
                                 calibration_guard = TRUE, alpha = 0.05,
                                 calib_mult = 2.0) {

  if (!requireNamespace("mclust", quietly = TRUE)) stop("Package 'mclust' is required.")
  lambda_cond_selection <- match.arg(lambda_cond_selection)

  # CRAN-safe RNG state management
  if (!is.null(seed)) {
    if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
      .old_seed <- get(".Random.seed", envir = globalenv(), inherits = FALSE)
      on.exit(assign(".Random.seed", .old_seed, envir = globalenv()), add = TRUE)
    }
    set.seed(seed)
  }

  # ---- Effect-size grids (each may be a scalar or a vector) ----
  de_grid  <- scenario_H1$DE_lfc            # SaLFC effect grid
  tpc_grid <- scenario_H1$target_prop_case  # LOR effect grid
  if (is.null(de_grid)) de_grid <- 0.0
  has_lor  <- !is.null(tpc_grid) && length(tpc_grid) > 0 && any(is.finite(tpc_grid))
  de_grid  <- sort(unique(as.numeric(de_grid)))
  if (has_lor) tpc_grid <- sort(unique(as.numeric(tpc_grid)))
  tpc_rep  <- if (has_lor) stats::median(tpc_grid) else NULL

  # Base scenario with the swept effect dimensions removed; set inside the loops.
  scen_base <- scenario_H1
  scen_base$DE_lfc <- NULL
  scen_base$target_prop_case <- NULL

  # Pilot baseline proportion p0 (LOR null: target_prop_case == p0  =>  delta == 0).
  lor_null_prop <- tryCatch({
    b0 <- object@params_composition$beta_binomial$beta0
    if (is.null(b0) || !is.finite(as.numeric(b0))) NULL else stats::plogis(as.numeric(b0))
  }, error = function(e) NULL)

  if (verbose) {
    message(">>> Joint H1 scenario split for orthogonal tuning:")
    message(sprintf("    - SaLFC: DE_lfc in {%s}; target_prop_case = NULL",
                    paste(sprintf("%.2f", de_grid), collapse = ", ")))
    if (has_lor) {
      message(sprintf("    - LOR  : DE_lfc = 0.0; target_prop_case in {%s}",
                      paste(sprintf("%.2f", tpc_grid), collapse = ", ")))
    } else {
      message("    - LOR  : skipped (no target_prop_case supplied)")
    }
    message(sprintf("    objective = snr | selection = %s | calibration = %s",
                    lambda_cond_selection,
                    if (calibration_guard) sprintf("on (alpha=%.3f)", alpha) else "off"))
  }

  # ===========================================================================
  # Auxiliary functions
  # ===========================================================================
  build_processed <- function(sim_data, lp) {
    syn_d_hat <- suppressMessages(pBANSKY(
      sim_data = sim_data, pilot_data = object@pilot_data, lambda = lp, k_neighbors = 15,
      verbose = FALSE, scale_pilot_coords = FALSE, scale_sim_coords = FALSE,
      hard_guidance = FALSE, n_pcs = 20, do_hungarian = TRUE, use_G = TRUE,
      balance_mode = "energy_proportional"
    ))
    hat_d_vec <- unlist(lapply(syn_d_hat, function(x) as.character(x$coords$hat_d)))
    synthetic_data <- stats::setNames(lapply(syn_d_hat, function(item) {
      list(counts = item$counts, logcounts = item$logcounts, coords = item$coords,
           meta = data.frame(sample_id = item$sample_id, group = item$group,
                             stringsAsFactors = FALSE))
    }), vapply(syn_d_hat, function(x) as.character(x$sample_id), character(1)))
    list(syn_data = synthetic_data, hat_d = hat_d_vec)
  }

  # Simulate once at the given scenario (K = K_syn) and return BOTH the
  # standardized test statistic and the raw p-value(s) for the requested endpoint.
  simulate_stat <- function(scen, de_genes, endpoint, lp) {
    sim <- suppressMessages(simulatespaCraft(
      object, n_sample_per_group = K_syn, scenario_settings = scen,
      target_domain = target_domain, genes_to_simulate = G_DE, de_genes = de_genes,
      spatial_mode = "original", graph_k = 15L, verbose = FALSE, n_cores = n_cores
    ))
    pd <- build_processed(sim, lp)$syn_data
    if (endpoint == "SaLFC") {
      res <- suppressMessages(SaLFC(
        processed_data = pd, object = object, genes = G_DE,
        target_domain = target_domain, reference_domain = reference_domain,
        domain_col = "hat_d", adjust_method = "bonferroni",
        test_method = "t", lfc_threshold = 0
      ))
      stat <- if (is.data.frame(res) && "t_stat" %in% names(res))
        mean(abs(res$t_stat), na.rm = TRUE) else NA_real_
      pval <- if (is.data.frame(res) && "p_value" %in% names(res))
        res$p_value else NA_real_
    } else {
      res <- suppressMessages(LOR(
        processed_data = pd, object = object, target_domain = target_domain,
        reference_domain = reference_domain, test_method = "t", null_type = "treat",
        tau = 0, eps = 0.01, verbose = FALSE
      ))
      # NOTE: the standardized statistic (delta_tilde / se), NOT the raw effect.
      stat <- if (is.list(res) && !is.null(res$stat) && is.finite(res$stat))
        as.numeric(res$stat) else NA_real_
      pval <- if (is.list(res) && !is.null(res$p_value)) as.numeric(res$p_value) else NA_real_
    }
    list(stat = stat, pval = pval)
  }

  # Set the effect dimension for an endpoint on a fresh copy of scen_base.
  set_effect <- function(value, endpoint) {
    scen <- scen_base
    if (endpoint == "SaLFC") { scen$DE_lfc <- value; scen$target_prop_case <- NULL }
    else                     { scen$DE_lfc <- 0.0;   scen$target_prop_case <- value }
    scen
  }

  # Regularized noncentrality on a vector of replicate statistics.
  reg_summary <- function(Tvals) {
    Tvals <- Tvals[is.finite(Tvals)]
    n <- length(Tvals)
    m <- if (n >= 1L) mean(Tvals) else NA_real_
    s <- if (n >= 2L) stats::sd(Tvals) else NA_real_
    c(mean_T = m, sd_T = s, n_valid = n)
  }

  # ===========================================================================
  # Stage 1: tune lambda_p ONCE (recovery weight), at the median effect
  # ===========================================================================
  lambda_c_0 <- stats::median(lambda_c_grid)
  if (verbose) message(sprintf("\n>>> [Stage 1] Tuning lambda_p (lambda_cond fixed at %.2f)", lambda_c_0))
  scen_s1 <- scen_base
  scen_s1$DE_lfc <- 0.0
  scen_s1$target_prop_case <- tpc_rep      # NULL when there is no LOR effect
  scen_s1$lambda_cond <- lambda_c_0

  stage1_summary <- data.frame()
  for (lp in lambda_p_grid) {
    if (verbose) cat(sprintf("  -> lambda_p = %.2f ... ", lp))
    hat_d_list <- vector("list", B)
    for (b in seq_len(B)) {
      sim_data <- suppressMessages(simulatespaCraft(
        object, n_sample_per_group = K_syn, scenario_settings = scen_s1,
        target_domain = target_domain, genes_to_simulate = G_DE, de_genes = character(0),
        spatial_mode = "original", graph_k = 15L, verbose = FALSE, n_cores = n_cores
      ))
      hat_d_list[[b]] <- build_processed(sim_data, lp)$hat_d
    }
    ari_vals <- c()
    for (b1 in seq_len(B - 1)) for (b2 in (b1 + 1):B)
      ari_vals <- c(ari_vals, mclust::adjustedRandIndex(hat_d_list[[b1]], hat_d_list[[b2]]))
    mean_ari <- mean(ari_vals, na.rm = TRUE)
    stage1_summary <- rbind(stage1_summary, data.frame(lambda_p = lp, ARI = mean_ari))
    if (verbose) cat(sprintf("ARI: %.3f\n", mean_ari))
  }
  best_lp <- stage1_summary$lambda_p[which.max(stage1_summary$ARI)]  # STRICT MAX (ARI)

  # ===========================================================================
  # Stage 2 engine (shared by SaLFC and LOR)
  # ===========================================================================
  tune_endpoint <- function(endpoint, effect_grid, effect_label,
                            de_genes_h1, null_prop) {

    # --- (a) Noncentrality surface over effect x lambda at K_syn ---------------
    curves <- data.frame()
    for (eff in effect_grid) {
      if (verbose) cat(sprintf("  -> %s = %.4g\n", effect_label, eff))
      for (lc in lambda_c_grid) {
        scen <- set_effect(eff, endpoint); scen$lambda_cond <- lc
        Tvals <- numeric(B)
        for (b in seq_len(B))
          Tvals[b] <- simulate_stat(scen, de_genes_h1, endpoint, best_lp)$stat
        rs <- reg_summary(Tvals)
        curves <- rbind(curves, data.frame(
          effect = eff, lambda_cond = lc,
          mean_T = rs[["mean_T"]], sd_T = rs[["sd_T"]], n_valid = rs[["n_valid"]]
        ))
      }
    }
    # pooled (robust) scale and regularized noncentrality
    s0 <- stats::median(curves$sd_T[is.finite(curves$sd_T)], na.rm = TRUE)
    if (!is.finite(s0) || s0 <= 0) s0 <- 1
    curves$SNR   <- curves$mean_T / (curves$sd_T + eps_floor * s0)
    curves$Var_T <- curves$sd_T^2

    effs <- sort(unique(curves$effect))
    nE <- length(effs); nL <- length(lambda_c_grid)
    cell <- function(col) {
      M <- matrix(NA_real_, nE, nL,
                  dimnames = list(as.character(effs), as.character(lambda_c_grid)))
      for (i in seq_len(nE)) for (j in seq_len(nL)) {
        v <- curves[[col]][curves$effect == effs[i] &
                             abs(curves$lambda_cond - lambda_c_grid[j]) < 1e-12]
        if (length(v) == 1L) M[i, j] <- v
      }
      M
    }
    SNR_mat <- cell("SNR")

    # --- (b) Null-calibration weight (per lambda) ------------------------------
    calib_tbl <- data.frame()
    cal <- stats::setNames(rep(1, nL), as.character(lambda_c_grid))
    do_calib <- calibration_guard && !(endpoint == "LOR" && is.null(null_prop))
    if (calibration_guard && endpoint == "LOR" && is.null(null_prop) && verbose)
      message("     [calibration] LOR null proportion p0 unavailable -> guard skipped for LOR.")
    if (do_calib) {
      if (verbose) message(sprintf("     [calibration] %s null type-I scan (alpha=%.3f)",
                                   endpoint, alpha))
      for (lc in lambda_c_grid) {
        scen <- scen_base
        if (endpoint == "SaLFC") { scen$DE_lfc <- 0.0; scen$target_prop_case <- NULL }
        else                     { scen$DE_lfc <- 0.0; scen$target_prop_case <- null_prop }
        scen$lambda_cond <- lc
        rej <- numeric(0)
        for (b in seq_len(B)) {
          pv <- simulate_stat(scen, character(0), endpoint, best_lp)$pval
          pv <- pv[is.finite(pv)]
          if (length(pv) > 0) rej <- c(rej, mean(pv < alpha))
        }
        typeI <- if (length(rej) > 0) mean(rej) else NA_real_
        if (is.finite(typeI)) {
          excess <- max(0, (typeI - alpha) / (alpha * (calib_mult - 1) + 1e-8))
          w <- exp(-0.5 * excess^2)            # ~1 if calibrated, ->0 if inflated
        } else {
          w <- NA_real_
        }
        cal[as.character(lc)] <- w
        calib_tbl <- rbind(calib_tbl, data.frame(lambda_cond = lc, typeI = typeI, weight = w))
      }
    }

    # --- (c) Select lambda_cond (per-effect by default, or a single shared) ----
    #   Validity-weighted SNR objective: each column (lambda) is down-weighted by
    #   its null-calibration weight, so noncentrality can never be purchased by
    #   type-I inflation. "per_effect" takes the pointwise argmax at every effect;
    #   "shared" takes the single lambda maximizing the worst-case (min-over-
    #   effects) objective.
    w_lambda <- cal
    w_lambda[!is.finite(w_lambda)] <- 0
    OBJ_adj <- sweep(SNR_mat, 2, w_lambda, `*`)   # nE x nL calibrated objective

    if (identical(lambda_cond_selection, "shared")) {
      agg <- apply(OBJ_adj, 2, function(col) {
        col <- col[is.finite(col)]; if (length(col) == 0) -Inf else min(col)  # worst-case
      })
      g_star <- which.max(agg)
      g_idx  <- rep(g_star, nE)                    # same lambda at every effect
    } else {
      g_idx <- vapply(seq_len(nE), function(i) {
        row <- OBJ_adj[i, ]
        if (all(!is.finite(row))) {                # degenerate: fall back to raw SNR
          row <- SNR_mat[i, ]
          if (all(!is.finite(row))) return(as.integer(ceiling(nL / 2)))  # then mid-grid
        }
        as.integer(which.max(replace(row, !is.finite(row), -Inf)))
      }, integer(1))
    }

    chosen_lc    <- lambda_c_grid[g_idx]                                         # per effect
    realized_snr <- vapply(seq_len(nE), function(i) SNR_mat[i, g_idx[i]], numeric(1))
    calib_w_eff  <- vapply(seq_len(nE), function(i) as.numeric(cal[g_idx[i]]),  numeric(1))

    # isotonic projection of realized noncentrality (nondecreasing in effect),
    # for reporting/diagnostics only -- it does NOT alter the selected lambda_cond.
    snr_iso <- realized_snr
    if (nE >= 2L) {
      y <- realized_snr
      if (any(!is.finite(y))) y[!is.finite(y)] <- min(y[is.finite(y)], na.rm = TRUE)
      snr_iso <- stats::isoreg(effs, y)$yf
    }

    by_effect <- data.frame(
      effect = effs, lambda_p = best_lp, lambda_cond = chosen_lc,
      SNR = realized_snr, SNR_iso = snr_iso,
      calib_weight = calib_w_eff
    )
    names(by_effect)[names(by_effect) == "effect"] <- effect_label
    names(curves)[names(curves) == "effect"]       <- effect_label

    # representative scalar + shared flag
    if (identical(lambda_cond_selection, "shared")) {
      shared_lc <- lambda_c_grid[g_idx[1]]         # the single shared weight
      rep_lc    <- shared_lc
    } else {
      shared_lc <- NA_real_                        # no single shared value
      med_i     <- which.min(abs(effs - stats::median(effs)))  # representative = median effect
      rep_lc    <- chosen_lc[med_i]
    }

    list(by_effect = by_effect, curves = curves, rep_lc = rep_lc,
         shared_lc = shared_lc,
         diagnostics = list(calibration = calib_tbl,
                            weights = data.frame(lambda_cond = lambda_c_grid,
                                                 calib_weight = as.numeric(cal))))
  }

  # ---- Stage 2A: SaLFC -------------------------------------------------------
  if (verbose) message(sprintf("\n>>> [Stage 2A] SaLFC lambda_cond (lambda_p = %.2f)", best_lp))
  salfc <- tune_endpoint("SaLFC", de_grid, "DE_lfc",
                         de_genes_h1 = G_DE, null_prop = NULL)

  # ---- Stage 2B: LOR ---------------------------------------------------------
  if (has_lor) {
    if (verbose) message(sprintf("\n>>> [Stage 2B] LOR lambda_cond (lambda_p = %.2f)", best_lp))
    lor <- tune_endpoint("LOR", tpc_grid, "target_prop_case",
                         de_genes_h1 = character(0), null_prop = lor_null_prop)
  } else {
    empty_be <- data.frame(target_prop_case = numeric(0), lambda_p = numeric(0),
                           lambda_cond = numeric(0), SNR = numeric(0), SNR_iso = numeric(0),
                           calib_weight = numeric(0))
    lor <- list(by_effect = empty_be, curves = data.frame(), rep_lc = NA_real_,
                shared_lc = NA_real_,
                diagnostics = list(calibration = data.frame(), weights = data.frame()))
  }

  res <- list(
    optimal_params = list(
      SaLFC = list(lambda_p = best_lp, lambda_cond = salfc$rep_lc,
                   shared_lambda_cond = salfc$shared_lc, by_effect = salfc$by_effect),
      LOR   = list(lambda_p = best_lp, lambda_cond = lor$rep_lc,
                   shared_lambda_cond = lor$shared_lc, by_effect = lor$by_effect)
    ),
    summary_stage1 = stage1_summary,
    summary_stage2 = list(SaLFC = salfc$curves, LOR = lor$curves),
    diagnostics    = list(SaLFC = salfc$diagnostics, LOR = lor$diagnostics),
    settings       = list(objective = "snr", selection = lambda_cond_selection,
                          eps_floor = eps_floor, alpha = alpha, calib_mult = calib_mult,
                          K_syn = K_syn, calibration_guard = calibration_guard),
    effect_grids   = list(DE_lfc = de_grid, target_prop_case = if (has_lor) tpc_grid else NULL),
    call = match.call()
  )
  class(res) <- "spaCraftTuning"
  res
}

# ==============================================================================
# S3 Methods (Summary & Plot)
# ==============================================================================

#' @title Summary Method for spaCraftTuning
#' @param object An object of class \code{spaCraftTuning}.
#' @param ... Additional arguments affecting the summary produced.
#' @method summary spaCraftTuning
#' @export
summary.spaCraftTuning <- function(object, ...) {
  cat("\nCall:\n"); print(object$call)
  st <- object$settings
  cat("\n========================================================")
  cat("\n Stability- and Validity-Driven Hyperparameter Tuning")
  cat("\n========================================================\n")
  cat(sprintf(" objective = %s  |  selection = %s\n", st$objective, st$selection))
  cat(sprintf(" calibration guard = %s (alpha = %.3f)\n",
              if (st$calibration_guard) "on" else "off", st$alpha))
  cat("--------------------------------------------------------\n")

  print_tbl <- function(df) {
    df[] <- lapply(df, function(z) if (is.numeric(z)) round(z, 3) else z)
    print(df, row.names = FALSE)
  }

  cat(" [ Endpoint: SaLFC (Gene Expression Power) ]\n")
  cat(sprintf("   - lambda_p (shared, Stage 1 ARI)          : %.2f\n",
              object$optimal_params$SaLFC$lambda_p))
  if (!is.na(object$optimal_params$SaLFC$shared_lambda_cond)) {
    cat(sprintf("   - lambda_cond (shared across DE_lfc)      : %.2f\n",
                object$optimal_params$SaLFC$shared_lambda_cond))
  } else {
    cat(sprintf("   - lambda_cond (per-effect; median rep.)   : %.2f\n",
                object$optimal_params$SaLFC$lambda_cond))
  }
  cat("   - per-DE_lfc optimal lambda_cond & noncentrality:\n")
  print_tbl(object$optimal_params$SaLFC$by_effect[,
                                                  c("DE_lfc", "lambda_cond", "SNR", "SNR_iso", "calib_weight")])
  cat("\n")

  if (nrow(object$optimal_params$LOR$by_effect) > 0) {
    cat(" [ Endpoint: LOR (Composition Power) ]\n")
    cat(sprintf("   - lambda_p (shared, Stage 1 ARI)          : %.2f\n",
                object$optimal_params$LOR$lambda_p))
    if (!is.na(object$optimal_params$LOR$shared_lambda_cond)) {
      cat(sprintf("   - lambda_cond (shared across tpc)         : %.2f\n",
                  object$optimal_params$LOR$shared_lambda_cond))
    } else {
      cat(sprintf("   - lambda_cond (per-effect; median rep.)   : %.2f\n",
                  object$optimal_params$LOR$lambda_cond))
    }
    cat("   - per-target_prop_case optimal lambda_cond & noncentrality:\n")
    print_tbl(object$optimal_params$LOR$by_effect[,
                                                  c("target_prop_case", "lambda_cond", "SNR", "SNR_iso", "calib_weight")])
  } else {
    cat(" [ Endpoint: LOR ] : skipped (no target_prop_case supplied)\n")
  }
  cat("========================================================\n")
  cat(" SNR = regularized noncentrality  E[T] / (sd[T] + eps*s0) = 1 / floored CV.\n")
  cat(" SNR_iso = isotonic (nondecreasing-in-effect) projection of SNR.\n")
  cat(" lambda_cond is selected per effect size (selection = 'per_effect') and\n")
  cat(" matched per effect by evaluatePowerSaLFC()/evaluatePowerLOR().\n")
  cat("========================================================\n\n")
  invisible(object)
}

# Silence R CMD check NOTES for ggplot2 non-standard evaluation variables
utils::globalVariables(c("lambda_p", "ARI", "lambda_cond", "SNR", "Var_T",
                         "mean_T", "sd_T", "DE_lfc", "target_prop_case"))

#' @title Plot Method for spaCraftTuning
#' @param x An object of class \code{spaCraftTuning}.
#' @param ... Additional arguments to be passed to plotting functions.
#' @method plot spaCraftTuning
#' @import ggplot2
#' @import patchwork
#' @export
plot.spaCraftTuning <- function(x, ...) {
  base_theme <- ggplot2::theme_bw() +
    ggplot2::theme(plot.title  = ggplot2::element_text(face = "bold", size = 11, hjust = 0.5),
                   axis.title  = ggplot2::element_text(size = 10),
                   legend.position = "bottom")

  sel_mode <- if (is.null(x$settings$selection)) "shared" else x$settings$selection
  sel_txt  <- if (identical(sel_mode, "per_effect")) "effect optimum" else "shared optimum"

  # Stage 1: ARI vs lambda_p (single shared recovery weight)
  p1 <- ggplot2::ggplot(x$summary_stage1, ggplot2::aes(x = lambda_p, y = ARI)) +
    ggplot2::geom_line(color = "black", linewidth = 1) +
    ggplot2::geom_point(size = 3) +
    ggplot2::geom_vline(xintercept = x$optimal_params$SaLFC$lambda_p,
                        color = "red", linetype = "dashed", linewidth = 1) +
    ggplot2::labs(title = expression(paste("Stage 1: Domain Recovery (", lambda[p], ")")),
                  subtitle = "Maximizing pairwise ARI",
                  x = expression(paste("pBANKSY ", lambda[p])), y = "Mean ARI") +
    base_theme

  # Stage 2A: SaLFC, one noncentrality curve per DE_lfc value; the selected
  # (per-effect) lambda_cond is marked with a dashed rule + a hollow point.
  s2s <- x$summary_stage2$SaLFC
  s2s$DE_lfc <- factor(s2s$DE_lfc)
  opt_s <- x$optimal_params$SaLFC$by_effect
  opt_s$DE_lfc <- factor(opt_s$DE_lfc)

  p2 <- ggplot2::ggplot(s2s, ggplot2::aes(x = lambda_cond, y = SNR, color = DE_lfc, group = DE_lfc)) +
    ggplot2::geom_line(linewidth = 1) + ggplot2::geom_point(size = 2) +
    ggplot2::geom_vline(data = opt_s,
                        ggplot2::aes(xintercept = lambda_cond, color = DE_lfc),
                        linetype = "dashed", linewidth = 0.8, alpha = 0.8, show.legend = FALSE) +
    ggplot2::geom_point(data = opt_s,
                        ggplot2::aes(x = lambda_cond, y = SNR, color = DE_lfc),
                        shape = 21, fill = "white", size = 3.6, stroke = 1.3,
                        show.legend = FALSE) +
    ggplot2::labs(title = expression(paste("Stage 2: Expr Power (", lambda[cond], ")")),
                  subtitle = bquote("Max noncentrality"~T[SaLFC]~"("*.(sel_txt)*")"),
                  x = expression(lambda[cond]),
                  y = expression(paste("SNR(", T[SaLFC], ") = 1 / CV")),
                  color = "DE_lfc") +
    base_theme

  panels <- p1 | p2

  # Stage 2B: LOR, one noncentrality curve per target_prop_case value (if present)
  if (nrow(x$optimal_params$LOR$by_effect) > 0) {
    s2l <- x$summary_stage2$LOR
    s2l$target_prop_case <- factor(s2l$target_prop_case)
    opt_l <- x$optimal_params$LOR$by_effect
    opt_l$target_prop_case <- factor(opt_l$target_prop_case)

    p3 <- ggplot2::ggplot(s2l, ggplot2::aes(x = lambda_cond, y = SNR,
                                            color = target_prop_case, group = target_prop_case)) +
      ggplot2::geom_line(linewidth = 1) + ggplot2::geom_point(size = 2) +
      ggplot2::geom_vline(data = opt_l,
                          ggplot2::aes(xintercept = lambda_cond, color = target_prop_case),
                          linetype = "dashed", linewidth = 0.8, alpha = 0.8, show.legend = FALSE) +
      ggplot2::geom_point(data = opt_l,
                          ggplot2::aes(x = lambda_cond, y = SNR, color = target_prop_case),
                          shape = 21, fill = "white", size = 3.6, stroke = 1.3,
                          show.legend = FALSE) +
      ggplot2::labs(title = expression(paste("Stage 2: Comp Power (", lambda[cond], ")")),
                    subtitle = bquote("Max noncentrality"~T[LOR]~"("*.(sel_txt)*")"),
                    x = expression(lambda[cond]),
                    y = expression(paste("SNR(", T[LOR], ") = 1 / CV")),
                    color = "target_prop_case") +
      base_theme
    panels <- p1 | p2 | p3
  }

  panels + patchwork::plot_annotation(
    title = "spaCraft Hyperparameter Tuning Diagnostics",
    theme = ggplot2::theme(plot.title = ggplot2::element_text(size = 14, face = "bold", hjust = 0.5))
  )
}
