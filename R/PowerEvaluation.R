#' @title Evaluate SaLFC Power over a Grid
#' @description
#' Grid simulation over sample sizes (K) and effect sizes (DE_lfc) to evaluate
#' the power and FDP of the SaLFC endpoint. Per replicate the pipeline is
#' \code{simulatespaCraft()} -> \code{pBANSKY()} ->
#' \code{rearrangeSyntheticToPilot()} -> \code{SaLFC()}. Synthetic coordinates
#' are morphed onto the pilot domain geometry by \code{eval_domain_col}
#' ("hat_d" by default), and SaLFC selects target/reference by the SAME column,
#' so the spatial variance \eqn{\tau^2} is computed on a consistent basis.
#' Pilot TREAT thresholds (\code{lfc_threshold = "auto"}) are estimated once per
#' DE_lfc and shared across replicates. A progress bar is shown and simulation
#' errors are skipped gracefully.
#'
#' @param object A spaCraft object (fitted: pilot + params).
#' @param K_grid Numeric vector of samples per group (e.g. \code{seq(3, 8)}).
#' @param lfc_grid Numeric vector of effect sizes (DE_lfc, on the count-model
#'   log scale used by \code{simulatespaCraft}).
#' @param n_sim Integer. Replications per grid point.
#' @param seed Integer. Master seed.
#' @param G_svg Character. All genes to simulate (\code{sets$G_svg}).
#' @param G_null Character. Null genes to evaluate (\code{sets$G_null}).
#' @param G_spike Character. True spike-in DE genes (\code{sets$G_spike}).
#' @param target_domain,reference_domain Domain labels.
#' @param eval_domain_col Label column for morphing AND testing: "hat_d"
#'   (recovered clusters, default) or "domain" (ground truth).
#' @param rearrange Logical. Morph synthetic blobs onto pilot geometry first.
#' @param jitter_frac,k_nn Rearrangement tuning (used when \code{rearrange}).
#' @param test_method "t" or "wilcoxon".
#' @param adjust_method "BH" or "bonferroni" (default).
#' @param alternative SaLFC direction: "greater" (default), "less", "two.sided".
#' @param alpha_fdr Numeric. Rejection threshold for power/FDP.
#' @param lfc_threshold "auto" (pilot plug-in), a scalar, or a named vector.
#' @param lfc_q Numeric. Quantile for the auto threshold.
#' @param spatial_mode "original" or "basis".
#' @param n_cores Integer. Cores for simulate/cluster (keep modest under nesting).
#' @param scenario_base List. Base scenario (composition kept neutral by default
#'   so the swept effect is purely expression DE).
#' @param clustering_args List. Passed to \code{pBANSKY} (must contain
#'   \code{lambda}).
#' @param tuning_obj Optional \code{spaCraftTuning} object from
#'   \code{\link{tuneSpaDesignLambdas}}. When supplied, the conditional-texture
#'   weight \code{lambda_cond} (in the simulation scenario) and the clustering
#'   weight \code{lambda} (in \code{clustering_args}) are taken \emph{per effect
#'   size} from the matching row of \code{optimal_params$SaLFC$by_effect}
#'   (matched on \code{DE_lfc}; nearest value with a warning if there is no exact
#'   match), overriding the values in \code{scenario_base}/\code{clustering_args}.
#'   Default \code{NULL}.
#' @param ... Forwarded simulate tuning: \code{graph_k}, \code{gp_rank},
#'   \code{gp_jitter}, \code{cond_k_nn}.
#' @return A data.frame of per-replicate power/FDP rows, tagged with class
#'   \code{"powerSaLFC"} (see \code{\link{print.powerSaLFC}},
#'   \code{\link{summary.powerSaLFC}}, \code{\link{plotPowerCurveSaLFC}}).
#'
#' @seealso \code{\link{evaluatePowerLOR}}, \code{\link{SaLFC}},
#'   \code{\link{tuneSpaDesignLambdas}}
#'
#' @importFrom utils txtProgressBar setTxtProgressBar modifyList
#' @importFrom stats setNames
#' @importFrom dplyr bind_rows
#' @export
#'
evaluatePowerSaLFC <- function(object,
                               K_grid,
                               lfc_grid,
                               n_sim = 30,
                               seed = 1,
                               G_svg,
                               G_null,
                               G_spike,
                               target_domain,
                               reference_domain,
                               eval_domain_col  = c("hat_d", "domain"),
                               rearrange        = TRUE,
                               jitter_frac      = 0.3,
                               k_nn             = 100L,
                               test_method      = c("t", "wilcoxon"),
                               adjust_method    = c("bonferroni", "BH"),
                               alternative      = c("greater", "less", "two.sided"),
                               alpha_fdr        = 0.05,
                               lfc_threshold    = "auto",
                               lfc_q            = 0.5,
                               spatial_mode     = c("original", "basis"),
                               n_cores          = 1L,
                               scenario_base = list(target_prop_case = NULL,
                                                    delta_pp = NULL, kappa_m = 1.0,
                                                    kappa_bio = 1.0, delta_rho = 1.0,
                                                    lambda_cond = 0.5),
                               clustering_args = list(lambda = 0.5, k_neighbors = 20,
                                                      scale_pilot_coords = FALSE,
                                                      scale_sim_coords = FALSE,
                                                      hard_guidance = FALSE, n_pcs = 20,
                                                      do_hungarian = TRUE, use_G = TRUE,
                                                      balance_mode = "energy_proportional"),
                               tuning_obj = NULL,
                               ...) {

  test_method     <- match.arg(test_method)
  adjust_method   <- match.arg(adjust_method)
  alternative     <- match.arg(alternative)
  spatial_mode    <- match.arg(spatial_mode)
  eval_domain_col <- match.arg(eval_domain_col)

  ea        <- list(...)
  graph_k   <- ea$graph_k   %||% 20L
  gp_rank   <- ea$gp_rank   %||% 20L
  gp_jitter <- ea$gp_jitter %||% 1e-8
  cond_k_nn <- ea$cond_k_nn %||% 20L

  # CRAN-safe RNG: restore the caller's stream on exit.
  if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
    .old_seed <- get(".Random.seed", envir = globalenv(), inherits = FALSE)
    on.exit(assign(".Random.seed", .old_seed, envir = globalenv()), add = TRUE)
  }

  set.seed(seed)
  total_steps <- length(lfc_grid) * length(K_grid) * n_sim

  # --------------------------------------------------------------------------
  # 1. Precompute pilot TREAT thresholds per DE_lfc
  # --------------------------------------------------------------------------
  threshold_cache <- list()
  if (is.character(lfc_threshold) && length(lfc_threshold) == 1L && lfc_threshold == "auto") {
    message(">>> Precomputing gene-specific pilot thresholds by DE_lfc ...")
    tau_gene <- setdiff(G_null, G_spike); if (length(tau_gene) == 0) tau_gene <- G_null
    for (lfc in lfc_grid) {
      cap <- if (is.finite(lfc) && lfc > 0) 0.90 * lfc else NULL
      threshold_cache[[sprintf("%.12g", lfc)]] <- estimate_lfc_threshold_from_pilot(
        object = object, genes = tau_gene,
        target_domain = target_domain, reference_domain = reference_domain, rho_TR = 0,
        adjust_method = adjust_method, test_method = test_method,
        alternative = alternative, q = lfc_q, cap_effect = cap, verbose = FALSE)
    }
  }

  # --------------------------------------------------------------------------
  # 2. Main grid
  # --------------------------------------------------------------------------
  thr_label <- if (is.character(lfc_threshold)) lfc_threshold else "scalar/vector"
  message(sprintf(">>> Power grid (steps: %d | test: %s | thr: %s | rearrange: %s | eval: %s)",
                  total_steps, test_method, thr_label, rearrange, eval_domain_col))
  pb <- utils::txtProgressBar(min = 0, max = total_steps, style = 3, width = 50)

  sum_list <- vector("list", total_steps)
  step <- 0L

  for (lfc in lfc_grid) {
    resolved_thr <- if (length(threshold_cache) > 0)
      threshold_cache[[sprintf("%.12g", lfc)]] else lfc_threshold

    # Per-effect spatial weights: if a tuning object is supplied, use the
    # lambda_cond / lambda tuned for THIS DE_lfc (matched in by_effect),
    # overriding scenario_base$lambda_cond and clustering_args$lambda.
    lc_use  <- scenario_base$lambda_cond
    lam_use <- clustering_args$lambda
    if (!is.null(tuning_obj)) {
      bd <- tuning_obj$optimal_params$SaLFC$by_effect
      j  <- which.min(abs(bd$DE_lfc - lfc))
      if (!isTRUE(all.equal(bd$DE_lfc[j], lfc)))
        warning(sprintf("evaluatePowerSaLFC: no tuned lambda for DE_lfc=%.4g; using nearest (%.4g).",
                        lfc, bd$DE_lfc[j]))
      lc_use  <- bd$lambda_cond[j]
      lam_use <- bd$lambda_p[j]
    }
    clustering_args_lfc <- utils::modifyList(clustering_args, list(lambda = lam_use))

    for (K in K_grid) {
      for (b in seq_len(n_sim)) {
        step <- step + 1L
        this_seed <- seed + 10000L * match(lfc, lfc_grid) + 100L * match(K, K_grid) + b
        set.seed(this_seed)

        one_res <- tryCatch({
          # A. Simulate (expression DE = lfc; composition neutral)
          scen <- utils::modifyList(scenario_base, list(DE_lfc = lfc, lambda_cond = lc_use))
          sim_res <- simulatespaCraft(
            object = object, n_sample_per_group = K, scenario_settings = scen,
            target_domain = target_domain, genes_to_simulate = G_svg,
            de_genes = G_spike, spatial_mode = spatial_mode, graph_k = graph_k,
            cond_k_nn = cond_k_nn, gp_rank = gp_rank, gp_jitter = gp_jitter,
            verbose = FALSE, n_cores = n_cores, seed_base = this_seed)

          # B. Recover domains (pBANSKY; lambda lives in clustering_args)
          syn_d_hat <- do.call(pBANSKY, c(
            list(sim_data = sim_res, pilot_data = object@pilot_data,
                 verbose = FALSE, n_cores = n_cores), clustering_args_lfc))

          syn_fmt <- stats::setNames(
            lapply(syn_d_hat, function(it)
              list(counts = it$counts, coords = it$coords,
                   meta = data.frame(sample_id = it$sample_id, group = it$group,
                                     stringsAsFactors = FALSE))),
            vapply(syn_d_hat, function(it) it$sample_id, character(1)))

          # B2. Rearrange onto pilot geometry (canonical test coordinate)
          if (isTRUE(rearrange)) {
            syn_fmt <- rearrangeSyntheticToPilot(
              synthetic_data = syn_fmt, pilot_data = object@pilot_data,
              match_by = eval_domain_col, jitter_frac = jitter_frac,
              k_nn = k_nn, seed = this_seed, verbose = FALSE)
          }

          # C. SaLFC on the (rearranged) coords, by eval_domain_col
          res_df <- SaLFC(
            processed_data = syn_fmt, object = object,
            genes = unique(c(G_null, G_spike)),
            target_domain = target_domain, reference_domain = reference_domain,
            domain_col = eval_domain_col, rho_TR = 0,
            adjust_method = adjust_method, test_method = test_method,
            alternative = alternative, lfc_threshold = resolved_thr)

          # D. Power / FDP
          res_df$is_DE_truth <- res_df$gene %in% G_spike
          res_df$reject <- is.finite(res_df$padj) & (res_df$padj <= alpha_fdr)
          n_DE  <- sum(res_df$is_DE_truth)
          n_TP  <- sum(res_df$reject &  res_df$is_DE_truth)
          n_FP  <- sum(res_df$reject & !res_df$is_DE_truth)
          n_rej <- sum(res_df$reject)
          strict_ok    <- (n_FP == 0)
          power_raw    <- if (n_DE > 0) n_TP / n_DE else NA_real_
          power_strict <- if (strict_ok) power_raw else 0
          fdp          <- if (n_rej > 0) n_FP / n_rej else 0

          data.frame(
            K_syn = K, DE_lfc = lfc, adjust_method = adjust_method,
            test_method = test_method, eval_domain_col = eval_domain_col,
            rearrange = rearrange, strict_ok = strict_ok,
            power_strict = power_strict, power_raw = power_raw, fdp = fdp,
            n_rej = n_rej, n_TP = n_TP, n_FP = n_FP, n_test = nrow(res_df),
            n_DE = n_DE, spatial_mode = spatial_mode, rep = b,
            seed_used = this_seed, stringsAsFactors = FALSE)
        }, error = function(e) NULL)

        if (!is.null(one_res)) sum_list[[step]] <- one_res
        utils::setTxtProgressBar(pb, step)
      }
    }
  }

  close(pb)
  message(">>> Grid completed.")
  out <- dplyr::bind_rows(sum_list)
  attr(out, "alpha") <- alpha_fdr            # 표시용 (선택)
  attr(out, "n_sim") <- n_sim                # 표시용 (선택)
  class(out) <- c("powerSaLFC", class(out))  # 핵심: 이 한 줄
  out
}

#' @title Evaluate Composition (LOR) Power over a Grid
#' @description
#' Grid simulation over sample sizes (K) and composition effect sizes to evaluate
#' the power and empirical size of the sample-level log-odds-ratio test
#' \code{\link{LOR}}. Per replicate the pipeline is \code{simulatespaCraft()} ->
#' \code{pBANSKY()} -> \code{LOR()}. The effect size is the case-vs-control shift
#' in the \code{target_domain} abundance, injected through the simulator's
#' composition scenario (\code{delta_pp}, in percentage points, or
#' \code{target_prop_case}); expression DE is switched off (\code{DE_lfc = 0},
#' \code{de_genes = character(0)}) so the swept signal is purely compositional.
#'
#' Unlike the per-gene SaLFC test, \code{LOR} returns a SINGLE test per
#' (target, reference) pair, so there is no multiple testing within a replicate:
#' \strong{power = rejection rate} across replicates at each (K, effect), and the
#' \code{effect = 0} row is the \strong{empirical size}. The composition TREAT
#' threshold \eqn{\tau} does not depend on the effect size, so under
#' \code{null_type = "treat"} the pilot \eqn{\tau} is estimated ONCE (a scalar)
#' and reused. A progress bar is shown and simulation errors are skipped.
#'
#' @param object A spaCraft object (fitted: pilot, params, \code{@params_composition}).
#' @param K_grid Numeric vector of samples per group (e.g. \code{seq(3, 8)}).
#' @param delta_grid Numeric vector of composition effect sizes. When
#'   \code{effect_type = "delta_pp"} these are percentage-point shifts of the
#'   target proportion (0 = null); when \code{"target_prop_case"} they are the
#'   absolute case-group target proportions in (0, 1).
#' @param n_sim Integer. Replications per grid point.
#' @param seed Integer. Master seed.
#' @param genes Character vector of genes to simulate (e.g. \code{sets$G_svg});
#'   needed only so \code{pBANSKY} has counts to cluster. Expression DE is off.
#' @param target_domain,reference_domain Domain labels.
#' @param effect_type "delta_pp" or "target_prop_case" (default).
#' @param eval_on Label column \code{LOR} counts by: "hat_d" (recovered clusters,
#'   default) or "domain" (ground-truth oracle).
#' @param null_type "treat" (default; centered at the pilot baseline
#'   \eqn{\delta_0=\beta_1} with margin \eqn{\tau}) or "two_sample"
#'   (\eqn{H_0:\ \bar L_1=\bar L_0}, \eqn{\tau} forced to 0).
#' @param tau "auto" (default; pilot plug-in scalar via
#'   \code{estimate_lor_tau_from_pilot}, built from \code{null_type},
#'   \code{tau_q}, \code{target_prop_min} and \code{tau_args}) or a numeric
#'   scalar. Used only when \code{null_type = "treat"}.
#' @param tau_q Numeric in (0,1). Quantile level for the auto \eqn{\tau}
#'   noise-floor multiplier \eqn{t_{q,df}}. Default 0.5. \strong{NB:}
#'   \code{tau_q = 0.5} gives \eqn{t_{0.5,df}=0}, hence a zero noise floor; with
#'   \code{target_prop_min = NULL} the auto \eqn{\tau} then collapses to 0.
#' @param target_prop_min Numeric in (0,1) or NULL (default). Forwarded to
#'   \code{estimate_lor_tau_from_pilot()} to set the scientific-anchor
#'   \eqn{\tau_{sci}} = the (direction-aware) log-odds gap of \code{target_prop_min}
#'   from the fitted case baseline \code{mu1_hat} = \code{plogis(beta0 + beta1)}
#'   (under \code{null_type = "treat"}). NULL leaves the auto \eqn{\tau} as the
#'   pilot noise floor alone. Only used when \code{tau = "auto"}.
#' @param tau_args Named list of extra arguments forwarded to
#'   \code{estimate_lor_tau_from_pilot()} for the auto \eqn{\tau}, e.g.
#'   \code{list(use_noise_floor = FALSE)} for a pure fitted anchor, or
#'   \code{list(cap_effect = 0.5, min_tau = 0.02, q = 0.9)}. Entries here
#'   override the defaults built from \code{null_type}/\code{tau_q}/\code{target_prop_min}.
#' @param test_method "t" (Welch, default) or "wilcoxon".
#' @param alternative "greater" (default), "two.sided", or "less".
#' @param alpha Numeric. Rejection threshold.
#' @param eps Numeric > 0. Continuity correction forwarded to \code{LOR}.
#' @param rearrange Logical. Morph synthetic coords onto pilot geometry first.
#'   Default FALSE: \code{LOR} counts by label and is invariant to the spatial
#'   rearrangement, so this only costs time.
#' @param jitter_frac,k_nn Rearrangement tuning (used only when \code{rearrange}).
#' @param n_cores Integer. Cores for simulate/cluster.
#' @param scenario_base List. Base scenario; any \code{delta_pp},
#'   \code{target_prop_case}, or \code{DE_lfc} here are dropped (they are
#'   controlled by the sweep).
#' @param clustering_args List passed to \code{pBANSKY} (must contain
#'   \code{lambda}). NB: for composition power the recovered \code{hat_d} must
#'   preserve per-sample target abundance, so avoid size-equalizing balance modes.
#' @param tuning_obj Optional \code{spaCraftTuning} object from
#'   \code{\link{tuneSpaDesignLambdas}}. When supplied, \code{lambda_cond}
#'   (scenario) and \code{lambda} (\code{clustering_args}) are taken \emph{per
#'   effect size} from \code{optimal_params$LOR$by_effect} (matched on
#'   \code{target_prop_case}; nearest value with a warning if no exact match).
#'   Only applied when \code{effect_type = "target_prop_case"}. Default \code{NULL}.
#' @param ... Forwarded simulate tuning: \code{graph_k}, \code{gp_rank},
#'   \code{gp_jitter}, \code{cond_k_nn}.
#' @return A data.frame of per-replicate rejection rows, tagged with class
#'   \code{"powerLOR"} (see \code{\link{print.powerLOR}},
#'   \code{\link{summary.powerLOR}}, \code{\link{plotPowerCurveLOR}}).
#'
#' @seealso \code{\link{evaluatePowerSaLFC}}, \code{\link{LOR}},
#'   \code{\link{estimate_lor_tau_from_pilot}}, \code{\link{tuneSpaDesignLambdas}}
#'
#' @importFrom utils txtProgressBar setTxtProgressBar modifyList
#' @importFrom stats setNames
#' @importFrom dplyr bind_rows
#' @export
#'
evaluatePowerLOR <- function(object,
                             K_grid,
                             delta_grid,
                             n_sim = 30,
                             seed = 1,
                             genes,
                             target_domain,
                             reference_domain,
                             effect_type = c("target_prop_case", "delta_pp"),
                             eval_on          = c("hat_d", "domain"),
                             null_type        = c("treat", "two_sample"),
                             tau              = "auto",
                             tau_q            = 0.5,
                             target_prop_min  = NULL,
                             tau_args         = list(),
                             test_method      = c("t", "wilcoxon"),
                             alternative      = c("greater", "two.sided", "less"),
                             alpha            = 0.05,
                             eps              = 0.5,
                             rearrange        = FALSE,
                             jitter_frac      = 0.3,
                             k_nn             = 50L,
                             n_cores          = 1L,
                             scenario_base = list(kappa_m = 1.0, kappa_bio = 1.0,
                                                  delta_rho = 1.0, lambda_cond = 0.5),
                             clustering_args = list(lambda = 0.5, k_neighbors = 20,
                                                    scale_pilot_coords = FALSE,
                                                    scale_sim_coords = FALSE,
                                                    hard_guidance = FALSE, n_pcs = 20,
                                                    do_hungarian = TRUE, use_G = TRUE,
                                                    balance_mode = "energy_proportional"),
                             tuning_obj = NULL,
                             ...) {

  effect_type <- match.arg(effect_type)
  eval_on     <- match.arg(eval_on)
  null_type   <- match.arg(null_type)
  test_method <- match.arg(test_method)
  alternative <- match.arg(alternative)

  ea        <- list(...)
  graph_k   <- ea$graph_k   %||% 20L
  gp_rank   <- ea$gp_rank   %||% 20L
  gp_jitter <- ea$gp_jitter %||% 1e-8
  cond_k_nn <- ea$cond_k_nn %||% 20L

  # the sweep owns the composition/expression effect keys
  scenario_base[c("delta_pp", "target_prop_case", "DE_lfc")] <- NULL

  # CRAN-safe RNG: restore the caller's stream on exit.
  if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
    .old_seed <- get(".Random.seed", envir = globalenv(), inherits = FALSE)
    on.exit(assign(".Random.seed", .old_seed, envir = globalenv()), add = TRUE)
  }

  set.seed(seed)
  total_steps <- length(delta_grid) * length(K_grid) * n_sim

  # --------------------------------------------------------------------------
  # 1. Composition TREAT threshold (effect-independent -> a single scalar)
  #    The auto path forwards null_type / tau_q (as q) / target_prop_min, plus
  #    any overrides in tau_args, to estimate_lor_tau_from_pilot().
  # --------------------------------------------------------------------------
  if (identical(null_type, "treat")) {
    if (is.numeric(tau) && length(tau) == 1L) {
      tau_fixed <- as.numeric(tau)
    } else if (is.character(tau) && length(tau) == 1L && tau == "auto") {
      message(">>> Precomputing pilot composition tau (single scalar) ...")
      tau_call <- utils::modifyList(
        list(object = object, target_domain = target_domain,
             reference_domain = reference_domain, null_type = null_type,
             test_method = test_method, alternative = alternative,
             q = tau_q, target_prop_min = target_prop_min, eps = eps,
             verbose = FALSE),
        tau_args)
      tau_fixed <- do.call(estimate_lor_tau_from_pilot, tau_call)
    } else stop("tau must be 'auto' or a numeric scalar.")
  } else {
    tau_fixed <- 0  # two_sample: H0 bar1 = bar0
  }

  # --------------------------------------------------------------------------
  # 2. Main grid
  # --------------------------------------------------------------------------
  message(sprintf(">>> LOR power grid (steps: %d | test: %s | null: %s | eval: %s | tau: %s)",
                  total_steps, test_method, null_type, eval_on,
                  if (null_type == "treat") sprintf("%.4f", tau_fixed) else "-"))
  pb <- utils::txtProgressBar(min = 0, max = total_steps, style = 3, width = 50)

  sum_list <- vector("list", total_steps)
  step <- 0L

  for (d in delta_grid) {
    # Per-effect spatial weights: if a tuning object is supplied, use the
    # lambda_cond / lambda tuned for THIS effect (matched on target_prop_case),
    # overriding scenario_base$lambda_cond and clustering_args$lambda.
    lc_use  <- scenario_base$lambda_cond
    lam_use <- clustering_args$lambda
    if (!is.null(tuning_obj)) {
      if (identical(effect_type, "target_prop_case")) {
        bd <- tuning_obj$optimal_params$LOR$by_effect
        j  <- which.min(abs(bd$target_prop_case - d))
        if (!isTRUE(all.equal(bd$target_prop_case[j], d)))
          warning(sprintf("evaluatePowerLOR: no tuned lambda for target_prop_case=%.4g; using nearest (%.4g).",
                          d, bd$target_prop_case[j]))
        lc_use  <- bd$lambda_cond[j]
        lam_use <- bd$lambda_p[j]
      } else {
        warning("evaluatePowerLOR: tuned lambdas are keyed by target_prop_case; with effect_type='delta_pp' they are not matched (using scenario_base/clustering_args).")
      }
    }
    clustering_args_d <- utils::modifyList(clustering_args, list(lambda = lam_use))

    for (K in K_grid) {
      for (b in seq_len(n_sim)) {
        step <- step + 1L
        this_seed <- seed + 10000L * match(d, delta_grid) + 100L * match(K, K_grid) + b
        set.seed(this_seed)

        one_res <- tryCatch({
          # A. Simulate composition shift; expression DE off
          scen <- utils::modifyList(
            scenario_base,
            if (effect_type == "delta_pp") list(delta_pp = d, lambda_cond = lc_use)
            else list(target_prop_case = d, lambda_cond = lc_use))
          scen$DE_lfc <- 0
          sim_res <- simulatespaCraft(
            object = object, n_sample_per_group = K, scenario_settings = scen,
            target_domain = target_domain, genes_to_simulate = genes,
            de_genes = character(0), spatial_mode = "original", graph_k = graph_k,
            cond_k_nn = cond_k_nn, gp_rank = gp_rank, gp_jitter = gp_jitter,
            verbose = FALSE, n_cores = n_cores, seed_base = this_seed)

          # B. Recover domains (hat_d)
          syn_d_hat <- do.call(pBANSKY, c(
            list(sim_data = sim_res, pilot_data = object@pilot_data,
                 verbose = FALSE, n_cores = n_cores), clustering_args_d))

          syn_fmt <- stats::setNames(
            lapply(syn_d_hat, function(it) {
              co <- it$coords
              if (identical(eval_on, "domain")) co$hat_d <- NULL  # force ground truth
              list(coords = co,
                   meta = data.frame(sample_id = it$sample_id, group = it$group,
                                     stringsAsFactors = FALSE))
            }),
            vapply(syn_d_hat, function(it) it$sample_id, character(1)))

          # B2. Optional rearrange (LOR is label-count invariant to coordinates)
          if (isTRUE(rearrange)) {
            syn_fmt <- rearrangeSyntheticToPilot(
              synthetic_data = syn_fmt, pilot_data = object@pilot_data,
              match_by = if (eval_on == "domain") "domain" else "hat_d",
              jitter_frac = jitter_frac, k_nn = k_nn, seed = this_seed, verbose = FALSE)
          }

          # C. Composition test (single p-value)
          res <- LOR(
            processed_data = syn_fmt, object = object,
            target_domain = target_domain, reference_domain = reference_domain,
            test_method = test_method, alternative = alternative,
            null_type = null_type, tau = tau_fixed, eps = eps, verbose = FALSE)

          # D. Reject?
          p   <- res$p_value %||% NA_real_
          rej <- is.finite(p) && (p <= alpha)

          data.frame(
            K_syn = K, effect = d, effect_type = effect_type,
            delta_pp = if (effect_type == "delta_pp") d else NA_real_,
            target_prop_case = if (effect_type == "target_prop_case") d else NA_real_,
            null_type = null_type, eval_on = eval_on, test_method = test_method,
            reject = rej, p_value = p,
            delta_hat = res$delta_hat %||% NA_real_, delta0 = res$delta0 %||% NA_real_,
            delta_tilde = res$delta_tilde %||% NA_real_, tau = res$tau %||% tau_fixed,
            K0 = res$K0 %||% NA_integer_, K1 = res$K1 %||% NA_integer_,
            rep = b, seed_used = this_seed, stringsAsFactors = FALSE)
        }, error = function(e) NULL)

        if (!is.null(one_res)) sum_list[[step]] <- one_res
        utils::setTxtProgressBar(pb, step)
      }
    }
  }

  close(pb)
  message(">>> Grid completed.")
  out <- dplyr::bind_rows(sum_list)
  attr(out, "alpha") <- alpha
  attr(out, "n_sim") <- n_sim
  class(out) <- c("powerLOR", class(out))
  out
}

# ============================================================================
# Power curves
# ============================================================================

#' @title Plot SCAM-Smoothed SaLFC Power Curve (soft power)
#' @description
#' Smooths the empirical SOFT-power surface (\eqn{\text{power} = \overline{n_{TP}/n_{DE}}})
#' with a shape-constrained additive model (monotone-increasing in K per effect
#' size, \code{bs = "mpi"}). Falls back to monotone isotonic interpolation with a
#' warning if \pkg{scam} is unavailable or there are < 4 distinct K. Each anchor
#' is labelled with the fitted power and the mean false-discovery proportion
#' (FDR), which are reported on independent axes. A dashed reference line marks
#' \code{target_power}.
#' @param summary_df Data frame from \code{evaluatePowerSaLFC}.
#' @param alpha Numeric. Significance level (title only).
#' @param k_spline Integer. Spline basis dimension (capped at #unique K).
#' @param target_power Numeric in (0,1). Dashed reference line.
#' @return A ggplot object.
#' @seealso \code{\link{evaluatePowerSaLFC}}
#' @import ggplot2
#' @importFrom dplyr group_by summarise mutate
#' @export
plotPowerCurveSaLFC <- function(summary_df, alpha = 0.05, k_spline = 5,
                                target_power = 0.8) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Package 'ggplot2' required.")
  if (!requireNamespace("mgcv",    quietly = TRUE)) stop("Package 'mgcv' (or 'scam') required.")
  `%>%` <- dplyr::`%>%`

  # 1. Aggregate empirical data  -- SOFT power = mean(power_raw) = mean(n_TP/n_DE)
  sum_df <- summary_df %>%
    dplyr::group_by(.data$K_syn, .data$DE_lfc, .data$adjust_method, .data$test_method) %>%
    dplyr::summarise(mean_power = mean(.data$power_raw, na.rm = TRUE),
                     mean_fdp   = mean(.data$fdp,       na.rm = TRUE), .groups = "drop") %>%
    dplyr::filter(is.finite(.data$mean_power)) %>%
    dplyr::mutate(lab_fdp = ifelse(is.finite(.data$mean_fdp),
                                   sprintf("%.2f", .data$mean_fdp), "NA"))
  sum_df$DE_lfc_fac <- factor(sum_df$DE_lfc)

  # 2-4. Monotone-in-K curve per effect size, then prediction grid + anchors.
  Kvals    <- sort(unique(sum_df$K_syn))
  Kgrid    <- seq(min(Kvals), max(Kvals), length.out = 100)
  use_scam <- requireNamespace("scam", quietly = TRUE) && length(Kvals) >= 4

  if (use_scam) {
    # shape-constrained additive model: monotone-increasing P-splines (bs="mpi")
    s <- mgcv::s
    k_use <- min(k_spline, length(Kvals))
    fml <- stats::as.formula(sprintf(
      "mean_power ~ DE_lfc_fac + s(K_syn, k = %d, bs = \"mpi\", by = DE_lfc_fac)",
      k_use))
    environment(fml) <- environment()
    fit <- scam::scam(fml, data = sum_df)
    pred_grid <- expand.grid(K_syn = Kgrid, DE_lfc_fac = levels(sum_df$DE_lfc_fac))
    pred_grid$smooth_power <- pmax(0, pmin(1, as.numeric(stats::predict(fit, pred_grid))))
    sum_df$fitted_power    <- pmax(0, pmin(1, as.numeric(stats::predict(fit, sum_df))))
  } else {
    if (!requireNamespace("scam", quietly = TRUE))
      warning("scam not installed; using monotone isotonic interpolation.")
    else
      message("Fewer than 4 distinct K values; using monotone isotonic interpolation instead of SCAM.")
    # per-level isotonic (nondecreasing) fit on the empirical means
    .iso <- function(d) { d <- d[order(d$K_syn), ]
    d$fitted_power <- pmax(0, pmin(1, stats::isoreg(d$K_syn, d$mean_power)$yf)); d }
    sum_df    <- do.call(rbind, lapply(split(sum_df, sum_df$DE_lfc_fac), .iso))
    pred_grid <- do.call(rbind, lapply(split(sum_df, sum_df$DE_lfc_fac), function(d) {
      d <- d[order(d$K_syn), ]
      data.frame(K_syn = Kgrid, DE_lfc_fac = d$DE_lfc_fac[1],
                 smooth_power = pmax(0, pmin(1,
                                             stats::approx(d$K_syn, d$fitted_power, xout = Kgrid, rule = 2)$y))) }))
  }
  pred_grid$DE_lfc <- as.numeric(as.character(pred_grid$DE_lfc_fac))
  sum_df$lab <- sprintf("Pwr=%.2f\nFDR=%s", sum_df$fitted_power, sum_df$lab_fdp)

  # 5. Plot
  p <- ggplot2::ggplot() +
    ggplot2::geom_hline(yintercept = target_power, linetype = "dashed", colour = "grey60") +
    ggplot2::geom_line(data = pred_grid,
                       ggplot2::aes(x = .data$K_syn, y = .data$smooth_power,
                                    color = factor(.data$DE_lfc), group = factor(.data$DE_lfc)),
                       linewidth = 1) +
    ggplot2::geom_point(data = sum_df,
                        ggplot2::aes(x = .data$K_syn, y = .data$fitted_power,
                                     color = factor(.data$DE_lfc)), size = 2.5)
  if (requireNamespace("ggrepel", quietly = TRUE))
    p <- p + ggrepel::geom_text_repel(data = sum_df,
                                      ggplot2::aes(x = .data$K_syn, y = .data$fitted_power,
                                                   color = factor(.data$DE_lfc), label = .data$lab),
                                      size = 3.5, direction = "y", box.padding = 0.5, point.padding = 0.3,
                                      min.segment.length = 0, max.overlaps = Inf, show.legend = FALSE)
  p +
    ggplot2::scale_x_continuous(breaks = sort(unique(sum_df$K_syn))) +
    ggplot2::scale_y_continuous(breaks = seq(0, 1, 0.2),
                                expand = ggplot2::expansion(mult = c(0.02, 0.10))) +
    ggplot2::coord_cartesian(ylim = c(0, 1.05)) +
    ggplot2::labs(
      title = "SaLFC Test Power Curve",
      x = "Number of synthetic samples per group (K)",
      y = "Statistical power", color = "Effect Size (LFC)") +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(plot.title = ggplot2::element_text(face = "bold"),
                   legend.position = "bottom")
}



#' @title Plot Composition (LOR) Power / Size Curve
#' @description
#' Smooths the empirical rejection rate (= power for effect > 0, = size at
#' effect = 0) against K, one curve per composition effect size, with a dashed
#' reference line at the nominal level \code{alpha}. Uses a shape-constrained
#' additive model (monotone-increasing in K, \code{scam}, \code{bs = "mpi"}) when
#' \pkg{scam} is installed and there are >= 4 distinct K; otherwise falls back to
#' monotone isotonic interpolation.
#' @param summary_df Data frame from \code{evaluatePowerLOR}.
#' @param alpha Numeric. Nominal level (reference line + title).
#' @param k_spline Integer. Spline basis dimension (capped at #unique K).
#' @return A ggplot object.
#' @seealso \code{\link{evaluatePowerLOR}}
#' @import ggplot2
#' @importFrom dplyr group_by summarise
#' @export
#'

plotPowerCurveLOR <- function(summary_df, alpha = 0.05, k_spline = 5) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Package 'ggplot2' required.")
  `%>%` <- dplyr::`%>%`

  eff_lab <- if (identical(unique(summary_df$effect_type)[1], "target_prop_case"))
    "Target prop." else "Effect (delta pp)"

  sum_df <- summary_df %>%
    dplyr::group_by(.data$K_syn, .data$effect, .data$test_method, .data$null_type) %>%
    dplyr::summarise(reject_rate = mean(.data$reject, na.rm = TRUE),
                     n_eff = sum(is.finite(.data$p_value)), .groups = "drop") %>%
    dplyr::filter(is.finite(.data$reject_rate))
  sum_df$effect_fac <- factor(sum_df$effect)

  # monotone-in-K fit per effect; SCAM if feasible, else isotonic interpolation
  Kvals    <- sort(unique(sum_df$K_syn))
  Kgrid    <- seq(min(Kvals), max(Kvals), length.out = 100)
  use_scam <- requireNamespace("scam", quietly = TRUE) && length(Kvals) >= 4

  if (use_scam) {
    s <- mgcv::s
    k_use <- min(k_spline, length(Kvals))
    fml <- stats::as.formula(sprintf(
      "reject_rate ~ effect_fac + s(K_syn, k = %d, bs = \"mpi\", by = effect_fac)", k_use))
    environment(fml) <- environment()
    fit <- scam::scam(fml, data = sum_df)
    pred_grid <- expand.grid(K_syn = Kgrid, effect_fac = levels(sum_df$effect_fac))
    pred_grid$smooth_rate <- pmax(0, pmin(1, as.numeric(stats::predict(fit, pred_grid))))
    sum_df$fitted_rate    <- pmax(0, pmin(1, as.numeric(stats::predict(fit, sum_df))))
  } else {
    if (!requireNamespace("scam", quietly = TRUE))
      warning("scam not installed; using monotone isotonic interpolation.")
    else
      message("Fewer than 4 distinct K values; using monotone isotonic interpolation instead of SCAM.")
    .iso <- function(d) { d <- d[order(d$K_syn), ]
    d$fitted_rate <- pmax(0, pmin(1, stats::isoreg(d$K_syn, d$reject_rate)$yf)); d }
    sum_df    <- do.call(rbind, lapply(split(sum_df, sum_df$effect_fac), .iso))
    pred_grid <- do.call(rbind, lapply(split(sum_df, sum_df$effect_fac), function(d) {
      d <- d[order(d$K_syn), ]
      data.frame(K_syn = Kgrid, effect_fac = d$effect_fac[1],
                 smooth_rate = pmax(0, pmin(1,
                                            stats::approx(d$K_syn, d$fitted_rate, xout = Kgrid, rule = 2)$y))) }))
  }
  pred_grid$effect <- as.numeric(as.character(pred_grid$effect_fac))
  sum_df$lab <- sprintf("%.2f", sum_df$fitted_rate)

  p <- ggplot2::ggplot() +
    ggplot2::geom_line(data = pred_grid,
                       ggplot2::aes(x = .data$K_syn, y = .data$smooth_rate,
                                    color = factor(.data$effect), group = factor(.data$effect)),
                       linewidth = 1) +
    ggplot2::geom_point(data = sum_df,
                        ggplot2::aes(x = .data$K_syn, y = .data$fitted_rate,
                                     color = factor(.data$effect)), size = 2.5)
  if (requireNamespace("ggrepel", quietly = TRUE))
    p <- p + ggrepel::geom_text_repel(data = sum_df,
                                      ggplot2::aes(x = .data$K_syn, y = .data$fitted_rate,
                                                   color = factor(.data$effect), label = .data$lab),
                                      size = 3.3, direction = "y", box.padding = 0.4, point.padding = 0.3,
                                      min.segment.length = 0, max.overlaps = Inf, show.legend = FALSE)
  p +
    ggplot2::scale_x_continuous(breaks = Kvals) +
    ggplot2::scale_y_continuous(breaks = seq(0, 1, 0.2),
                                expand = ggplot2::expansion(mult = c(0.02, 0.10))) +
    ggplot2::coord_cartesian(ylim = c(0, 1.05)) +
    ggplot2::labs(
      title = bquote(paste("LOR Composition Power / Size Curve (", alpha, " = ", .(alpha), ")")),
      x = "Number of synthetic samples per group (K)",
      y = "Statistical Power", color = eff_lab) +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(plot.title = ggplot2::element_text(face = "bold"),
                   legend.position = "bottom")
}


# ============================================================================
# Class taggers
# ----------------------------------------------------------------------------
# evaluatePowerSaLFC() / evaluatePowerLOR() already tag their output, so these
# are only needed to (re)tag a plain data.frame, e.g. one round-tripped through
# rbind/subset. `alpha` is stored for display only.
# ============================================================================

#' @title Tag a SaLFC power grid as class "powerSaLFC"
#' @param x Data frame from \code{evaluatePowerSaLFC}.
#' @param alpha Numeric. The (FDR) level used, stored for display.
#' @return \code{x} with class \code{c("powerSaLFC", ...)} and an \code{alpha} attribute.
#' @export
as_powerSaLFC <- function(x, alpha = 0.05) {
  req <- c("K_syn", "DE_lfc", "power_strict", "fdp", "rep")
  miss <- setdiff(req, names(x))
  if (length(miss)) stop("Not a SaLFC power result; missing columns: ",
                         paste(miss, collapse = ", "))
  attr(x, "alpha") <- alpha; attr(x, "endpoint") <- "SaLFC"
  class(x) <- c("powerSaLFC", setdiff(class(x), "powerSaLFC")); x
}

#' @title Tag a LOR power grid as class "powerLOR"
#' @param x Data frame from \code{evaluatePowerLOR}.
#' @param alpha Numeric. The level used, stored for display.
#' @return \code{x} with class \code{c("powerLOR", ...)} and an \code{alpha} attribute.
#' @export
as_powerLOR <- function(x, alpha = 0.05) {
  req <- c("K_syn", "effect", "reject", "rep")
  miss <- setdiff(req, names(x))
  if (length(miss)) stop("Not a LOR power result; missing columns: ",
                         paste(miss, collapse = ", "))
  attr(x, "alpha") <- alpha; attr(x, "endpoint") <- "LOR"
  class(x) <- c("powerLOR", setdiff(class(x), "powerLOR")); x
}

# ---- internal helpers ------------------------------------------------------
.as_df <- function(x) { d <- as.data.frame(x, stringsAsFactors = FALSE)
class(d) <- "data.frame"; d }

# smallest K reaching power >= target (grid + monotone-interpolated continuous)
.min_k_for_power <- function(Kv, pv, target) {
  o <- order(Kv); Kv <- Kv[o]; pv <- pv[o]
  reached <- any(is.finite(pv) & pv >= target)
  i1 <- if (reached) which(pv >= target)[1] else NA_integer_
  pm <- tryCatch(stats::isoreg(Kv, pv)$yf, error = function(e) pv)  # monotone envelope
  k_cont <- NA_real_
  if (any(is.finite(pm) & pm >= target)) {
    j <- which(pm >= target)[1]
    k_cont <- if (j == 1L) Kv[1]
    else { y0 <- pm[j - 1]; y1 <- pm[j]
    if (is.finite(y1) && is.finite(y0) && y1 > y0)
      Kv[j - 1] + (target - y0) / (y1 - y0) * (Kv[j] - Kv[j - 1])
    else Kv[j] }
  }
  list(reached = reached,
       k_grid = if (reached) Kv[i1] else NA_real_,
       p_at   = if (reached) pv[i1] else NA_real_,
       k_cont = k_cont,
       max_power = max(pv, na.rm = TRUE),
       k_at_max  = Kv[which.max(pv)])
}

# build K x effect numeric matrices (value + se) from a per-cell aggregate
.cell_matrices <- function(tab, effcol, valcol, secol) {
  Kv <- sort(unique(tab$K_syn)); Ev <- sort(unique(tab[[effcol]]))
  vm <- sm <- matrix(NA_real_, length(Kv), length(Ev),
                     dimnames = list(paste0("K=", Kv), format(Ev, trim = TRUE)))
  for (r in seq_len(nrow(tab))) {
    i <- match(tab$K_syn[r], Kv); j <- match(tab[[effcol]][r], Ev)
    vm[i, j] <- tab[[valcol]][r]; sm[i, j] <- tab[[secol]][r]
  }
  list(value = vm, se = sm, Kv = Kv, Ev = Ev)
}

# print a K x effect rate matrix with a '*' marker for cells reaching target
.print_rate_matrix <- function(vm, target) {
  cm <- matrix("", nrow(vm), ncol(vm), dimnames = dimnames(vm))
  for (i in seq_len(nrow(vm))) for (j in seq_len(ncol(vm))) {
    v <- vm[i, j]
    cm[i, j] <- if (!is.finite(v)) "   -- "
    else sprintf("%.2f%s", v, if (v >= target) "*" else " ")
  }
  print(noquote(cm))
}


# ============================================================================
# SaLFC (expression endpoint)
# ============================================================================

#' @title Print a SaLFC power grid (compact, soft power)
#' @param x A \code{powerSaLFC} object.
#' @param digits Rounding for the power matrix.
#' @param ... Ignored.
#' @return \code{x}, invisibly.
#' @export
print.powerSaLFC <- function(x, digits = 3, ...) {
  df <- .as_df(x); Kv <- sort(unique(df$K_syn)); Lv <- sort(unique(df$DE_lfc))
  nsim <- attr(x, "n_sim") %||% suppressWarnings(max(df$rep))
  failed <- max(0L, length(Kv) * length(Lv) * nsim - nrow(df))
  cat("spaCraft SaLFC power result\n")
  cat(sprintf("  test %s | adjust %s | eval %s | spatial %s%s\n",
              df$test_method[1] %||% NA, df$adjust_method[1] %||% NA,
              df$eval_domain_col[1] %||% NA, df$spatial_mode[1] %||% NA,
              if (isTRUE(df$rearrange[1])) " (rearranged)" else ""))
  cat(sprintf("  K in {%s} | DE_lfc in {%s} | %s reps/cell | %d rows%s\n",
              paste(Kv, collapse = ", "), paste(Lv, collapse = ", "), format(nsim),
              nrow(df), if (failed > 0) sprintf(" (%d failed)", failed) else ""))
  m <- tapply(df$power_raw, list(K = df$K_syn, DE_lfc = df$DE_lfc), mean)
  cat("\nmean soft power (TP rate) [rows = K, cols = DE_lfc]:\n"); print(round(m, digits))
  cat("\nUse summary() for the full report (SE, K-to-power, FDR control).\n")
  invisible(x)
}


#' @title Summarize a SaLFC power analysis (lm-style report, soft power)
#' @description Produces an interpretive summary on the SOFT-power metric
#'   (\code{power = mean(power_raw) = mean(n_TP/n_DE)}): a power table over
#'   K x DE_lfc, the minimum K reaching \code{target_power} per effect (grid value
#'   plus a monotone-interpolated continuous estimate), and a separate FDR-control
#'   diagnostic. The all-or-nothing \code{power_strict} is retained in
#'   \code{$by_cell} for reference only.
#' @param object A \code{powerSaLFC} object.
#' @param target_power Numeric in (0,1).
#' @param ... Ignored.
#' @return An object of class \code{summary.powerSaLFC} (a list with
#'   \code{$by_cell}, \code{$power}, \code{$se}, \code{$min_k}, \code{$fdr},
#'   \code{$design}); printed as a report. \code{$power} is the soft-power matrix.
#' @importFrom stats sd isoreg median
#' @export
summary.powerSaLFC <- function(object, target_power = 0.8, ...) {
  df  <- .as_df(object)
  key <- interaction(df$K_syn, df$DE_lfc, drop = TRUE, lex.order = TRUE)
  tab <- do.call(rbind, lapply(split(seq_len(nrow(df)), key), function(ix) {
    d <- df[ix, , drop = FALSE]; n <- nrow(d); pr <- d$power_raw   # SOFT power
    data.frame(DE_lfc = d$DE_lfc[1], K_syn = d$K_syn[1], n_rep = n,
               power = mean(pr, na.rm = TRUE),
               power_se = if (n > 1) stats::sd(pr, na.rm = TRUE) / sqrt(n) else NA_real_,
               power_strict = mean(d$power_strict, na.rm = TRUE),   # reference only
               fdp = mean(d$fdp, na.rm = TRUE), pct_fp_free = mean(d$strict_ok),
               mean_rej = mean(d$n_rej), mean_TP = mean(d$n_TP),
               mean_FP = mean(d$n_FP), stringsAsFactors = FALSE)
  })); rownames(tab) <- NULL; tab <- tab[order(tab$DE_lfc, tab$K_syn), ]

  mats <- .cell_matrices(tab, "DE_lfc", "power", "power_se")
  min_k <- do.call(rbind, lapply(mats$Ev, function(L) {
    s <- tab[tab$DE_lfc == L, ]; mk <- .min_k_for_power(s$K_syn, s$power, target_power)
    data.frame(DE_lfc = L, reached = mk$reached, K_grid = mk$k_grid,
               power_at_K = mk$p_at, K_continuous = mk$k_cont,
               max_power = mk$max_power, K_at_max = mk$k_at_max,
               stringsAsFactors = FALSE)
  }))

  design <- list(
    test = if (identical(df$test_method[1], "t")) "Welch t" else df$test_method[1],
    adjust = df$adjust_method[1] %||% NA, eval = df$eval_domain_col[1] %||% NA,
    spatial = df$spatial_mode[1] %||% NA, rearrange = isTRUE(df$rearrange[1]),
    n_sim = attr(object, "n_sim") %||% suppressWarnings(max(df$rep)),
    n_cells = nrow(tab), n_K = length(mats$Kv), n_eff = length(mats$Ev),
    n_total = nrow(df), alpha = attr(object, "alpha") %||% NA,
    Kv = mats$Kv, Ev = mats$Ev)
  failed <- max(0L, design$n_K * design$n_eff * design$n_sim - design$n_total)

  fdr <- list(mean_fdp = mean(df$fdp, na.rm = TRUE), max_fdp = max(df$fdp, na.rm = TRUE),
              pct_fp_free = mean(df$strict_ok, na.rm = TRUE))

  structure(list(by_cell = tab, power = mats$value, se = mats$se, min_k = min_k,
                 fdr = fdr, design = design, target_power = target_power,
                 failed = failed),
            class = "summary.powerSaLFC")
}

#' @method print summary.powerSaLFC
#' @export
print.summary.powerSaLFC <- function(x, ...) {
  d <- x$design
  cat("SaLFC power analysis\n"); cat(strrep("=", 20), "\n\n", sep = "")
  cat("Design:\n")
  cat("  Endpoint   : SaLFC (expression TREAT; soft power = mean TP rate)\n")
  cat(sprintf("  Test       : %s  |  adjust: %s  |  eval: %s  |  spatial: %s%s\n",
              d$test, d$adjust, d$eval, d$spatial, if (d$rearrange) " (rearranged)" else ""))
  cat(sprintf("  Replicates : %s per cell  |  %d K x %d effects = %d cells  |  %d sims%s\n",
              format(d$n_sim), d$n_K, d$n_eff, d$n_cells, d$n_total,
              if (x$failed > 0) sprintf(" (%d failed)", x$failed) else ""))
  cat(sprintf("  alpha (FDR): %s\n\n", format(d$alpha)))

  cat(sprintf("Soft power by K (rows) and DE_lfc (cols)   [* = power >= %.2f]:\n", x$target_power))
  .print_rate_matrix(x$power, x$target_power)
  se <- x$se[is.finite(x$se)]
  if (length(se)) cat(sprintf("  Monte-Carlo SE: median %.3f, max %.3f\n",
                              stats::median(se), max(se)))

  cat(sprintf("\nMinimum K for soft power >= %.2f:\n", x$target_power))
  mk <- x$min_k
  for (i in seq_len(nrow(mk))) {
    if (isTRUE(mk$reached[i]))
      cat(sprintf("  DE_lfc = %-5s:  K = %g  (power %.2f%s)\n",
                  format(mk$DE_lfc[i]), mk$K_grid[i], mk$power_at_K[i],
                  if (is.finite(mk$K_continuous[i]))
                    sprintf("; interp %.1f", mk$K_continuous[i]) else ""))
    else
      cat(sprintf("  DE_lfc = %-5s:  not reached  (max power %.2f at K=%g)\n",
                  format(mk$DE_lfc[i]), mk$max_power[i], mk$K_at_max[i]))
  }

  cat(sprintf("\nFDR control (separate axis): mean FDP = %.3f (nominal %s) | max FDP = %.3f | FP-free reps: %.0f%%\n",
              x$fdr$mean_fdp, format(d$alpha), x$fdr$max_fdp, 100 * x$fdr$pct_fp_free))
  cat("---\npower = soft power (mean n_TP/n_DE), reported independently of FDR; '*' marks cells at/above target; interp = continuous K (monotone fit).\n")
  invisible(x)
}


# ============================================================================
# LOR (composition endpoint)
# ============================================================================

#' @title Print a LOR composition power grid (compact)
#' @param x A \code{powerLOR} object.
#' @param digits Rounding for the rate matrix.
#' @param ... Ignored.
#' @return \code{x}, invisibly.
#' @export
print.powerLOR <- function(x, digits = 3, ...) {
  df <- .as_df(x); Kv <- sort(unique(df$K_syn)); Ev <- sort(unique(df$effect))
  etype <- df$effect_type[1] %||% NA
  ecol  <- if (identical(etype, "target_prop_case")) "target_prop_case" else "delta_pp"
  nsim <- attr(x, "n_sim") %||% suppressWarnings(max(df$rep))
  failed <- max(0L, length(Kv) * length(Ev) * nsim - nrow(df))
  cat("spaCraft LOR composition power result\n")
  cat(sprintf("  test %s | null %s | eval %s\n",
              df$test_method[1] %||% NA, df$null_type[1] %||% NA, df$eval_on[1] %||% NA))
  cat(sprintf("  K in {%s} | %s in {%s} | %s reps/cell | %d rows%s\n",
              paste(Kv, collapse = ", "), ecol, paste(Ev, collapse = ", "), format(nsim),
              nrow(df), if (failed > 0) sprintf(" (%d failed)", failed) else ""))
  m <- tapply(df$reject, list(K = df$K_syn, effect = df$effect), mean)
  cat(sprintf("\nrejection rate [rows = K, cols = %s]:\n", ecol)); print(round(m, digits))
  if (identical(etype, "delta_pp") && 0 %in% Ev)
    cat("  (column effect = 0 is the empirical size)\n")
  cat("\nUse summary() for the full report (SE, K-to-power, size check).\n")
  invisible(x)
}

#' @title Summarize a LOR composition power analysis (lm-style report)
#' @description Produces an interpretive summary: a rejection-rate table over
#'   K x effect, the minimum K reaching \code{target_power} per (non-null)
#'   effect (grid + continuous estimate), and the empirical size at the null.
#' @param object A \code{powerLOR} object.
#' @param target_power Numeric in (0,1).
#' @param ... Ignored.
#' @return An object of class \code{summary.powerLOR}.
#' @importFrom stats isoreg median
#' @export
summary.powerLOR <- function(object, target_power = 0.8, ...) {
  df  <- .as_df(object)
  key <- interaction(df$K_syn, df$effect, drop = TRUE, lex.order = TRUE)
  tab <- do.call(rbind, lapply(split(seq_len(nrow(df)), key), function(ix) {
    d <- df[ix, , drop = FALSE]; n <- nrow(d)
    neff <- sum(is.finite(d$p_value)); r <- mean(d$reject, na.rm = TRUE)
    data.frame(effect = d$effect[1], K_syn = d$K_syn[1], n_rep = n, n_eff = neff,
               reject_rate = r,
               rate_se = if (neff > 0) sqrt(r * (1 - r) / neff) else NA_real_,
               mean_p = mean(d$p_value, na.rm = TRUE),
               mean_delta_hat = mean(d$delta_hat, na.rm = TRUE),
               mean_delta_tilde = mean(d$delta_tilde, na.rm = TRUE),
               stringsAsFactors = FALSE)
  })); rownames(tab) <- NULL; tab <- tab[order(tab$effect, tab$K_syn), ]

  etype       <- df$effect_type[1] %||% NA
  null_effect <- if (identical(etype, "delta_pp")) 0 else NA_real_

  mats <- .cell_matrices(tab, "effect", "reject_rate", "rate_se")
  eff_eval <- setdiff(mats$Ev, null_effect)
  min_k <- if (length(eff_eval)) do.call(rbind, lapply(eff_eval, function(E) {
    s <- tab[tab$effect == E, ]; mk <- .min_k_for_power(s$K_syn, s$reject_rate, target_power)
    data.frame(effect = E, reached = mk$reached, K_grid = mk$k_grid,
               power_at_K = mk$p_at, K_continuous = mk$k_cont,
               max_power = mk$max_power, K_at_max = mk$k_at_max,
               stringsAsFactors = FALSE)
  })) else NULL

  size <- if (!is.na(null_effect) && null_effect %in% mats$Ev)
    tab[tab$effect == null_effect, c("K_syn", "reject_rate", "rate_se")] else NULL

  design <- list(
    test = if (identical(df$test_method[1], "t")) "Welch t" else df$test_method[1],
    null = df$null_type[1] %||% NA, eval = df$eval_on[1] %||% NA, effect_type = etype,
    tau = suppressWarnings(stats::median(df$tau, na.rm = TRUE)),
    delta0 = suppressWarnings(stats::median(df$delta0, na.rm = TRUE)),
    n_sim = attr(object, "n_sim") %||% suppressWarnings(max(df$rep)),
    n_cells = nrow(tab), n_K = length(mats$Kv), n_eff = length(mats$Ev),
    n_total = nrow(df), alpha = attr(object, "alpha") %||% NA,
    ecol = if (identical(etype, "target_prop_case")) "target_prop_case" else "delta_pp",
    Kv = mats$Kv, Ev = mats$Ev)
  failed <- max(0L, design$n_K * design$n_eff * design$n_sim - design$n_total)

  structure(list(by_cell = tab, reject_rate = mats$value, se = mats$se,
                 min_k = min_k, size = size, design = design,
                 target_power = target_power, null_effect = null_effect,
                 failed = failed),
            class = "summary.powerLOR")
}

#' @method print summary.powerLOR
#' @export
print.summary.powerLOR <- function(x, ...) {
  d <- x$design
  cat("LOR composition power analysis\n"); cat(strrep("=", 30), "\n\n", sep = "")
  cat("Design:\n")
  cat("  Endpoint   : LOR (composition log-odds-ratio; rejection rate)\n")
  cat(sprintf("  Test       : %s  |  null: %s  |  eval: %s\n", d$test, d$null, d$eval))
  cat(sprintf("  Effect     : %s\n", d$effect_type))
  if (identical(d$null, "treat") && is.finite(d$tau))
    cat(sprintf("  TREAT      : margin tau = %.4f | pilot baseline delta0 (beta1) = %.4f\n",
                d$tau, d$delta0))
  cat(sprintf("  Replicates : %s per cell  |  %d K x %d effects = %d cells  |  %d sims%s\n",
              format(d$n_sim), d$n_K, d$n_eff, d$n_cells, d$n_total,
              if (x$failed > 0) sprintf(" (%d failed)", x$failed) else ""))
  cat(sprintf("  alpha      : %s\n\n", format(d$alpha)))

  cat(sprintf("Rejection rate by K (rows) and %s (cols)   [* = power >= %.2f]:\n",
              d$ecol, x$target_power))
  .print_rate_matrix(x$reject_rate, x$target_power)
  se <- x$se[is.finite(x$se)]
  if (length(se)) cat(sprintf("  Monte-Carlo (binomial) SE: median %.3f, max %.3f\n",
                              stats::median(se), max(se)))

  cat("\nEmpirical size (null): ")
  if (!is.null(x$size)) {
    rr <- x$size$reject_rate
    cat(sprintf("%s = %g -> rate %.3f-%.3f across K (nominal %s)\n",
                d$ecol, x$null_effect, min(rr, na.rm = TRUE), max(rr, na.rm = TRUE),
                format(d$alpha)))
  } else if (identical(d$effect_type, "target_prop_case")) {
    cat("not in grid (add target_prop_case = p0 = plogis(beta0) to evaluate size)\n")
  } else cat("not in grid (add delta_pp = 0 to evaluate size)\n")

  cat(sprintf("\nMinimum K for power >= %.2f:\n", x$target_power))
  mk <- x$min_k
  if (is.null(mk) || nrow(mk) == 0) cat("  (no non-null effect in grid)\n") else
    for (i in seq_len(nrow(mk))) {
      if (isTRUE(mk$reached[i]))
        cat(sprintf("  %s = %-6s:  K = %g  (rate %.2f%s)\n", d$ecol,
                    format(mk$effect[i]), mk$K_grid[i], mk$power_at_K[i],
                    if (is.finite(mk$K_continuous[i]))
                      sprintf("; interp %.1f", mk$K_continuous[i]) else ""))
      else
        cat(sprintf("  %s = %-6s:  not reached  (max rate %.2f at K=%g)\n", d$ecol,
                    format(mk$effect[i]), mk$max_power[i], mk$K_at_max[i]))
    }
  cat("---\nrate = P(reject) = power for non-null effects, size at the null; '*' marks cells at/above target.\n")
  invisible(x)
}
