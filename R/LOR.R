#' @title Composition Test (K-driven, Sample-Level): Naive Two-Sample or TREAT
#' @description
#' Computes per-sample \eqn{L_k = \log((y_k + \epsilon)/(b_k + \epsilon))} and
#' compares Case vs Control at the sample level using a Welch t-test or a
#' Wilcoxon rank-sum test, for a single (target, reference) domain pair. Here
#' \eqn{y_k} and \eqn{b_k} are the spot counts of the target and reference
#' domains in sample \eqn{k}.
#'
#' Two null modes, controlled by \code{null_type}:
#' \describe{
#'   \item{\code{"treat"} (default)}{TREAT-style test (Smyth, 2009) centered on
#'     the pilot baseline \eqn{\delta_0 = \beta_1}.
#'     \itemize{
#'       \item \code{tau = 0} (centered null): \eqn{H_0:\ \hat\delta = \delta_0}.
#'       \item \code{tau > 0} (TREAT): \eqn{H_0:\ |\hat\delta - \delta_0| \le \tau}.
#'     }}
#'   \item{\code{"two_sample"}}{Naive two-sample test; the pilot baseline is
#'     ignored and \eqn{\delta_0 = 0}, i.e. \eqn{H_0:\ \bar\delta_1 = \bar\delta_0}.}
#' }
#'
#' The threshold \code{tau} may be a numeric scalar or the string \code{"pilot"},
#' in which case it is estimated from \code{object@pilot_data} via
#' \code{estimate_lor_tau_from_pilot()} for this single (target, reference) pair.
#'
#' @param processed_data List of samples. Each element is a list with
#'   \code{coords} (containing \code{hat_d} or \code{domain}) and \code{meta}
#'   (with \code{group} in \{0, 1\} and optionally \code{sample_id}).
#' @param object A spaCraft object with \code{@params_composition} (and
#'   \code{@pilot_data} when \code{tau = "pilot"}).
#' @param target_domain Character. Target domain label.
#' @param reference_domain Character. Reference domain label.
#' @param test_method "t" (Welch, default) or "wilcoxon".
#' @param alternative "greater" (default), "two.sided", or "less".
#' @param null_type "treat" (default) or "two_sample".
#' @param tau Numeric scalar >= 0, or the string \code{"pilot"}. TREAT threshold
#'   (used only when \code{null_type = "treat"}); ignored (forced to 0) for
#'   "two_sample". If \code{"pilot"}, the threshold is estimated on the fly.
#' @param pilot_tau_args Named list of extra arguments forwarded to
#'   \code{estimate_lor_tau_from_pilot()} when \code{tau = "pilot"}, e.g.
#'   \code{list(target_prop_min = 0.12, use_noise_floor = FALSE)} to anchor
#'   \eqn{\tau} to a minimum meaningful CASE target proportion from the fitted
#'   composition baseline. \code{object}, \code{target_domain},
#'   \code{reference_domain}, \code{null_type}, \code{test_method},
#'   \code{alternative}, and \code{eps} are supplied automatically.
#' @param eps Numeric > 0. Continuity correction (default 0.5).
#' @param verbose Logical (default TRUE).
#'
#' @return A named list with elements \code{null_type, target, reference, K0, K1,
#'   delta_hat, delta0, delta_tilde, tau, eps, se, Z_bar0, Z_bar1, y, m, sample_id,
#'   group, stat, df, p_value, method} (and \code{stat_treat} in TREAT mode).
#'   \code{se} is the Welch standard error of \eqn{\hat\delta} ("t") or \code{NA}
#'   ("wilcoxon").
#'
#' @references Smyth, G.K. (2009). Testing significance relative to a fold-change
#'   threshold is a TREAT. \emph{Bioinformatics}, 25(6), 765-771.
#'
#' @importFrom stats wilcox.test var pt
#' @export
LOR <- function(processed_data,
                object,
                target_domain    = "WM",
                reference_domain = "Layer6",
                test_method      = c("t", "wilcoxon"),
                alternative      = c("greater", "two.sided", "less"),
                null_type        = c("treat", "two_sample"),
                tau              = 0,
                pilot_tau_args   = list(),
                eps              = 0.5,
                verbose          = TRUE) {

  test_method <- match.arg(test_method)
  alternative <- match.arg(alternative)
  null_type   <- match.arg(null_type)

  # ── Input validation ────────────────────────────────────────────────────────
  if (!is.list(processed_data) || length(processed_data) == 0)
    stop("processed_data must be a non-empty list.")
  if (is.null(object@params_composition) || length(object@params_composition) == 0)
    stop("object@params_composition is empty.")
  if (!is.finite(eps) || eps <= 0) stop("eps must be > 0.")

  # ── Optional: estimate tau from the pilot on the fly (single pair) ───────────
  # The estimator runs LOR on pilot with a numeric tau = 0, so no recursion.
  if (is.character(tau) && length(tau) == 1L && tolower(tau) == "pilot") {
    if (null_type == "two_sample") {
      # two_sample forces tau = 0 anyway; skip the (pointless) pilot estimation.
      tau <- 0
    } else {
      est_defaults <- list(
        object = object, target_domain = target_domain,
        reference_domain = reference_domain, null_type = null_type,
        test_method = test_method, alternative = alternative,
        eps = eps, verbose = verbose, test_fn = sys.function()
      )
      est_args <- utils::modifyList(est_defaults, pilot_tau_args)
      tau <- do.call(estimate_lor_tau_from_pilot, est_args)
    }
  }

  if (!is.finite(tau) || tau < 0) stop("tau must be >= 0 (or the string \"pilot\").")
  if (null_type == "two_sample" && tau > 0) {
    warning("tau > 0 is ignored when null_type = 'two_sample'. Setting tau = 0.")
    tau <- 0
  }

  # ── Pilot baseline delta0 ────────────────────────────────────────────────────
  # treat      : delta0 = pilot beta1 (null center)
  # two_sample : delta0 = 0           (H0: bar1 = bar0)
  if (null_type == "two_sample") {
    delta0 <- 0
    if (verbose) message(">>> [two_sample] delta0 = 0. Pilot baseline ignored.")
  } else {
    delta0 <- object@params_composition$beta_binomial$beta1
    if (!is.finite(delta0)) {
      warning("delta0 (beta1) is non-finite; defaulting to 0.")
      delta0 <- 0
    }
    if (verbose) message(sprintf(
      ">>> [treat] delta0 = %.4f, tau = %.4f  (%s)",
      delta0, tau, if (tau == 0) "centered null" else "TREAT threshold"))
  }

  # ── Per-sample L_k (vectorized) ──────────────────────────────────────────────
  nm <- names(processed_data)
  if (is.null(nm) || any(nm == "")) nm <- paste0("s", seq_along(processed_data))

  stats_mat <- vapply(seq_along(processed_data), function(i) {
    samp   <- processed_data[[i]]
    na_out <- c(L_k = NA_real_, grp = NA_real_, m_k = NA_real_, y_k = NA_real_)

    if (is.null(samp$coords)) return(na_out)
    preds <- as.character(samp$coords$hat_d %||% samp$coords$domain)
    if (is.null(preds)) return(na_out)

    g <- NA_real_
    if (!is.null(samp$meta) && "group" %in% names(samp$meta))
      g <- as.numeric(samp$meta$group)
    if (!is.finite(g) || !(g %in% c(0, 1))) return(na_out)

    y_k <- sum(preds == target_domain,    na.rm = TRUE)
    b_k <- sum(preds == reference_domain, na.rm = TRUE)
    m_k <- y_k + b_k
    if (!is.finite(m_k) || m_k <= 0) return(na_out)

    c(L_k = log(y_k + eps) - log(b_k + eps), grp = g, m_k = m_k, y_k = y_k)
  }, numeric(4))

  # ── Filter valid samples ─────────────────────────────────────────────────────
  valid <- is.finite(stats_mat["L_k", ]) & is.finite(stats_mat["grp", ])
  if (!any(valid)) return(list(p_value = NA_real_, msg = "No valid samples."))

  Z     <- stats_mat["L_k", valid]
  grp   <- stats_mat["grp", valid]
  m_vec <- stats_mat["m_k", valid]
  y_vec <- stats_mat["y_k", valid]

  sid <- vapply(seq_along(processed_data)[valid], function(i) {
    samp <- processed_data[[i]]
    if (!is.null(samp$meta) && "sample_id" %in% names(samp$meta))
      return(as.character(samp$meta$sample_id))
    nm[i]
  }, character(1))

  if (length(unique(grp)) < 2)
    return(list(p_value = NA_real_, msg = "Only one group present."))

  Z0 <- Z[grp == 0]; Z1 <- Z[grp == 1]
  K0 <- length(Z0);  K1 <- length(Z1)
  if (K0 < 1 || K1 < 1)
    return(list(p_value = NA_real_, msg = "Insufficient samples per group."))

  bar0        <- mean(Z0); bar1 <- mean(Z1)
  delta_hat   <- bar1 - bar0
  delta_tilde <- delta_hat - delta0  # two_sample: delta0 = 0 -> delta_tilde = delta_hat

  out <- list(
    null_type   = null_type,
    target      = target_domain,    reference   = reference_domain,
    K0          = K0,               K1          = K1,
    delta_hat   = delta_hat,        delta0      = delta0,
    delta_tilde = delta_tilde,      tau         = tau,
    eps         = eps,              se          = NA_real_,
    Z_bar0      = bar0,             Z_bar1      = bar1,
    y           = y_vec,            m           = m_vec,
    sample_id   = sid,              group       = grp
  )

  # ===========================================================================
  # (i) Welch t-test
  # ===========================================================================
  if (test_method == "t") {
    s0 <- if (K0 >= 2) stats::var(Z0) else NA_real_
    s1 <- if (K1 >= 2) stats::var(Z1) else NA_real_

    if (!is.finite(s0) || !is.finite(s1)) {
      out$stat <- out$df <- out$p_value <- NA_real_
      out$method <- "welch_t"
      out$msg    <- "Welch t-test requires >= 2 samples per group."
      return(out)
    }

    se <- sqrt(s1 / K1 + s0 / K0)
    df <- (s1/K1 + s0/K0)^2 / ((s1/K1)^2/(K1 - 1) + (s0/K0)^2/(K0 - 1))

    if (!is.finite(se) || se <= 0 || !is.finite(df) || df <= 0) {
      out$stat <- out$df <- out$p_value <- NA_real_
      out$method <- "welch_t"
      out$msg    <- "Non-finite SE or df."
      return(out)
    }

    # t_stat = delta_tilde / se  (treat: (delta_hat-delta0)/se; two_sample: delta_hat/se)
    t_stat <- delta_tilde / se

    if (null_type == "treat" && tau > 0) {
      if (alternative == "two.sided") {
        t_shift    <- (abs(delta_tilde) - tau) / se
        p_val      <- 2 * stats::pt(t_shift, df = df, lower.tail = FALSE)
        out$method <- "welch_t_treat_two_sided"
      } else if (alternative == "greater") {
        t_shift    <- (delta_tilde - tau) / se
        p_val      <- stats::pt(t_shift, df = df, lower.tail = FALSE)
        out$method <- "welch_t_treat_greater"
      } else {
        t_shift    <- (-delta_tilde - tau) / se
        p_val      <- stats::pt(t_shift, df = df, lower.tail = FALSE)
        out$method <- "welch_t_treat_less"
      }
      out$stat_treat <- t_shift
      out$p_value    <- min(1, p_val)
    } else {
      if (alternative == "two.sided") {
        p_val <- 2 * stats::pt(abs(t_stat), df = df, lower.tail = FALSE)
      } else if (alternative == "greater") {
        p_val <- stats::pt(t_stat, df = df, lower.tail = FALSE)
      } else {
        p_val <- stats::pt(t_stat, df = df, lower.tail = TRUE)
      }
      out$p_value <- p_val
      out$method  <- if (null_type == "two_sample")
        paste0("welch_t_two_sample_", alternative)
      else
        paste0("welch_t_treat_centered_", alternative)  # tau = 0
    }

    out$stat <- t_stat
    out$df   <- df
    out$se   <- se

    if (verbose) message(sprintf(
      ">>> Composition K-test (t) [%s]: alt=%s, delta_hat=%.4f, delta0=%.4f, delta_tilde=%.4f, tau=%.4f, p=%.3e",
      null_type, alternative, delta_hat, delta0, delta_tilde, tau, out$p_value))
    return(out)
  }

  # ===========================================================================
  # (ii) Wilcoxon rank-sum
  # ===========================================================================
  if (test_method == "wilcoxon") {
    # mu: two_sample -> 0; treat tau=0 -> delta0; treat tau>0 -> delta0 +/- tau
    mu_base <- delta0

    if (null_type == "treat" && tau > 0) {
      if (alternative == "two.sided") {
        p_one <- if (delta_tilde >= 0) {
          tryCatch(stats::wilcox.test(Z1, Z0, alternative = "greater",
                                      mu = mu_base + tau, exact = FALSE)$p.value,
                   error = function(e) NA_real_)
        } else {
          tryCatch(stats::wilcox.test(Z1, Z0, alternative = "less",
                                      mu = mu_base - tau, exact = FALSE)$p.value,
                   error = function(e) NA_real_)
        }
        p_val      <- if (is.finite(p_one)) min(1, 2 * p_one) else NA_real_
        out$method <- "wilcoxon_treat_two_sided"
      } else if (alternative == "greater") {
        p_val      <- tryCatch(stats::wilcox.test(Z1, Z0, alternative = "greater",
                                                  mu = mu_base + tau, exact = FALSE)$p.value,
                               error = function(e) NA_real_)
        out$method <- "wilcoxon_treat_greater"
      } else {
        p_val      <- tryCatch(stats::wilcox.test(Z1, Z0, alternative = "less",
                                                  mu = mu_base - tau, exact = FALSE)$p.value,
                               error = function(e) NA_real_)
        out$method <- "wilcoxon_treat_less"
      }
    } else {
      p_val <- tryCatch(stats::wilcox.test(Z1, Z0, alternative = alternative,
                                           mu = mu_base, exact = FALSE)$p.value,
                        error = function(e) NA_real_)
      out$method <- if (null_type == "two_sample")
        paste0("wilcoxon_two_sample_", alternative)
      else
        paste0("wilcoxon_treat_centered_", alternative)  # tau = 0
    }

    out$stat    <- NA_real_
    out$p_value <- p_val

    if (verbose) message(sprintf(
      ">>> Composition K-test (wilcoxon) [%s]: alt=%s, delta_hat=%.4f, delta0=%.4f, delta_tilde=%.4f, tau=%.4f, p=%.3e",
      null_type, alternative, delta_hat, delta0, delta_tilde, tau, out$p_value))
    return(out)
  }

  out
}


#' @title Estimate the Composition TREAT Threshold (tau) from Pilot Data
#' @description
#' Plug-in estimator of the scalar TREAT threshold \eqn{\tau} for the
#' sample-level composition test \code{\link{LOR}}, on the log-odds-ratio (LOR)
#' scale. The threshold is the larger of a scientific anchor and a finite-sample
#' noise floor:
#' \deqn{\tau = \max(\tau_{\mathrm{sci}},\; \tau_{\mathrm{noise}},\; \mathrm{min\_tau}).}
#'
#' \strong{Scientific anchor} \eqn{\tau_{\mathrm{sci}}} (from fitted values; no
#' pilot simulation). The direction-aware log-odds distance between a minimum
#' meaningful CASE target proportion \code{target_prop_min} (written \eqn{p^*})
#' and the fitted baseline proportion \eqn{p_0}:
#' \deqn{\tau_{\mathrm{sci}} = |\mathrm{logit}(p^*) - \mathrm{logit}(p_0)|.}
#' Under \code{null_type = "treat"} the baseline is
#' \eqn{p_0 = \mathrm{plogis}(\beta_0 + \beta_1)} (the fitted case mean, stored as
#' \code{mu1_hat}), because the test centers at \eqn{\delta_0 = \beta_1}; under
#' \code{"two_sample"} it is \eqn{p_0 = \mathrm{plogis}(\beta_0)} (stored as
#' \code{mu0_hat}). This term is 0 when \code{target_prop_min = NULL}.
#'
#' \strong{Noise floor} \eqn{\tau_{\mathrm{noise}}} (optional; finite-sample). A
#' pilot band \eqn{t_{q,df}\, se}: \code{LOR} is run once on the pilot (centered
#' null, \code{tau = 0}) purely to read the Welch standard error \code{se} and
#' degrees of freedom \code{df}; the multiplier is \code{qt(q, df)}, or
#' \code{qnorm(q)} when \code{df} is unavailable (e.g. \code{wilcoxon}). Set
#' \code{use_noise_floor = FALSE} to drop this term and obtain \eqn{\tau} as a
#' deterministic function of the fitted values alone (no pilot run).
#' Returns one scalar per (target, reference) pair. Normally invoked
#' automatically by \code{LOR(..., tau = "pilot")} or
#' \code{evaluatePowerLOR(..., tau = "auto")}; can also be called directly. The
#' pilot run uses the ground-truth \code{coords$domain} (pilot data has no
#' \code{hat_d}).
#'
#' @param object A spaCraft object with \code{@params_composition} (and, when
#'   \code{use_noise_floor = TRUE}, \code{@pilot_data}).
#' @param target_domain,reference_domain Character domain labels.
#' @param null_type "treat" (default) or "two_sample". Selects the fitted
#'   baseline for the anchor: \code{mu1_hat} (treat) or \code{mu0_hat}
#'   (two_sample). NB: \code{LOR}/\code{evaluatePowerLOR} force \eqn{\tau = 0}
#'   under \code{"two_sample"}, so in that pipeline \eqn{\tau} matters only for
#'   \code{"treat"}.
#' @param test_method "t" (default) or "wilcoxon". Only "t" yields an se for the
#'   noise floor; "wilcoxon" uses the anchor alone (noise term = 0).
#' @param alternative "greater" (default), "less", or "two.sided". Sets the sign
#'   convention of the anchor: \code{greater} keeps only positive departures of
#'   \code{target_prop_min} above the baseline, \code{less} only negative,
#'   \code{two.sided} the absolute gap.
#' @param target_prop_min Numeric in (0,1) or NULL (default). Minimum meaningful
#'   CASE target proportion driving \eqn{\tau_{\mathrm{sci}}}; NULL disables the
#'   anchor (\eqn{\tau_{\mathrm{sci}} = 0}).
#' @param q Numeric in (0,1). Quantile level for the noise-floor multiplier.
#' @param use_noise_floor Logical (default TRUE). Add the pilot noise floor
#'   \eqn{\tau_{\mathrm{noise}}}. FALSE returns
#'   \eqn{\max(\tau_{\mathrm{sci}}, \mathrm{min\_tau})} with no pilot run.
#' @param min_tau Numeric. Lower bound for \eqn{\tau}.
#' @param cap_effect Numeric or NULL. Optional upper cap on \eqn{\tau}.
#' @param eps Numeric > 0. Continuity correction forwarded to \code{LOR}.
#' @param verbose Logical.
#' @param test_fn Function. The composition test to run on the pilot (used only
#'   when \code{use_noise_floor = TRUE}). If NULL, resolves \code{LOR} by name.
#'   When invoked via \code{LOR(..., tau = "pilot")}, the running test function
#'   is passed automatically.
#'
#' @return A numeric scalar \eqn{\tau} (>= \code{min_tau}).
#'
#' @seealso \code{\link{LOR}}, \code{\link{evaluatePowerLOR}}
#'
#' @importFrom stats qnorm qt plogis qlogis
#' @export
estimate_lor_tau_from_pilot <- function(object,
                                        target_domain    = "WM",
                                        reference_domain = "Layer6",
                                        null_type        = c("treat", "two_sample"),
                                        test_method      = c("t", "wilcoxon"),
                                        alternative      = c("greater", "less", "two.sided"),
                                        target_prop_min  = NULL,
                                        q                = 0.5,
                                        use_noise_floor  = TRUE,
                                        min_tau          = 0,
                                        cap_effect       = NULL,
                                        eps              = 0.5,
                                        verbose          = TRUE,
                                        test_fn          = NULL) {
  null_type   <- match.arg(null_type)
  test_method <- match.arg(test_method)
  alternative <- match.arg(alternative)

  if (!is.finite(q) || q <= 0 || q >= 1) stop("q must be in (0,1).")
  if (!is.finite(min_tau) || min_tau < 0) stop("min_tau must be >= 0.")
  if (!is.null(target_prop_min) &&
      (!is.finite(target_prop_min) || target_prop_min <= 0 || target_prop_min >= 1))
    stop("target_prop_min must be NULL or a proportion in (0, 1).")

  # =========================================================================
  # (A) Scientific anchor from the FITTED composition baseline (no pilot run).
  #     tau is on the log-odds (LOR) scale, matching delta.
  #       treat      : baseline = plogis(beta0 + beta1) = mu1_hat  (delta0 = beta1)
  #       two_sample : baseline = plogis(beta0)         = mu0_hat  (delta0 = 0)
  #     tau_sci = direction-aware log-odds gap of target_prop_min from baseline.
  # =========================================================================
  tau_sci <- 0
  p_base  <- NA_real_
  if (!is.null(target_prop_min)) {
    bb <- object@params_composition$beta_binomial
    if (is.null(bb)) stop("object@params_composition$beta_binomial is missing.")
    b0 <- suppressWarnings(as.numeric(bb$beta0))
    b1 <- suppressWarnings(as.numeric(bb$beta1))
    p_base <- switch(null_type,
                     treat      = if (is.finite(b0) && is.finite(b1)) stats::plogis(b0 + b1)
                     else suppressWarnings(as.numeric(bb$mu1_hat)),
                     two_sample = if (is.finite(b0)) stats::plogis(b0)
                     else suppressWarnings(as.numeric(bb$mu0_hat)))
    if (!is.finite(p_base) || p_base <= 0 || p_base >= 1)
      stop("Cannot resolve the pilot baseline proportion from beta_binomial ",
           "(need finite beta0[/beta1] or mu0_hat/mu1_hat).")
    d_star  <- stats::qlogis(target_prop_min) - stats::qlogis(p_base)  # signed gap
    tau_sci <- switch(alternative,
                      greater   = max( d_star, 0),
                      less      = max(-d_star, 0),
                      two.sided = abs(d_star))
    if (!is.finite(tau_sci)) tau_sci <- 0
  }

  # =========================================================================
  # (B) Finite-sample NOISE FLOOR from the pilot (optional; best-effort).
  # =========================================================================
  tau_noise <- 0
  se <- NA_real_; df_p <- Inf
  if (isTRUE(use_noise_floor)) {
    if (is.null(object@pilot_data) || length(object@pilot_data) == 0)
      stop("use_noise_floor = TRUE requires a non-empty object@pilot_data.")

    if (is.null(test_fn)) {
      if (exists("LOR", mode = "function")) {
        test_fn <- match.fun("LOR")
      } else if (exists("test_composition_K_treat", mode = "function")) {
        test_fn <- match.fun("test_composition_K_treat")
      } else stop("No composition test function found; pass it via `test_fn`.")
    }
    if (!is.function(test_fn)) stop("`test_fn` must be a function.")

    # Build pilot processed_data (ground-truth domains; LOR uses domain via %||%)
    pilot_pd <- lapply(seq_along(object@pilot_data), function(i) {
      s <- object@pilot_data[[i]]
      if (is.null(s$coords) || is.null(s$group)) return(NULL)
      if (!("domain" %in% names(s$coords))) return(NULL)
      sid <- if (!is.null(names(object@pilot_data)) && names(object@pilot_data)[i] != "")
        names(object@pilot_data)[i] else paste0("pilot_", i)
      list(coords = s$coords,
           meta   = data.frame(sample_id = sid, group = as.numeric(s$group),
                               stringsAsFactors = FALSE))
    })
    pilot_pd <- Filter(Negate(is.null), pilot_pd)
    if (length(pilot_pd) == 0) stop("No valid pilot samples with coords/group.")
    names(pilot_pd) <- vapply(pilot_pd, function(x) x$meta$sample_id, character(1))

    # Run the pilot test once (centered null, tau = 0) purely for se and df.
    pilot_res <- tryCatch(
      test_fn(processed_data   = pilot_pd, object = object,
              target_domain    = target_domain, reference_domain = reference_domain,
              test_method      = test_method, alternative = alternative,
              null_type        = "treat", tau = 0, eps = eps, verbose = FALSE),
      error = function(e) { message("\n[ERROR in LOR pilot] ", e$message); NULL })

    if (!is.null(pilot_res)) {
      dt <- pilot_res$delta_tilde %||% NA_real_
      se <- pilot_res$se %||% NA_real_               # explicit se if available
      if ((!is.finite(se)) && !is.null(pilot_res$stat) &&
          is.finite(pilot_res$stat) && pilot_res$stat != 0 && is.finite(dt))
        se <- abs(dt / pilot_res$stat)               # else recover: stat = dt / se
      if (!is.null(pilot_res$df) && is.finite(pilot_res$df)) df_p <- pilot_res$df
    }
    if (!is.finite(se) || se < 0) se <- 0

    t_q <- if (is.finite(df_p)) stats::qt(q, df = df_p) else stats::qnorm(q)
    if (!is.finite(t_q)) t_q <- 0
    t_q <- max(t_q, 0)
    tau_noise <- if (test_method == "t" && se > 0) t_q * se else 0
  }

  # =========================================================================
  # (C) Combine: tau is at least the practical margin AND the noise floor.
  # =========================================================================
  if (is.null(target_prop_min) && !isTRUE(use_noise_floor))
    warning("estimate_lor_tau_from_pilot: target_prop_min = NULL and ",
            "use_noise_floor = FALSE  ->  tau = min_tau (", format(min_tau), ").")

  tau <- max(tau_sci, tau_noise, min_tau, 0)
  if (!is.null(cap_effect) && is.finite(cap_effect) && cap_effect > 0)
    tau <- min(tau, cap_effect)

  if (verbose) {
    se_txt <- if (is.finite(se)) sprintf("%.4f", se) else "NA"
    df_txt <- if (is.finite(df_p)) sprintf("%.1f", df_p) else "NA"
    p_txt  <- if (!is.null(target_prop_min))
      sprintf(" | p*=%.4f (base=%.4f)", target_prop_min, p_base) else ""
    message(sprintf(
      ">>> [pilot tau] %s | tgt=%s ref=%s alt=%s%s | tau_sci=%.4f | tau_noise=%.4f (q=%.2f, se=%s, df=%s) -> tau=%.4f",
      null_type, target_domain, reference_domain, alternative, p_txt,
      tau_sci, tau_noise, q, se_txt, df_txt, tau))
  }

  tau
}
