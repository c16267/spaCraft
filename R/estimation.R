#' @title Estimate Expression Parameters (Multi-Sample Optimized)
#' @description
#' Estimates gene expression parameters across multiple samples.
#' Uses Rcpp/OpenMP to pre-compute biological variance, grand means,
#' and residuals, massively speeding up the spatial estimation loop.
#'
#' @param object A \code{spaCraft} object.
#' @param target_domain Character. Target domain name.
#' @param reference_domain Character. Reference domain name.
#' @param genes_to_use Character vector or NULL. Genes to evaluate.
#' @param n_neighbors Integer. Number of neighbors for spatial modeling (default: 15).
#' @param n_cores Integer. Number of cores for parallel processing.
#' @param verbose Logical. Print progress messages.
#'
#' @return Updated \code{spaCraft} object with \code{params_expression} populated.
#' @importFrom stats var runif
#' @useDynLib spaCraft, .registration = TRUE
#' @export
#'
#' @examples
#' \dontrun{
#' data("mini_obj_features", package = "spaCraft")
#' spaCraft_obj <- estimateExpressionParams(
#'   object = mini_obj_features,
#'   target_domain = "WM",
#'   reference_domain = "Layer6",
#'   genes_to_use = sets$G_svg,
#'   n_cores = 8
#' )
#' }
estimateExpressionParams <- function(object,
                                        target_domain,
                                        reference_domain,
                                        genes_to_use = NULL,
                                        n_neighbors = 15,
                                        n_cores = 1,
                                        verbose = TRUE) {

  # CRAN Safe: Save and restore environment variables
  old_env <- Sys.getenv(c("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS"))
  on.exit({
    for (nm in names(old_env)) {
      if (old_env[nm] == "") Sys.unsetenv(nm) else do.call(Sys.setenv, stats::setNames(list(old_env[nm]), nm))
    }
  }, add = TRUE)

  Sys.setenv(
    OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1",
    VECLIB_MAXIMUM_THREADS = "1", NUMEXPR_NUM_THREADS = "1"
  )

  if (.Platform$OS.type == "windows" && n_cores > 1L) {
    if (verbose) message(">>> Parallelisation limited to 1 core on Windows.")
    n_cores <- 1L
  }
  n_cores <- max(1L, as.integer(n_cores))

  pilot_data <- object@pilot_data
  groups <- unique(vapply(pilot_data, function(x) as.character(x$group), character(1)))
  all_genes <- rownames(pilot_data[[1]]$counts)

  expr_params <- list()
  if (!is.null(object@params_expression$top_genes)) expr_params$top_genes <- object@params_expression$top_genes
  if (!is.null(object@params_expression$stable_genes)) expr_params$stable_genes <- object@params_expression$stable_genes
  if (!is.null(object@params_expression$stable_gene_stats)) expr_params$stable_gene_stats <- object@params_expression$stable_gene_stats

  if (is.null(genes_to_use)) {
    genes_to_use <- if (!is.null(expr_params$top_genes)) expr_params$top_genes else all_genes
  }
  genes_to_use <- intersect(genes_to_use, all_genes)
  if (length(genes_to_use) == 0) stop("No valid genes to estimate.")

  if (verbose) {
    message(sprintf(">>> Starting Multi-Sample Estimation for %d genes across %d groups (n_cores=%d)...",
                    length(genes_to_use), length(groups), n_cores))
  }

  for (grp in groups) {
    grp_char <- as.character(grp)
    if (verbose) message(sprintf("    Processing Group %s...", grp_char))

    grp_samples <- pilot_data[vapply(pilot_data, function(x) as.character(x$group) == grp, logical(1))]
    K_grp <- length(grp_samples)

    if (K_grp <= 1) stop(sprintf("Group %s has K=%d. This function requires K>1.", grp_char, K_grp))

    raw_domains <- unlist(lapply(grp_samples, function(x) as.character(x$coords$domain)))
    all_domains <- sort(unique(raw_domains))
    if (requireNamespace("stringr", quietly = TRUE)) {
      all_domains <- stringr::str_sort(all_domains, numeric = TRUE)
    }

    if (verbose) message("    -> Pre-processing Coordinates...")
    sample_cache <- vector("list", K_grp)
    for (k in seq_len(K_grp)) {
      coords_mat <- as.matrix(grp_samples[[k]]$coords[, c("x", "y"), drop = FALSE])
      n_obs <- nrow(coords_mat)
      set.seed(k)
      coords_mat <- matrix(as.double(coords_mat), ncol = 2)
      coords_mat[, 1] <- coords_mat[, 1] + stats::runif(n_obs, -1e-6, 1e-6)
      coords_mat[, 2] <- coords_mat[, 2] + stats::runif(n_obs, -1e-6, 1e-6)
      intercept_x <- matrix(1, nrow = n_obs, ncol = 1)
      storage.mode(intercept_x) <- "double"

      if (any(!is.finite(coords_mat))) { sample_cache[[k]] <- NULL; next }
      tryCatch({
        ord <- BRISC::BRISC_order(coords_mat, order = "AMMD", verbose = FALSE)
        sample_cache[[k]] <- list(
          coords = coords_mat, x_cov = intercept_x, order = as.integer(ord),
          doms = as.character(grp_samples[[k]]$coords$domain)
        )
      }, error = function(e) { sample_cache[[k]] <<- NULL })
    }

    logcount_cache <- lapply(seq_len(K_grp), function(k) {
      if (is.null(sample_cache[[k]])) return(NULL)
      raw_counts <- grp_samples[[k]]$counts
      available <- intersect(genes_to_use, rownames(raw_counts))
      as.matrix(log1p(suppressWarnings(as.matrix(raw_counts[available, , drop = FALSE]))))
    })

    if (verbose) message("    -> Pre-computing statistics and residuals via Rcpp...")
    dom_idx_list <- lapply(seq_len(K_grp), function(k) match(sample_cache[[k]]$doms, all_domains) - 1L)

    cpp_res <- compute_expr_stats_cpp(
      logcounts_list = logcount_cache, dom_idx_list = dom_idx_list,
      n_genes = length(genes_to_use), n_domains = length(all_domains), n_cores = n_cores
    )

    pre_mu_grand <- cpp_res$mu_grand
    colnames(pre_mu_grand) <- all_domains
    pre_sigma_bio <- cpp_res$sigma_bio
    pre_residuals <- cpp_res$residuals

    if (verbose) message(sprintf("    -> Fitting Spatial Model for %d genes...", length(genes_to_use)))

    .fit_one_gene <- function(g_idx) {
      theta_matrix <- matrix(NA_real_, nrow = K_grp, ncol = 3)
      for (k in seq_len(K_grp)) {
        if (is.null(sample_cache[[k]])) next
        res_vec <- pre_residuals[[k]][g_idx, ]
        if (stats::var(res_vec) < 1e-8) next
        tryCatch({
          fit <- BRISC::BRISC_estimation(
            coords = sample_cache[[k]]$coords, y = res_vec, x = sample_cache[[k]]$x_cov,
            neighbor = NULL, ordering = sample_cache[[k]]$order, n.neighbors = as.integer(n_neighbors),
            n_omp = 1, verbose = FALSE, cov.model = "exponential"
          )
          th <- fit$Theta
          theta_matrix[k, ] <- c(
            alpha    = if ("sigma.sq" %in% names(th)) th["sigma.sq"] else th[1],
            rho      = if ("phi" %in% names(th))      th["phi"]      else th[2],
            sigma.sq = if ("tau.sq" %in% names(th))   th["tau.sq"]   else th[3]
          )
        }, error = function(e) NULL)
      }
      pooled_theta <- colMeans(theta_matrix, na.rm = TRUE)
      if (any(is.na(pooled_theta))) return(NULL)
      list(theta = pooled_theta, sigma_bio = pre_sigma_bio[g_idx], mu_grand = pre_mu_grand[g_idx, ])
    }

    if (n_cores > 1L) {
      gene_results_list <- parallel::mclapply(seq_along(genes_to_use), .fit_one_gene, mc.cores = n_cores, mc.set.seed = TRUE)
    } else {
      gene_results_list <- lapply(seq_along(genes_to_use), function(i) {
        if (verbose && i %% 5 == 0) cat(sprintf("\r      [Fitting] %d/%d", i, length(genes_to_use)))
        .fit_one_gene(i)
      })
    }
    names(gene_results_list) <- genes_to_use
    gene_results_list <- gene_results_list[!vapply(gene_results_list, is.null, logical(1))]

    if (verbose) message(sprintf("\n    -> Group %s: %d / %d genes fitted successfully.", grp_char, length(gene_results_list), length(genes_to_use)))
    expr_params[[grp_char]] <- gene_results_list
  }

  expr_params$K_regime <- "multi"
  expr_params$calibration <- list(target_domain = target_domain, reference_domain = reference_domain)
  object@params_expression <- expr_params

  if (verbose) message(">>> Estimation Completed Successfully.")
  return(object)
}

# Internal Helper for NULL coalescing
#' @keywords internal
#' @noRd
`%||%` <- function(a, b) if (!is.null(a)) a else b

#' @title Estimate Spatial Geometry Parameters (Multi-Sample)
#' @description
#' Estimates parameters for domain geometry via FGKMM pooling across multiple samples.
#'
#' @param object A \code{spaCraft} object.
#' @param M_candidates Integer vector. Candidate M values for FGKMM.
#' @param iter_max Integer. Max EM iterations.
#' @param tol Numeric. EM convergence tolerance.
#' @param n_cores Integer. Number of PSOCK cores.
#' @param verbose Logical. Print progress.
#' @param min_spots Integer. Minimum spots per domain per sample.
#' @param trim_outliers Logical. Density-based outlier trimming before FGKMM fitting.
#' @param trim_k Integer. KNN for density estimation.
#' @param trim_quantile Numeric in (0,1).
#' @param gc_each_fit Logical. Trigger garbage collection each fit.
#' @param store_criteria Logical. Store BIC criteria table.
#' @param worker_sources Character vector. Optional R scripts to source on PSOCK workers.
#' @param max_pool_spots Integer. Max pooled centered spots for M selection.
#'
#' @return Updated \code{spaCraft} object with \code{params_geometry}.
#' @importFrom stats cov quantile complete.cases
#' @export
#'
estimateGeometryParams <- function(object,
                                      M_candidates = 3:6,
                                      iter_max = 500,
                                      tol = 1e-2,
                                      n_cores = 1,
                                      verbose = TRUE,
                                      min_spots = 30,
                                      trim_outliers = TRUE,
                                      trim_k = 10L,
                                      trim_quantile = 0.98,
                                      gc_each_fit = FALSE,
                                      store_criteria = TRUE,
                                      worker_sources = NULL,
                                      max_pool_spots = 5000L) {

  # CRAN Safe: Save and restore environment variables
  old_env <- Sys.getenv(c("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"))
  on.exit({
    for (nm in names(old_env)) {
      if (old_env[nm] == "") Sys.unsetenv(nm) else do.call(Sys.setenv, stats::setNames(list(old_env[nm]), nm))
    }
  }, add = TRUE)

  Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1")

  if (verbose) message(">>> Pre-normalizing all pilot coordinates to [0,1]...")
  normalized_pilot <- lapply(object@pilot_data, function(ps) {
    coords <- ps$coords
    mx <- max(coords$x, na.rm = TRUE); mnx <- min(coords$x, na.rm = TRUE)
    my <- max(coords$y, na.rm = TRUE); mny <- min(coords$y, na.rm = TRUE)
    dx <- if ((mx - mnx) < 1e-6) 1 else (mx - mnx)
    dy <- if ((my - mny) < 1e-6) 1 else (my - mny)
    coords$x <- (coords$x - mnx) / dx
    coords$y <- (coords$y - mny) / dy
    ps$coords <- coords
    ps
  })

  groups <- sort(unique(vapply(normalized_pilot, function(x) as.character(x$group), character(1))))
  FIT_SCALE <- 1.0

  .sort_domains <- function(dom_vec) {
    doms <- unique(as.character(dom_vec))
    if (requireNamespace("stringr", quietly = TRUE)) stringr::str_sort(doms, numeric = TRUE) else sort(doms)
  }

  .rescale_model <- function(model, factor) {
    if (is.null(model)) return(NULL)
    model$c <- model$c * factor
    model$r <- model$r * factor
    model$sigma_sq <- model$sigma_sq * (factor^2)
    model
  }

  .trim_outlier_spots <- function(coords_mat, k_nn = 10L, q_cut = 0.95) {
    n <- nrow(coords_mat)
    if (n <= k_nn + 1L) return(seq_len(n))
    if (!requireNamespace("FNN", quietly = TRUE)) {
      warning("FNN not available; skipping outlier trimming.", call. = FALSE)
      return(seq_len(n))
    }
    k_use <- min(as.integer(k_nn), n - 1L)
    nn <- FNN::get.knn(coords_mat, k = k_use)
    avg_dist <- rowMeans(nn$nn.dist)
    threshold <- stats::quantile(avg_dist, q_cut)
    keep <- which(avg_dist <= threshold)
    if (length(keep) < n * 0.8) keep <- which(avg_dist <= stats::quantile(avg_dist, 0.80))
    if (length(keep) < min_spots) keep <- order(avg_dist)[seq_len(min(n, max(min_spots, 30L)))]
    keep
  }

  .anchor_global_theta <- function(theta_k, omega_k, p = 0.15, a0_val = 50) {
    M <- length(theta_k)
    global_theta <- lapply(seq_len(M), function(m) {
      c_m <- as.numeric(theta_k[[m]]$c); r_m <- as.numeric(theta_k[[m]]$r)
      phi_m <- as.numeric(theta_k[[m]]$phi); tau_m <- as.numeric(theta_k[[m]]$tau)
      sigma2_m <- as.numeric(theta_k[[m]]$sigma2)
      mu_c <- c_m; Sigma_c <- diag((p * pmax(abs(c_m), 1e-5))^2, 2)
      mu_r <- r_m; var_r <- (p * pmax(abs(r_m), 1e-5))^2
      R_bar <- sqrt(sum(phi_m^2))
      mu_phi <- if (R_bar > 0) phi_m / R_bar else c(1, 0)
      kappa_phi_val <- estimate_kappa(min(max(R_bar, 1e-3), 0.999), d = 2)
      mean_tau <- tau_m; var_tau <- (p * pmax(abs(tau_m), 1e-3))^2
      a_tau <- mean_tau^2 / pmax(var_tau, 1e-8); b_tau <- mean_tau / pmax(var_tau, 1e-8)
      inv_sigma <- 1 / pmax(sigma2_m, 1e-12); mean_inv <- inv_sigma
      var_inv <- (p * pmax(abs(inv_sigma), 1e-3))^2
      a_sigma <- mean_inv^2 / pmax(var_inv, 1e-8); b_sigma <- mean_inv / pmax(var_inv, 1e-8)
      list(mu_c = mu_c, Sigma_c = Sigma_c, mu_r = mu_r, var_r = var_r,
           mu_phi = mu_phi, kappa_phi = kappa_phi_val, a_tau = a_tau, b_tau = b_tau,
           a_sigma = a_sigma, b_sigma = b_sigma,
           center_cloud = list(c_m), orientation_cloud = list(phi_m))
    })
    alpha_hat <- pmax(a0_val * as.numeric(omega_k), 1e-6)
    list(global_theta = global_theta, alpha_hat = alpha_hat)
  }

  .fit_select_M_once <- function(x_fitting, downscale_factor) {
    if (gc_each_fit) gc(verbose = FALSE)
    out <- tryCatch({
      sel <- select_best_M(x = x_fitting, M_candidates = M_candidates, iter_max = iter_max, tol = tol)
      model_rescaled <- .rescale_model(sel$best_model_BIC, downscale_factor)
      list(ok = TRUE, model = model_rescaled, best_M = sel$best_M_BIC, criteria = sel$criteria_table)
    }, error = function(e) list(ok = FALSE, error = e$message))
    out
  }

  .fit_forced_M <- function(x_fitting, forced_M, downscale_factor) {
    if (gc_each_fit) gc(verbose = FALSE)
    out <- tryCatch({
      fit <- FG_EM(x = x_fitting, M = forced_M, iter_max = iter_max, tol = tol)
      model_rescaled <- .rescale_model(fit, downscale_factor)
      list(ok = TRUE, model = model_rescaled, best_M = forced_M)
    }, error = function(e) list(ok = FALSE, error = e$message))
    out
  }

  cl <- NULL
  if (n_cores > 1) {
    if (verbose) message(sprintf(">>> Initializing worker cluster (n_cores = %d)...", n_cores))
    cl <- parallel::makeCluster(n_cores, type = "PSOCK")
    on.exit(if (!is.null(cl)) parallel::stopCluster(cl), add = TRUE)

    parallel::clusterEvalQ(cl, {
      library(stats); library(matrixStats); library(Directional)
    })

    if (!is.null(worker_sources)) {
      parallel::clusterExport(cl, varlist = "worker_sources", envir = environment())
      parallel::clusterEvalQ(cl, { for (f in worker_sources) source(f); NULL })
    }

    parallel::clusterExport(cl, varlist = c(
      "%||%", "select_best_M", "FG_EM", "FG_EM_with_criteria",
      "initial_kmeans", "estimate_kappa", "M_candidates", "iter_max", "tol", "gc_each_fit",
      ".fit_select_M_once", ".fit_forced_M", ".rescale_model"
    ), envir = environment())
  }

  if (.Platform$OS.type == "windows" && n_cores > 1L) {
    if (verbose) message(">>> Parallelisation limited to 1 core on Windows.")
    n_cores <- 1L
  }

  .lapply_cluster <- function(X, FUN) {
    if (n_cores > 1L) parallel::mclapply(X, FUN, mc.cores = n_cores, mc.set.seed = TRUE) else lapply(X, FUN)
  }

  geometry_params <- list()
  if (verbose && trim_outliers) message(sprintf(">>> Outlier trimming ENABLED (k=%d, quantile=%.2f)", trim_k, trim_quantile))

  for (grp in groups) {
    grp_char <- as.character(grp)
    if (verbose) message(sprintf(">>> Processing Geometry for Group %s...", grp_char))

    grp_samples <- normalized_pilot[vapply(normalized_pilot, function(x) as.character(x$group) == grp, logical(1))]
    K_grp <- length(grp_samples)
    if (K_grp <= 1) stop(sprintf("Group %s has K=%d. This function requires K>1.", grp_char, K_grp))

    raw_domains <- unlist(lapply(grp_samples, function(x) as.character(x$coords$domain)))
    all_domains <- .sort_domains(raw_domains)
    group_results <- list()

    for (dom in all_domains) {
      if (verbose) message(sprintf("    -> Domain %s", dom))

      centroids <- matrix(NA_real_, nrow = K_grp, ncol = 2L); colnames(centroids) <- c("x", "y")
      per_sample_items <- vector("list", K_grp)
      pooled_centered_list <- list()
      n_trimmed_total <- 0L; n_original_total <- 0L

      for (k in seq_len(K_grp)) {
        dom_coords <- as.matrix(grp_samples[[k]]$coords[as.character(grp_samples[[k]]$coords$domain) == dom, c("x", "y"), drop = FALSE])
        if (nrow(dom_coords) < min_spots) { per_sample_items[[k]] <- NULL; next }

        centroid <- colMeans(dom_coords); centroids[k, ] <- centroid
        x_centered <- sweep(dom_coords, 2, centroid, "-", check.margin = FALSE) * FIT_SCALE
        n_original <- nrow(x_centered); n_original_total <- n_original_total + n_original

        if (trim_outliers && n_original > trim_k + 1L) {
          keep_idx <- .trim_outlier_spots(x_centered, k_nn = trim_k, q_cut = trim_quantile)
          x_fitting <- x_centered[keep_idx, , drop = FALSE]
          n_trimmed_total <- n_trimmed_total + (n_original - length(keep_idx))
        } else {
          x_fitting <- x_centered
        }
        per_sample_items[[k]] <- list(k = k, x_fitting = x_fitting, centroid = centroid)
        pooled_centered_list[[length(pooled_centered_list) + 1L]] <- x_fitting
      }

      valid_items <- per_sample_items[!vapply(per_sample_items, is.null, logical(1))]
      if (length(valid_items) == 0L) {
        if (verbose) message(sprintf("      ! skip: no sample has >= %d spots", min_spots))
        next
      }

      if (verbose && trim_outliers && n_trimmed_total > 0L) {
        message(sprintf("      -> Trimmed %d / %d outlier spots (%.1f%%)",
                        n_trimmed_total, n_original_total, 100 * n_trimmed_total / max(n_original_total, 1)))
      }

      centroids_valid <- centroids[stats::complete.cases(centroids), , drop = FALSE]
      mu_G <- colMeans(centroids_valid)
      Sigma_G <- stats::cov(centroids_valid) + diag(1e-8, 2L)
      placement <- list(mu_G = mu_G, Sigma_G = Sigma_G)

      pooled_x <- do.call(rbind, pooled_centered_list)
      if (nrow(pooled_x) > max_pool_spots) {
        set.seed(123)
        pooled_x <- pooled_x[sample(seq_len(nrow(pooled_x)), max_pool_spots), , drop = FALSE]
      }

      if (verbose) message(sprintf("      -> Selecting M on pooled coords (n=%d)...", nrow(pooled_x)))
      dom_select <- .fit_select_M_once(pooled_x, 1.0 / FIT_SCALE)
      if (!isTRUE(dom_select$ok)) {
        if (verbose) message("      ! model selection failed: ", dom_select$error); next
      }
      consensus_M <- dom_select$best_M
      criteria_tbl <- dom_select$criteria %||% NULL
      if (verbose) message(sprintf("      -> Using consensus M = %d", consensus_M))

      fit_fun_wrapper <- function(item) .fit_forced_M(item$x_fitting, consensus_M, 1.0 / FIT_SCALE)
      fit_results <- .lapply_cluster(valid_items, fit_fun_wrapper)
      ok_idx <- vapply(fit_results, function(z) isTRUE(z$ok), logical(1))
      if (!any(ok_idx)) { if (verbose) message("      ! all forced-M fits failed"); next }

      fit_ok <- fit_results[ok_idx]
      items_ok <- valid_items[ok_idx]

      theta_list <- lapply(fit_ok, function(f) {
        lapply(seq_len(consensus_M), function(m) {
          list(c = f$model$c[m, ], r = f$model$r[m], phi = f$model$phi[m, ],
               tau = f$model$tau[m], sigma2 = f$model$sigma_sq)
        })
      })
      omega_list <- lapply(fit_ok, function(f) f$model$pi)

      pooled <- tryCatch(pool_FGKMM_parameters(theta_list, omega_list), error = function(e) NULL)
      if (is.null(pooled)) {
        if (verbose) message("      -> pooling failed; using old-style anchoring")
        pooled <- .anchor_global_theta(theta_list[[1]], omega_list[[1]])
      }

      pooled$per_sample_sigma2 <- vapply(fit_ok, function(f) as.numeric(f$model$sigma_sq), numeric(1))

      fit_meta <- lapply(seq_along(fit_ok), function(j) {
        list(
          k = items_ok[[j]]$k, centroid = items_ok[[j]]$centroid,
          best_M = consensus_M, sigma2 = as.numeric(fit_ok[[j]]$model$sigma_sq),
          criteria = if (store_criteria) criteria_tbl else NULL
        )
      })

      group_results[[as.character(dom)]] <- list(
        placement = placement, pattern = list(M = consensus_M, pooled = pooled, scale = c(1.0, 1.0)),
        fit_meta = fit_meta
      )
    }
    geometry_params[[grp_char]] <- group_results
  }

  geometry_params$K_regime <- "multi"
  object@params_geometry <- geometry_params
  if (verbose) message(">>> Geometry estimation completed successfully.")
  return(object)
}

#' @title Estimate Composition Parameters (Multi-Sample)
#' @description
#' Estimates parameters for the composition model via Binomial GLM and empirical
#' pilot moments. Requires K>1 replicates per group.
#'
#' @param object A \code{spaCraft} object containing \code{pilot_data}.
#' @param target_domain Character. Target domain name.
#' @param reference_domain Character. Reference domain name.
#' @param prob_floor Numeric. Probability clamp for numerical safety.
#' @param verbose Logical. Print progress.
#'
#' @return Updated \code{spaCraft} object with \code{params_composition}.
#' @importFrom stats glm binomial coef qlogis plogis var
#' @export
#'
estimateCompositionParams <- function(object,
                                         target_domain,
                                         reference_domain,
                                         prob_floor = 1e-6,
                                         verbose = TRUE) {

  msg <- function(...) if (isTRUE(verbose)) message(...)
  clamp01 <- function(p, eps = prob_floor) pmin(pmax(p, eps), 1 - eps)

  pilot_data <- object@pilot_data
  if (is.null(pilot_data) || !is.list(pilot_data) || length(pilot_data) == 0) {
    stop("object@pilot_data must be a non-empty list.")
  }

  has_coords <- vapply(pilot_data, function(s) !is.null(s$coords), logical(1))
  if (!all(has_coords)) stop("Each pilot_data sample must contain $coords.")

  has_domain <- vapply(pilot_data, function(s) "domain" %in% colnames(s$coords), logical(1))
  if (!all(has_domain)) stop("Each pilot_data coords must have a 'domain' column.")

  has_group <- vapply(pilot_data, function(s) !is.null(s$group), logical(1))
  if (!all(has_group)) stop("Each pilot_data sample must contain $group (0/1).")

  all_domains <- sort(unique(unlist(lapply(pilot_data, function(s) as.character(s$coords$domain)))))
  if (!target_domain %in% all_domains) stop(sprintf("Target domain '%s' not found.", target_domain))
  if (!reference_domain %in% all_domains) stop(sprintf("Reference domain '%s' not found.", reference_domain))

  background_domains <- setdiff(all_domains, c(target_domain, reference_domain))

  grp_raw <- vapply(pilot_data, function(s) as.character(s$group), character(1))
  grp_num <- suppressWarnings(as.numeric(grp_raw))
  if (any(!is.finite(grp_num))) stop("Group labels must be coercible to numeric (typically 0/1).")
  if (!all(grp_num %in% c(0, 1))) stop("Group labels must be 0/1 for the current implementation.")

  K_per_group <- table(grp_num)
  if (any(K_per_group <= 1)) stop("This function strictly requires K > 1 per group.")

  msg(">>> Composition Parameter Estimation (Multi-Sample)")
  msg(sprintf("    Target: %s / Reference: %s", target_domain, reference_domain))
  msg(sprintf("    Background domains: %d | Regime: multi", length(background_domains)))

  sample_ids <- names(pilot_data)
  if (is.null(sample_ids) || any(sample_ids == "")) sample_ids <- paste0("s", seq_along(pilot_data))

  comp_list <- vector("list", length(pilot_data))
  names(comp_list) <- sample_ids

  for (i in seq_along(pilot_data)) {
    s <- pilot_data[[i]]
    dom <- as.character(s$coords$domain)
    tab <- table(dom)

    n_target <- sum(dom == target_domain)
    n_ref    <- sum(dom == reference_domain)

    bg_counts <- if (length(background_domains) > 0) {
      tmp <- as.numeric(tab[background_domains])
      tmp[is.na(tmp)] <- 0
      names(tmp) <- background_domains
      tmp
    } else numeric(0)

    comp_list[[i]] <- list(
      sample_id = sample_ids[i], y = n_target, m = n_target + n_ref,
      grp = grp_num[i], bg_counts = bg_counts, total_N = length(dom)
    )
  }

  df_binom <- do.call(rbind.data.frame, lapply(comp_list, function(x) {
    data.frame(sample_id = x$sample_id, y = as.integer(x$y), m = as.integer(x$m), grp = as.numeric(x$grp), stringsAsFactors = FALSE)
  }))
  rownames(df_binom) <- NULL
  df_valid <- df_binom[df_binom$m > 0, , drop = FALSE]

  est_beta0 <- NA_real_; est_beta1 <- NA_real_; fit_method <- NULL
  has_both_groups <- length(unique(df_valid$grp)) == 2

  if (nrow(df_valid) < 2 || !has_both_groups) {
    msg("   -> Insufficient group contrast; falling back to pooled proportions.")
    p_pool <- if (sum(df_valid$m) > 0) sum(df_valid$y) / sum(df_valid$m) else 0.5
    est_beta0 <- stats::qlogis(clamp01(p_pool))
    est_beta1 <- 0
    fit_method <- "pooled_mean_fallback"
  } else {
    fit <- tryCatch(stats::glm(cbind(y, m - y) ~ grp, data = df_valid, family = stats::binomial()), error = function(e) NULL)
    if (is.null(fit)) {
      msg("   -> GLM failed; falling back to pooled proportions.")
      p_pool <- if (sum(df_valid$m) > 0) sum(df_valid$y) / sum(df_valid$m) else 0.5
      est_beta0 <- stats::qlogis(clamp01(p_pool))
      est_beta1 <- 0
      fit_method <- "glm_failed_pooled_fallback"
    } else {
      cc <- stats::coef(fit)
      est_beta0 <- unname(cc[1])
      est_beta1 <- unname(if ("grp" %in% names(cc)) cc["grp"] else 0)
      fit_method <- "binomial_glm_mle"
    }
  }

  msg(sprintf("   -> Fit: %s", fit_method))
  msg(sprintf("      beta0 (Baseline Log-odds): %.4f", est_beta0))
  msg(sprintf("      beta1 (Effect Size): %.4f (OR=%.3f)", est_beta1, exp(est_beta1)))

  mu0_hat <- stats::plogis(est_beta0)
  mu1_hat <- stats::plogis(est_beta0 + est_beta1)

  bg_props <- list()
  if (length(background_domains) > 0) {
    grps <- sort(unique(df_binom$grp))
    for (g in grps) {
      idx <- which(vapply(comp_list, function(x) x$grp == g, logical(1)))
      samps <- comp_list[idx]
      props_mat <- do.call(rbind, lapply(samps, function(x) {
        R <- x$total_N - x$m
        if (R <= 0) return(rep(NA_real_, length(background_domains)))
        as.numeric(x$bg_counts) / R
      }))
      colnames(props_mat) <- background_domains
      avg <- colMeans(props_mat, na.rm = TRUE)
      avg[!is.finite(avg)] <- 0
      if (sum(avg) <= 0) {
        avg <- rep(1 / length(background_domains), length(background_domains))
        names(avg) <- background_domains
      } else {
        avg <- avg / sum(avg)
      }
      bg_props[[as.character(g)]] <- avg
    }
  }

  m_vals <- vapply(comp_list, function(x) x$m, numeric(1))
  contrast_stats <- list(
    mu = mean(m_vals), sigma2 = stats::var(m_vals), fixed = FALSE,
    by_group = lapply(split(m_vals, df_binom$grp), function(v) list(mu = mean(v), sigma2 = stats::var(v), fixed = FALSE))
  )

  object@params_composition <- list(
    beta_binomial = list(
      beta0 = est_beta0, beta1 = est_beta1, phi = NA_real_,
      mu0_hat = mu0_hat, mu1_hat = mu1_hat, target = target_domain,
      reference = reference_domain, fit_method = fit_method, model = "binomial_logistic"
    ),
    background = list(domains = background_domains, proportions = bg_props),
    contrast_size = contrast_stats, K_regime = "multi"
  )

  msg(">>> Composition Parameter Estimation Completed.")
  return(object)
}
