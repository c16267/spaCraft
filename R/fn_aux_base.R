# ==============================================================================
# File: fn_aux_base.R
# Internal helper functions for spatial generation, extraction, and evaluation
# ==============================================================================
utils::globalVariables(c("geneSummary"))


extract_theta_triplet <- function(theta) {
  if (is.null(theta)) {
    return(c(alpha = 0.2, rho = 10, nugget = 0.05))
  }
  nm <- names(theta)
  alpha_sp <- NA_real_; rho_sp <- NA_real_; nugget2 <- NA_real_

  if ("alpha" %in% nm) alpha_sp <- as.numeric(theta[["alpha"]])
  if ("rho"   %in% nm) rho_sp   <- as.numeric(theta[["rho"]])

  if ("nugget" %in% nm) {
    nugget2 <- as.numeric(theta[["nugget"]])
  } else if ("alpha" %in% nm && "sigma.sq" %in% nm) {
    nugget2 <- as.numeric(theta[["sigma.sq"]])
  }

  if (!is.finite(alpha_sp) && "sigma.sq" %in% nm) alpha_sp <- as.numeric(theta[["sigma.sq"]])
  if (!is.finite(rho_sp)   && "phi"      %in% nm) rho_sp   <- as.numeric(theta[["phi"]])
  if (!is.finite(nugget2)  && "tau.sq"   %in% nm) nugget2  <- as.numeric(theta[["tau.sq"]])

  if (!is.finite(alpha_sp) || alpha_sp < 0) alpha_sp <- 0
  if (!is.finite(rho_sp)   || rho_sp <= 0)  rho_sp   <- 10
  if (!is.finite(nugget2)  || nugget2 < 0)  nugget2  <- 0

  c(alpha = alpha_sp, rho = rho_sp, nugget = nugget2)
}

low_rank_gp_draw_from_cov <- function(K, rank_max = 20L, jitter = 1e-8) {
  n <- nrow(K)
  if (n == 0L) return(numeric(0))
  if (n == 1L) {
    v <- max(K[1, 1], 0)
    return(stats::rnorm(1, mean = 0, sd = sqrt(v)))
  }
  diag(K) <- diag(K) + jitter
  out <- tryCatch({
    R <- chol(K, pivot = TRUE)
    piv <- attr(R, "pivot")
    rk  <- attr(R, "rank")
    if (is.null(rk) || !is.finite(rk)) rk <- n
    r_use <- min(as.integer(rank_max), as.integer(rk), n)
    r_use <- max(r_use, 1L)
    R_r <- R[seq_len(r_use), , drop = FALSE]
    z_r <- stats::rnorm(r_use)
    x_piv <- as.numeric(crossprod(R_r, z_r))
    x <- numeric(n)
    x[piv] <- x_piv
    x
  }, error = function(e) {
    tryCatch({
      U <- chol(K, pivot = FALSE)
      as.numeric(t(U) %*% stats::rnorm(n))
    }, error = function(e2) {
      stats::rnorm(n, mean = 0, sd = sqrt(max(mean(diag(K)), 0)))
    })
  })
  out[!is.finite(out)] <- 0
  out
}

rescale_spatial_field <- function(x, alpha) {
  x <- as.numeric(x)
  x[!is.finite(x)] <- 0
  x <- x - mean(x)
  sx <- stats::sd(x)
  if (!is.finite(sx) || sx < 1e-10 || !is.finite(alpha) || alpha <= 0) {
    return(rep(0, length(x)))
  }
  sqrt(alpha) * (x / sx)
}

draw_graph_gp_whole_tissue <- function(spatial_cache, alpha, rho, include_self = TRUE) {
  n <- spatial_cache$n_cells
  if (n <= 0L || alpha <= 0) return(numeric(n))
  if (!is.finite(rho) || rho <= 0) rho <- 10
  nn_index <- spatial_cache$nn_index
  nn_dist  <- spatial_cache$nn_dist
  k_use <- ncol(nn_index)
  z <- stats::rnorm(n)

  if (k_use == 0L) return(rescale_spatial_field(z, alpha))
  W <- exp(-nn_dist / rho)

  if (include_self) {
    z_nb <- cbind(z, matrix(z[nn_index], nrow = n))
    W_all <- cbind(rep(1, n), W)
  } else {
    z_nb <- matrix(z[nn_index], nrow = n)
    W_all <- W
  }
  rs <- rowSums(W_all)
  rs[!is.finite(rs) | rs <= 0] <- 1
  W_all <- W_all / rs
  eta <- rowSums(W_all * z_nb)
  rescale_spatial_field(eta, alpha)
}

build_basis_matrix <- function(coords_mat, knots, rho,
                               basis_kernel = c("gaussian", "exponential"),
                               normalize_columns = TRUE) {
  basis_kernel <- match.arg(basis_kernel)
  n <- nrow(coords_mat); m <- nrow(knots)
  if (n == 0L || m == 0L) return(matrix(0, n, m))
  if (!is.finite(rho) || rho <= 0) rho <- 10
  dx <- outer(coords_mat[, 1], knots[, 1], "-")
  dy <- outer(coords_mat[, 2], knots[, 2], "-")
  d  <- sqrt(dx^2 + dy^2)
  Phi <- switch(basis_kernel,
                gaussian    = exp(-0.5 * (d / rho)^2),
                exponential = exp(-d / rho))
  if (normalize_columns && ncol(Phi) > 0L) {
    cn <- sqrt(colSums(Phi^2))
    cn[!is.finite(cn) | cn <= 0] <- 1
    Phi <- sweep(Phi, 2, cn, "/", check.margin = FALSE)
  }
  Phi
}

draw_basis_gp_whole_tissue <- function(spatial_cache, alpha, rho,
                                       basis_kernel = c("gaussian", "exponential")) {
  basis_kernel <- match.arg(basis_kernel)
  n <- spatial_cache$n_cells
  if (n <= 0L || !is.finite(alpha) || alpha <= 0) return(numeric(n))
  coords_mat <- spatial_cache$coords_mat
  knots <- spatial_cache$basis_knots
  if (is.null(knots) || nrow(knots) == 0L) return(rescale_spatial_field(stats::rnorm(n), alpha))
  Phi <- build_basis_matrix(coords_mat = coords_mat, knots = knots, rho = rho,
                            basis_kernel = basis_kernel, normalize_columns = TRUE)
  m <- ncol(Phi)
  if (m == 0L) return(rescale_spatial_field(stats::rnorm(n), alpha))
  z <- stats::rnorm(m)
  eta <- as.numeric(Phi %*% z)
  eta[!is.finite(eta)] <- 0
  rescale_spatial_field(eta, alpha)
}

gene_worker <- function(g, sim_coords_df, sim_coords_mat, sim_domains,
                        pilot_bundle = NULL, p_info_override = NULL, params_expr_grp = NULL,
                        is_case, target_domain, de_genes, scenario, k_nn = 20L,
                        spatial_cache = NULL, spatial_mode = c("original", "basis"),
                        gp_rank = 40L, gp_jitter = 1e-8,
                        basis_kernel = c("gaussian", "exponential"),
                        pilot_sample = NULL, pilot_coords_mat = NULL) {
  if (!requireNamespace("FNN", quietly = TRUE)) stop("Package 'FNN' required.")
  spatial_mode <- match.arg(spatial_mode)
  basis_kernel <- match.arg(basis_kernel)
  sim_coords_mat <- as.matrix(sim_coords_mat)
  n_cells <- nrow(sim_coords_mat)
  if (n_cells == 0L) return(numeric(0))

  domain_vec <- if ("tilde_d" %in% names(sim_coords_df)) {
    as.character(sim_coords_df$tilde_d)
  } else {
    as.character(sim_coords_df$domain)
  }

  if (!is.null(p_info_override)) {
    mu_grand  <- p_info_override$mu_grand
    theta     <- p_info_override$theta
    sigma_bio <- p_info_override$sigma_bio
  } else {
    p_ex <- if (!is.null(params_expr_grp)) params_expr_grp[[g]] else NULL
    if (is.null(p_ex)) {
      mu_grand  <- stats::setNames(rep(0.1, length(sim_domains)), sim_domains)
      theta     <- c(alpha = 0.2, rho = 10, sigma.sq = 0.05)
      sigma_bio <- 0
    } else {
      mu_grand  <- p_ex$mu_grand
      theta     <- p_ex$theta
      sigma_bio <- if (!is.null(p_ex$sigma_bio)) p_ex$sigma_bio else 0
    }
  }

  mu_grand <- as.numeric(mu_grand)
  if (!is.null(p_info_override) && !is.null(names(p_info_override$mu_grand))) {
    names(mu_grand) <- names(p_info_override$mu_grand)
  } else if (!is.null(params_expr_grp) && !is.null(params_expr_grp[[g]]) && !is.null(names(params_expr_grp[[g]]$mu_grand))) {
    names(mu_grand) <- names(params_expr_grp[[g]]$mu_grand)
  } else {
    names(mu_grand) <- sim_domains
  }

  th <- extract_theta_triplet(theta)
  alpha_sp <- as.numeric(th["alpha"]); rho_sp <- as.numeric(th["rho"]); nugget2 <- as.numeric(th["nugget"])
  if (!is.finite(alpha_sp) || alpha_sp < 0) alpha_sp <- 0
  if (!is.finite(rho_sp)   || rho_sp <= 0)  rho_sp   <- 10
  if (!is.finite(nugget2)  || nugget2 < 0)  nugget2  <- 0

  rho_eff <- rho_sp
  if (is_case) {
    delta_rho <- as.numeric(scenario$delta_rho)
    if (is.finite(delta_rho) && delta_rho > 0) rho_eff <- rho_sp * delta_rho
  }

  current_pilot_data <- NULL; current_pilot_coords <- NULL
  if (!is.null(pilot_bundle) && !is.null(pilot_bundle$pilot_samples) && length(pilot_bundle$pilot_samples) > 0L) {
    w <- pilot_bundle$pilot_mix_w
    if (is.null(w) || length(w) != length(pilot_bundle$pilot_samples) || any(!is.finite(w)) || sum(w) <= 0) {
      w <- rep(1, length(pilot_bundle$pilot_samples))
    }
    w <- w / sum(w)
    sel <- sample(seq_along(pilot_bundle$pilot_samples), size = 1L, prob = w)
    current_pilot_data <- pilot_bundle$pilot_samples[[sel]]
    if (!is.null(pilot_bundle$pilot_coords_mat_list) && length(pilot_bundle$pilot_coords_mat_list) >= sel) {
      current_pilot_coords <- as.matrix(pilot_bundle$pilot_coords_mat_list[[sel]])
    }
  } else {
    current_pilot_data <- pilot_sample
    current_pilot_coords <- if (!is.null(pilot_coords_mat)) as.matrix(pilot_coords_mat) else NULL
  }

  mu_log <- numeric(n_cells)
  mu_fallback <- mean(mu_grand[is.finite(mu_grand)], na.rm = TRUE)
  if (!is.finite(mu_fallback)) mu_fallback <- 0.1
  de_shift <- as.numeric(scenario$DE_lfc)
  if (!is.finite(de_shift)) de_shift <- 0

  for (dom in unique(domain_vec)) {
    base <- if (!is.null(names(mu_grand)) && dom %in% names(mu_grand)) mu_grand[[dom]] else mu_fallback
    if (!is.finite(base)) base <- 0.1
    if (is_case && identical(as.character(dom), as.character(target_domain)) && (g %in% de_genes)) {
      base <- base + de_shift
    }
    mu_log[domain_vec == dom] <- base
  }

  kappa_bio <- as.numeric(scenario$kappa_bio %||% 1)
  if (!is.finite(kappa_bio) || kappa_bio < 0) kappa_bio <- 1
  sigma_bio_new <- as.numeric(sigma_bio) * kappa_bio
  if (!is.finite(sigma_bio_new) || sigma_bio_new < 0) sigma_bio_new <- 0
  B_shift <- stats::rnorm(1, mean = 0, sd = sqrt(sigma_bio_new))

  lambda <- as.numeric(scenario$lambda_cond %||% 0)
  if (!is.finite(lambda)) lambda <- 0
  lambda <- max(0, min(1, lambda))

  eta_cond <- rep(0, n_cells)
  if (lambda > 1e-8 && !is.null(current_pilot_data) && !is.null(current_pilot_coords) && nrow(current_pilot_coords) > 0L) {
    gene_in_logcounts <- !is.null(current_pilot_data$logcounts) && (g %in% rownames(current_pilot_data$logcounts))
    gene_in_counts    <- !is.null(current_pilot_data$counts) && (g %in% rownames(current_pilot_data$counts))

    if (gene_in_logcounts || gene_in_counts) {
      p_log <- if (gene_in_logcounts) as.numeric(current_pilot_data$logcounts[g, ]) else log1p(as.numeric(current_pilot_data$counts[g, ]))
      pilot_dom <- if (!is.null(current_pilot_data$coords) && "domain" %in% names(current_pilot_data$coords)) {
        as.character(current_pilot_data$coords$domain)
      } else rep(NA_character_, length(p_log))

      p_mu <- vapply(pilot_dom, function(d) {
        vv <- if (!is.null(names(mu_grand)) && !is.na(d) && d %in% names(mu_grand)) mu_grand[[d]] else mu_fallback
        if (!is.finite(vv)) vv <- 0.1
        vv
      }, numeric(1))

      p_res <- p_log - p_mu
      p_res[!is.finite(p_res)] <- 0
      p_res <- p_res - mean(p_res)

      for (dom in unique(domain_vec)) {
        syn_idx <- which(domain_vec == dom)
        if (length(syn_idx) == 0L) next
        pilot_idx <- which(pilot_dom == dom)
        if (length(pilot_idx) < 3L) pilot_idx <- seq_along(pilot_dom)
        k_use <- min(as.integer(k_nn), length(pilot_idx))
        k_use <- max(k_use, 1L)

        nn <- FNN::get.knnx(data = current_pilot_coords[pilot_idx, , drop = FALSE],
                            query = sim_coords_mat[syn_idx, , drop = FALSE], k = k_use)
        d2 <- nn$nn.dist^2
        h  <- stats::median(d2[is.finite(d2)], na.rm = TRUE)
        if (!is.finite(h) || h <= 0) h <- mean(d2[is.finite(d2)], na.rm = TRUE)
        if (!is.finite(h) || h <= 0) h <- 1
        W  <- exp(-d2 / (h + 1e-8))
        rs <- rowSums(W)
        rs[!is.finite(rs) | rs <= 0] <- 1
        W  <- W / rs
        p_res_dom <- p_res[pilot_idx]
        eta_cond[syn_idx] <- rowSums(W * matrix(p_res_dom[nn$nn.index], nrow = length(syn_idx)))
      }
      eta_cond[!is.finite(eta_cond)] <- 0
      eta_cond <- eta_cond - mean(eta_cond)
    }
  }

  eta_para <- rep(0, n_cells)
  if (alpha_sp > 0) {
    used <- FALSE
    if (!is.null(spatial_cache)) {
      if (spatial_mode == "original" && exists("draw_graph_gp_whole_tissue", mode = "function")) {
        eta_para <- draw_graph_gp_whole_tissue(spatial_cache, alpha_sp, rho_eff, include_self = TRUE)
        used <- TRUE
      }
      if (spatial_mode == "basis" && exists("draw_basis_gp_whole_tissue", mode = "function")) {
        eta_para <- draw_basis_gp_whole_tissue(spatial_cache, alpha_sp, rho_eff, basis_kernel)
        used <- TRUE
      }
    }
    if (!used) {
      if (n_cells <= 1L) {
        eta_para <- stats::rnorm(n_cells, 0, sqrt(alpha_sp))
      } else {
        D_all <- as.matrix(stats::dist(sim_coords_mat))
        K_all <- alpha_sp * exp(-D_all / rho_eff)
        diag(K_all) <- diag(K_all) + gp_jitter
        eta_para <- low_rank_gp_draw_from_cov(K_all, rank_max = min(as.integer(gp_rank), n_cells), jitter = gp_jitter)
      }
      eta_para <- rescale_spatial_field(eta_para, alpha_sp)
    }
  }
  eta_para[!is.finite(eta_para)] <- 0

  eta <- lambda * eta_cond + (1 - lambda) * eta_para
  eta[!is.finite(eta)] <- 0

  eps_noise <- if (nugget2 > 0) stats::rnorm(n_cells, 0, sqrt(nugget2)) else rep(0, n_cells)

  log1p_y <- mu_log + B_shift + eta + eps_noise
  log1p_y[!is.finite(log1p_y)] <- 0
  log1p_y[log1p_y < 0] <- 0
  return(log1p_y)
}

get_conditional_parameters <- function(mu_full, cov_full, observed_indices, observed_values) {
  mu_obs <- mu_full[observed_indices]
  mu_unobs <- mu_full[-observed_indices]
  Sigma_obs_obs <- cov_full[observed_indices, observed_indices]
  Sigma_unobs_unobs <- cov_full[-observed_indices, -observed_indices]
  Sigma_unobs_obs <- cov_full[-observed_indices, observed_indices]
  diag(Sigma_obs_obs) <- diag(Sigma_obs_obs) + 1e-8

  Sigma_obs_obs_inv <- tryCatch(solve(Sigma_obs_obs), error = function(e) {
    warning(paste("Singular matrix encountered for Sigma_obs_obs:", e$message, "Returning NULL."))
    return(NULL)
  })
  if (is.null(Sigma_obs_obs_inv)) return(list(mu_conditional = NULL, Sigma_conditional_chol = NULL))

  Sigma_uo_Sigma_oo_inv <- Sigma_unobs_obs %*% Sigma_obs_obs_inv
  mu_conditional <- mu_unobs + Sigma_uo_Sigma_oo_inv %*% (observed_values - mu_obs)
  Sigma_conditional <- Sigma_unobs_unobs - Sigma_uo_Sigma_oo_inv %*% t(Sigma_unobs_obs)
  diag(Sigma_conditional) <- diag(Sigma_conditional) + 1e-8

  Sigma_conditional_chol <- tryCatch(chol(Sigma_conditional), error = function(e) {
    warning(paste("Non-positive definite conditional covariance:", e$message, "Returning NULL."))
    return(NULL)
  })
  list(mu_conditional = mu_conditional, Sigma_conditional_chol = Sigma_conditional_chol)
}

generate_fitted_data <- function(mu_est, alpha_spatial_est, rho_est, sigma2_B_est, sigma2_noise_est,
                                 coords, K_syn, Y_list_original, subratio=0.5, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  n <- nrow(coords); p <- length(mu_est); D <- as.matrix(stats::dist(coords))
  if (length(mu_est) == 1) mu_est <- rep(mu_est, p)
  if (length(rho_est) == 1) rho_est <- rep(rho_est, p)
  if (length(alpha_spatial_est) == 1) alpha_spatial_est <- rep(alpha_spatial_est, p)
  if (length(sigma2_B_est) == 1) sigma2_B_est <- rep(sigma2_B_est, p)
  if (length(sigma2_noise_est) == 1) sigma2_noise_est <- rep(sigma2_noise_est, p)

  n_obs <- round(subratio * n)
  observed_indices <- sample(1:n, n_obs)
  unobserved_indices <- setdiff(1:n, observed_indices)
  Y_original_array <- simplify2array(Y_list_original)
  Y_bar_original <- apply(Y_original_array, c(1, 2), mean, na.rm = TRUE)
  Y_bar_sub <- Y_bar_original[observed_indices, ]

  precomputed_cond_params <- lapply(1:p, function(g) {
    Sigma_spatial_full_g <- alpha_spatial_est[g] * exp(-D / rho_est[g])
    inflation_factor <- 1
    Sigma_Y_full_g <- Sigma_spatial_full_g + diag(sigma2_B_est[g] + inflation_factor*sigma2_noise_est[g], n)
    observed_y_bar_g <- Y_bar_sub[, g]
    mu_full_g <- rep(mu_est[g], n)
    get_conditional_parameters(mu_full = mu_full_g, cov_full = Sigma_Y_full_g,
                               observed_indices = observed_indices, observed_values = observed_y_bar_g)
  })

  Y_fitted_list <- lapply(1:K_syn, function(k_syn) {
    Y_k_syn <- matrix(NA, nrow = n, ncol = p)
    set.seed(seed+k_syn)
    Z_unobs_all_genes <- matrix(rnorm(length(unobserved_indices) * p), nrow = length(unobserved_indices), ncol = p)
    for (g in 1:p) {
      current_cond_params <- precomputed_cond_params[[g]]
      if (is.null(current_cond_params$mu_conditional) || is.null(current_cond_params$Sigma_conditional_chol)) {
        Y_k_syn[, g] <- NA
        warning(paste("Skipping gene", g, "due to failed conditional parameter computation."))
      } else {
        unobserved_values_g_sampled <- current_cond_params$mu_conditional +
          t(current_cond_params$Sigma_conditional_chol) %*% Z_unobs_all_genes[, g]
        temp_Y_g <- numeric(n)
        temp_Y_g[observed_indices] <- Y_bar_sub[, g]
        temp_Y_g[unobserved_indices] <- unobserved_values_g_sampled
        Y_k_syn[, g] <- temp_Y_g
      }
    }
    return(Y_k_syn)
  })
  return(Y_fitted_list)
}

estimate_parameters_by_group <- function(Y_list_by_group, groups, n, p, coords, p_B, sigma_B_sq_vec=NULL,
                                         effective_gene_indices, non_effective_gene_indices) {
  if (length(groups) == 0 || length(Y_list_by_group[[as.character(groups[1])]]) == 0) {
    stop("Y_list_by_group is empty or malformed. Cannot determine K.")
  }
  K <- length(Y_list_by_group[[as.character(groups[1])]])
  if (K > 1) {
    result <- estimate_parameters_by_group_multiple_K(Y_list_by_group, groups, n, p, coords, p_B, sigma_B_sq_vec,
                                                      effective_gene_indices, non_effective_gene_indices)
  } else if (K == 1) {
    result <- estimate_parameters_by_group_K1(Y_list_by_group, groups, n, p, coords, p_B,
                                              effective_gene_indices, non_effective_gene_indices)
  } else {
    stop("Invalid value for K (number of pilot samples). K must be >= 1.")
  }
  result$Y_list_original_for_conditional_sampling <- Y_list_by_group
  return(result)
}

estimate_parameters_by_group_multiple_K <- function(Y_list_by_group, groups, n, p, coords, p_B, sigma_B_sq_vec,
                                                    effective_gene_indices, non_effective_gene_indices) {
  results_by_group <- list()
  for (c in seq_along(groups)) {
    group <- groups[c]
    Y_list      <- Y_list_by_group[[group]]
    p_B_group   <- p_B[c]
    sigma_B_sq  <- sigma_B_sq_vec[group]
    K <- length(Y_list)
    brisc_param_array <- array(NA, dim = c(p, 4, K),
                               dimnames = list(paste0("gene", 1:p), c("mu", "alpha_spatial", "rho", "sigma2_noise"), paste0("sample", 1:K)))
    n_neighbors = 10
    task_grid <- expand.grid(g = 1:p, k = 1:K)

    results_list <- lapply(1:nrow(task_grid), function(i) {
      g <- task_grid$g[i]; k <- task_grid$k[i]; Y_gk <- Y_list[[k]][, g]
      if (any(!is.finite(Y_gk))) return(rep(NA, 4))
      if (var(Y_gk, na.rm = TRUE) == 0) Y_gk <- Y_gk + rnorm(length(Y_gk), 0, 1e-6)
      brisc_fit <- suppressMessages(try(BRISC::BRISC_estimation(coords = coords, y = Y_gk, x = matrix(1, nrow = n),
                                                                cov.model = "exponential", n.neighbors = n_neighbors, order = "AMMD",
                                                                verbose = FALSE, nugget_status = 1), silent = TRUE))
      if (!inherits(brisc_fit, "try-error")) {
        return(c(mu = as.numeric(brisc_fit$Beta), alpha_spatial = brisc_fit$Theta["sigma.sq"],
                 rho = brisc_fit$Theta["phi"], sigma2_noise = brisc_fit$Theta["tau.sq"]))
      } else return(rep(NA, 4))
    })

    for (i in 1:nrow(task_grid)) {
      g <- task_grid$g[i]; k <- task_grid$k[i]; brisc_param_array[g, , k] <- results_list[[i]]
    }

    brisc_param_avg <- apply(brisc_param_array, c(1, 2), mean, na.rm = TRUE)
    sigma2_B_data <- apply(brisc_param_array[, "mu", ], 1, var, na.rm = TRUE)
    Y_concat        <- do.call(rbind, Y_list)
    total_var_vec   <- apply(Y_concat, 2, var, na.rm = TRUE)
    sigma2_B_prior  <- p_B_group * total_var_vec
    w_Kp            <- min(1, (K - 1) / 10)
    sigma2_B_hybrid <- (1 - w_Kp) * sigma2_B_prior + w_Kp * sigma2_B_data

    results_by_group[[group]] <- list(
      brisc_param_avg  = brisc_param_avg, sigma2_B_data    = sigma2_B_data,
      sigma2_B_prior   = sigma2_B_prior,  sigma2_B_hybrid  = sigma2_B_hybrid,
      sigma2_B_true    = rep(sigma_B_sq, p)
    )
  }
  return(results_by_group)
}

estimate_parameters_by_group_K1 <- function(Y_list_by_group, groups, n, p, coords, p_B,
                                            effective_gene_indices, non_effective_gene_indices) {
  brisc_results_by_group <- list()
  for (group in groups) {
    Y_mat <- Y_list_by_group[[group]][[1]]
    brisc_param_matrix <- matrix(NA, nrow = p, ncol = 4, dimnames = list(paste0("gene", 1:p), c("mu", "alpha_spatial", "rho", "sigma2_noise")))
    n_neighbors = 15
    for (g in 1:p) {
      Y_g <- Y_mat[, g]
      if (var(Y_g, na.rm = TRUE) == 0) Y_g <- Y_g + rnorm(length(Y_g), 0, 1e-6)
      brisc_fit <- suppressMessages(try(BRISC::BRISC_estimation(coords = coords, y = Y_g, x = matrix(1, nrow = n),
                                                                cov.model = "exponential", n.neighbors = n_neighbors, order = "AMMD",
                                                                verbose = FALSE, nugget_status = 1), silent = TRUE))
      if (!inherits(brisc_fit, "try-error")) {
        brisc_param_matrix[g, ] <- c(as.numeric(brisc_fit$Beta), brisc_fit$Theta["sigma.sq"],
                                     brisc_fit$Theta["phi"], brisc_fit$Theta["tau.sq"])
      }
    }
    brisc_results_by_group[[group]] <- brisc_param_matrix
  }
  alpha_hat_pooled <- (brisc_results_by_group$case[, "alpha_spatial"] + brisc_results_by_group$control[, "alpha_spatial"]) / 2
  rho_hat_pooled   <- (brisc_results_by_group$case[, "rho"] + brisc_results_by_group$control[, "rho"]) / 2

  final_results <- list()
  sigma2_B_hybrid_by_group <- list()
  for(c_idx in seq_along(groups)){
    group <- groups[c_idx]
    estimated_nuggets_g <- brisc_results_by_group[[group]][, "sigma2_noise"]
    sigma2_B_hybrid_by_group[[group]] <- mean(p_B[c_idx] * estimated_nuggets_g, na.rm = TRUE)
  }
  for (c_idx in seq_along(groups)) {
    group <- groups[c_idx]
    brisc_param_avg <- brisc_results_by_group[[group]]
    brisc_param_avg[, "alpha_spatial"] <- alpha_hat_pooled
    brisc_param_avg[, "rho"] <- rho_hat_pooled
    final_results[[group]] <- list(brisc_param_avg = brisc_param_avg, sigma2_B_hybrid = sigma2_B_hybrid_by_group[[group]])
  }
  return(final_results)
}

get_param_table <- function(results_by_group, groups, p) {
  final_param_estimates <- do.call(rbind, lapply(groups, function(group) {
    brisc_param_avg <- results_by_group[[group]]$brisc_param_avg
    sigma2_B_hybrid <- results_by_group[[group]]$sigma2_B_hybrid
    data.frame(Gene = paste0("Gene", 1:p), mu_est = brisc_param_avg[, "mu"],
               alpha_spatial_est = brisc_param_avg[, "alpha_spatial"], rho_est = brisc_param_avg[, "rho"],
               sigma2_B_hybrid = sigma2_B_hybrid, sigma2_noise_est = brisc_param_avg[, "sigma2_noise"], Group = group)
  }))
  rownames(final_param_estimates) <- NULL
  return(final_param_estimates)
}

get_sample_means <- function(Y_list) {
  K <- length(Y_list); p <- ncol(Y_list[[1]])
  means <- matrix(NA, nrow = p, ncol = K)
  for (k in 1:K) means[, k] <- colMeans(Y_list[[k]])
  rownames(means) <- paste0("Gene", 1:p); colnames(means) <- paste0("Sample", 1:K)
  means
}

get_pooled_variances <- function(Y_list) {
  s2_matrix <- do.call(rbind, lapply(Y_list, function(mat) apply(mat, 2, stats::var, na.rm = TRUE)))
  colMeans(s2_matrix, na.rm = TRUE)
}

make_alternative_from_true_beta <- function(beta_mat, tol = 0) {
  if (!all(c("case","control") %in% colnames(beta_mat))) stop("beta_mat must have columns 'case' and 'control'.")
  diff <- beta_mat[, "case"] - beta_mat[, "control"]
  alt  <- character(length(diff))
  alt[diff >  tol] <- "greater"
  alt[diff < -tol] <- "less"
  alt[abs(diff) <= tol] <- "two.sided"
  alt
}

moderated_ttest <- function(Y_case, Y_control, gene_effective = NULL, alpha = 0.05, alternative = "auto", direction = "greater") {
  .match_alt <- function(x) { x <- tolower(x); if (!x %in% c("greater", "less", "two.sided")) stop("invalid alternative"); x }
  case_means    <- get_sample_means(Y_case)
  control_means <- get_sample_means(Y_control)
  p <- nrow(case_means); K_case <- ncol(case_means); K_control <- ncol(control_means)
  s2_case_all    <- get_pooled_variances(Y_case)
  s2_control_all <- get_pooled_variances(Y_control)

  bar_s2_case    <- mean(s2_case_all, na.rm = TRUE); bar_s2_control <- mean(s2_control_all, na.rm = TRUE)
  v_case         <- stats::var(s2_case_all, na.rm = TRUE); v_control <- stats::var(s2_control_all, na.rm = TRUE)
  d0_case        <- if (!is.na(v_case) && v_case > 0) 2 * bar_s2_case^2 / v_case else 10
  d0_control     <- if (!is.na(v_control) && v_control > 0) 2 * bar_s2_control^2 / v_control else 10

  if (length(alternative) == 1L) {
    if (tolower(alternative) == "auto") {
      alt_vec <- rep("two.sided", p)
      if (!is.null(gene_effective) && length(gene_effective)) {
        idx <- as.integer(gene_effective); idx <- idx[idx >= 1 & idx <= p]; alt_vec[idx] <- .match_alt(direction)
      }
    } else alt_vec <- rep(.match_alt(alternative), p)
  } else {
    alt_vec <- tolower(alternative); if (length(alt_vec) != p || !all(alt_vec %in% c("greater","less","two.sided"))) stop("invalid per-gene alternative")
  }

  results <- lapply(seq_len(p), function(g) {
    var_case <- s2_case_all[g]; var_control <- s2_control_all[g]
    s2_tilde_case    <- (d0_case * bar_s2_case    + (K_case    - 1) * var_case   ) / (d0_case    + K_case    - 1)
    s2_tilde_control <- (d0_control * bar_s2_control + (K_control - 1) * var_control) / (d0_control + K_control - 1)
    mean_case_g      <- mean(case_means[g, ]); mean_control_g  <- mean(control_means[g, ])
    mean_diff        <- mean_case_g - mean_control_g
    denom            <- sqrt(s2_tilde_case / K_case + s2_tilde_control / K_control)

    if (!is.finite(denom) || denom <= 0) return(c(mean_case = mean_case_g, mean_control = mean_control_g, mean_diff = mean_diff, sd = denom, t_moderated = NA_real_, p_moderated = NA_real_, df_moderated = NA_real_))

    t_mod <- mean_diff / denom
    denom_df <- ((s2_tilde_case / K_case)^2 / (d0_case + K_case - 1)) + ((s2_tilde_control / K_control)^2 / (d0_control + K_control - 1))
    df_mod <- ((s2_tilde_case / K_case + s2_tilde_control / K_control)^2) / denom_df
    p_mod <- if (alt_vec[g] == "greater") stats::pt(t_mod, df = df_mod, lower.tail = FALSE) else if (alt_vec[g] == "less") stats::pt(t_mod, df = df_mod, lower.tail = TRUE) else 2 * stats::pt(-abs(t_mod), df = df_mod)
    c(mean_case = mean_case_g, mean_control = mean_control_g, mean_diff = mean_diff, sd = denom, t_moderated = t_mod, p_moderated = p_mod, df_moderated = df_mod)
  })

  t_test_results <- as.data.frame(do.call(rbind, results))
  t_test_results <- cbind(Gene = paste0("Gene", seq_len(p)), Alternative = alt_vec, t_test_results)
  t_test_results$q_value    <- stats::p.adjust(t_test_results$p_moderated, method = "BH")
  t_test_results$p_value_bf <- stats::p.adjust(t_test_results$p_moderated, method = "bonferroni")
  t_test_results$significant_fdr  <- t_test_results$q_value    < alpha
  t_test_results$significant_fwer <- t_test_results$p_value_bf < alpha
  t_test_results
}

run_wilcoxon_test_per_gene <- function(case_means, control_means,
                                       gene_effective = NULL) {
  p <- nrow(case_means); p_values <- numeric(p)
  for (g in seq_len(p)) {
    x_case <- as.numeric(case_means[g, ]); x_control <- as.numeric(control_means[g, ])
    if (stats::var(x_case, na.rm = TRUE) == 0 &&
        stats::var(x_control, na.rm = TRUE) == 0) { p_values[g] <- 1.0; next }
    alt_hypothesis <- if (!is.null(gene_effective) && g %in% gene_effective)
      "greater" else "two.sided"
    p_values[g] <- stats::wilcox.test(x_case, x_control,
                                      alternative = alt_hypothesis)$p.value
  }
  p_values
}

compute_directional_neighbors <- function(loc_df, radius = 0.1) {
  coords <- as.matrix(loc_df[, c("x", "y")]); domain <- as.character(loc_df$domain); unique_domains <- sort(unique(domain))
  domain_neighbors <- stats::setNames(vector("list", length(unique_domains)), unique_domains)
  for (d1 in unique_domains) {
    idx1 <- which(domain == d1)
    for (d2 in unique_domains) {
      if (d1 == d2) next
      idx2 <- which(domain == d2)
      dist_mat <- as.matrix(stats::dist(rbind(coords[idx1, , drop=FALSE], coords[idx2, , drop=FALSE])))
      dist_sub <- dist_mat[1:length(idx1), (length(idx1)+1):(length(idx1)+length(idx2))]
      if (any(dist_sub <= radius)) domain_neighbors[[d1]] <- c(domain_neighbors[[d1]], d2)
    }
  }
  for (d in unique_domains) domain_neighbors[[d]] <- sort(unique(domain_neighbors[[d]]))
  domain_neighbors
}

compute_domain_neighbors <- function(loc_df, radius = 0.1) {
  coords <- as.matrix(loc_df[, c("x", "y")]); domain <- as.character(loc_df$domain); unique_domains <- sort(unique(domain))
  domain_neighbors <- stats::setNames(vector("list", length(unique_domains)), unique_domains)
  for (d1 in unique_domains) {
    idx1 <- which(domain == d1); pts1 <- coords[idx1, , drop = FALSE]; neighbor_set <- c()
    for (d2 in unique_domains) {
      if (d1 == d2) next
      idx2 <- which(domain == d2); pts2 <- coords[idx2, , drop = FALSE]
      dist_mat <- as.matrix(stats::dist(rbind(pts1, pts2)))
      dist_sub <- dist_mat[1:nrow(pts1), (nrow(pts1)+1):(nrow(pts1)+nrow(pts2))]
      if (any(dist_sub <= radius)) neighbor_set <- c(neighbor_set, d2)
    }
    domain_neighbors[[d1]] <- sort(unique(neighbor_set))
  }
  domain_neighbors
}

generate_pattern <- function(location_df, perturb_domains = NULL, perturb_probs = NULL, seed = 123, boundary_radius = 0.05, boundary_threshold = 0.5) {
  set.seed(seed)
  loc_df <- location_df; coords <- as.matrix(loc_df[, c("x", "y")]); domain <- as.character(loc_df$domain); original_domains <- sort(unique(domain))
  domain_neighbors <- compute_directional_neighbors(loc_df, radius = boundary_radius)
  domain_centroids <- stats::aggregate(cbind(x, y) ~ domain, data = loc_df, FUN = stats::median)
  for (d_char in sort(original_domains, decreasing = TRUE)) {
    if (!is.null(perturb_domains) && !(d_char %in% perturb_domains)) next
    idx_d <- which(domain == d_char); n_d <- length(idx_d)
    prob_d <- if (!is.null(perturb_probs)) perturb_probs[[as.character(d_char)]] else 0.05
    if (is.na(prob_d) || prob_d <= 0) next
    n_change <- floor(n_d * prob_d); neighbors_of_d <- domain_neighbors[[d_char]]
    if (length(neighbors_of_d) == 0) next
    is_boundary_point <- function(i) {
      pt <- matrix(coords[i, ], nrow = 1); dists <- sqrt(rowSums(sweep(coords, 2, pt)^2)); nearby_idx <- which(dists > 0 & dists <= boundary_radius)
      if (length(nearby_idx) == 0) return(FALSE)
      mean(domain[nearby_idx] %in% neighbors_of_d) >= boundary_threshold
    }
    boundary_idx <- idx_d[sapply(idx_d, is_boundary_point)]
    if (length(boundary_idx) == 0) next
    boundary_dists_matrix <- sapply(neighbors_of_d, function(target_d) {
      sapply(boundary_idx, function(i) { pt <- matrix(coords[i, ], nrow = 1); target_coords <- coords[domain == target_d, , drop = FALSE]; min(sqrt(rowSums(sweep(target_coords, 2, pt)^2)), na.rm = TRUE) })
    })
    min_dists <- apply(boundary_dists_matrix, 1, min); sorted_idx <- boundary_idx[order(min_dists)]; change_idx <- utils::head(sorted_idx, n_change)
    if (length(change_idx) > 0) {
      neighbor_centroids <- domain_centroids[domain_centroids$domain %in% neighbors_of_d, ]
      for (i in change_idx) {
        point_coords <- coords[i, ]
        distances_to_neighbors <- sqrt(rowSums(sweep(neighbor_centroids[, c("x", "y")], 2, point_coords)^2))
        domain[i] <- as.character(neighbor_centroids$domain[which.min(distances_to_neighbors)])
      }
    }
  }
  loc_df$domain <- factor(domain, levels = original_domains); loc_df <- loc_df[order(loc_df$domain), ]; rownames(loc_df) <- NULL; loc_df
}

estimate_dirichlet_alpha_mom <- function(counts, epsilon=1e-2) {
  props <- (counts + epsilon) / rowSums(counts + epsilon * ncol(counts))
  pi_bar <- colMeans(props); v_bar <- matrixStats::colVars(props); v_bar[v_bar < 1e-6] <- 1e-6
  alpha0_ests <- (pi_bar * (1 - pi_bar)) / v_bar - 1
  alpha0 <- mean(alpha0_ests, na.rm = TRUE)
  alpha0 <- min(max(alpha0, 1e-2), 100)
  alpha0 * pi_bar
}

dm_loglik <- function(alpha, counts) {
  K <- nrow(counts); D <- ncol(counts); alpha0 <- sum(alpha)
  term1 <- K * lgamma(alpha0); term2 <- sum(lgamma(rowSums(counts) + alpha0))
  term3 <- sum(lgamma(sweep(counts, 2, alpha, "+"))); term4 <- K * sum(lgamma(alpha))
  term1 - term2 + term3 - term4
}

run_global_composition_test <- function(counts_case, counts_ctrl) {
  D <- ncol(counts_case); counts_pooled <- rbind(counts_case, counts_ctrl)
  alpha_null <- estimate_dirichlet_alpha_mom(counts_pooled); ll_null <- dm_loglik(alpha_null, counts_pooled)
  alpha_case <- estimate_dirichlet_alpha_mom(counts_case); alpha_ctrl <- estimate_dirichlet_alpha_mom(counts_ctrl)
  ll_alt <- dm_loglik(alpha_case, counts_case) + dm_loglik(alpha_ctrl, counts_ctrl)
  lrt_stat <- 2 * (ll_alt - ll_null); p_val <- stats::pchisq(lrt_stat, df = D, lower.tail = FALSE)
  list(statistic = lrt_stat, p_value = p_val, df = D)
}

extract_geometry_features <- function(sample_list, D) {
  K <- length(sample_list); R_mat <- matrix(NA, K, D); A_mat <- matrix(NA, K, D)
  for (k in 1:K) {
    df <- sample_list[[k]]
    if(!"domain" %in% colnames(df)) next
    for (d in 1:D) {
      coords <- df[df$domain == d, c("x", "y")]
      if (nrow(coords) < 5) next
      Sigma <- stats::cov(coords); eig <- eigen(Sigma)$values; lambda1 <- max(eig); lambda2 <- max(min(eig), 1e-6)
      R_mat[k, d] <- 0.5 * log(lambda1 + lambda2); A_mat[k, d] <- log(lambda1 / lambda2)
    }
  }
  list(R = R_mat, A = A_mat)
}

estimate_prior_params_mle <- function(sample_variances, df1) {
  sample_variances <- sample_variances[!is.na(sample_variances)]
  neg_log_likelihood <- function(params) {
    s0_sq <- params[1]; d0 <- params[2]
    if (s0_sq <= 0 || d0 <= 4.0001) return(1e100)
    vals <- sample_variances / s0_sq
    log_f <- stats::df(vals, df1 = df1, df2 = d0, log = TRUE)
    log_lik_terms <- log_f - log(s0_sq)
    if (any(!is.finite(log_lik_terms))) return(1e100)
    -sum(log_lik_terms)
  }
  initial_s0_sq <- mean(sample_variances); v_s2 <- stats::var(sample_variances)
  initial_d0 <- if (is.na(v_s2) || v_s2 <= 1e-8) 20 else min(max(2 * initial_s0_sq^2 / v_s2 + 4, 5), 50)
  if (initial_s0_sq <= 0) initial_s0_sq <- 1e-3
  initial_params <- c(initial_s0_sq, initial_d0)
  opt_result <- tryCatch({ stats::optim(par = initial_params, fn = neg_log_likelihood, method = "L-BFGS-B", lower = c(1e-6, 4.1), upper = c(Inf, 1e5)) }, error = function(e) list(par = initial_params))
  list(s0_sq = opt_result$par[1], d0 = opt_result$par[2])
}

run_joint_geometry_test <- function(feat_case, feat_ctrl, alpha = 0.05) {
  n_domains <- ncol(feat_case$R)
  results_df <- data.frame(domain = 1:n_domains, p_value = NA, T_R = NA, T_A = NA, Q = NA)
  calc_stats <- function(mat) list(mean = colMeans(mat, na.rm=TRUE), var = matrixStats::colVars(mat, na.rm=TRUE), n = colSums(!is.na(mat)))
  stats_R_case <- calc_stats(feat_case$R); stats_R_ctrl <- calc_stats(feat_ctrl$R)
  stats_A_case <- calc_stats(feat_case$A); stats_A_ctrl <- calc_stats(feat_ctrl$A)
  prior_R_case <- estimate_prior_params_mle(stats_R_case$var, stats::median(stats_R_case$n)-1)
  prior_R_ctrl <- estimate_prior_params_mle(stats_R_ctrl$var, stats::median(stats_R_ctrl$n)-1)
  prior_A_case <- estimate_prior_params_mle(stats_A_case$var, stats::median(stats_A_case$n)-1)
  prior_A_ctrl <- estimate_prior_params_mle(stats_A_ctrl$var, stats::median(stats_A_ctrl$n)-1)
  for (d in 1:n_domains) {
    df1 <- stats_R_case$n[d] - 1; df2 <- stats_R_ctrl$n[d] - 1
    if(df1 < 1 || df2 < 1) next
    s2_R_case_mod <- (prior_R_case$d0 * prior_R_case$s0_sq + df1 * stats_R_case$var[d]) / (prior_R_case$d0 + df1)
    s2_R_ctrl_mod <- (prior_R_ctrl$d0 * prior_R_ctrl$s0_sq + df2 * stats_R_ctrl$var[d]) / (prior_R_ctrl$d0 + df2)
    se_R <- sqrt(s2_R_case_mod/stats_R_case$n[d] + s2_R_ctrl_mod/stats_R_ctrl$n[d])
    t_R <- (stats_R_case$mean[d] - stats_R_ctrl$mean[d]) / se_R
    df1_A <- stats_A_case$n[d] - 1; df2_A <- stats_A_ctrl$n[d] - 1
    s2_A_case_mod <- (prior_A_case$d0 * prior_A_case$s0_sq + df1_A * stats_A_case$var[d]) / (prior_A_case$d0 + df1_A)
    s2_A_ctrl_mod <- (prior_A_ctrl$d0 * prior_A_ctrl$s0_sq + df2_A * stats_A_ctrl$var[d]) / (prior_A_ctrl$d0 + df2_A)
    se_A <- sqrt(s2_A_case_mod/stats_A_case$n[d] + s2_A_ctrl_mod/stats_A_ctrl$n[d])
    t_A <- (stats_A_case$mean[d] - stats_A_ctrl$mean[d]) / se_A
    Q <- t_R^2 + t_A^2
    results_df$p_value[d] <- stats::pchisq(Q, df = 2, lower.tail = FALSE)
    results_df$T_R[d] <- t_R; results_df$T_A[d] <- t_A; results_df$Q[d] <- Q
  }
  results_df$q_value <- stats::p.adjust(results_df$p_value, method = "BH")
  results_df$rejected <- results_df$q_value < alpha
  results_df
}

process_single_visium <- function(folder_path, group_label, domain_col="benmarklabel") {
  sample_id <- basename(folder_path)
  message(sprintf("Processing Sample: %s (Group: %s)...", sample_id, group_label))
  rdata_path <- file.path(folder_path, "brain_processed.RData")
  meta_csv_path <- file.path(folder_path, "meta.csv")
  if(!file.exists(rdata_path) || !file.exists(meta_csv_path)) { warning(paste("Missing files in:", folder_path)); return(NULL) }
  e <- new.env(); load(rdata_path, envir = e); obj_name <- ls(e)[1]; brain <- e[[obj_name]]
  if(!inherits(brain, "Seurat")) { warning(paste("Object in", rdata_path, "is not a Seurat object.")); return(NULL) }
  assay_use <- "Spatial"
  if(!"Spatial" %in% names(brain@assays)) assay_use <- Seurat::DefaultAssay(brain)
  count_matrix <- Seurat::GetAssayData(brain, assay = assay_use, layer = "counts")
  if(length(brain@images) == 0) { warning(paste("No images found in Seurat object for sample", sample_id)); return(NULL) }
  image_name <- names(brain@images)[1]; coords <- brain@images[[image_name]]@coordinates
  loc0 <- data.frame(x = as.numeric(coords$col), y = as.numeric(coords$row), stringsAsFactors = FALSE); rownames(loc0) <- rownames(coords)
  x_min <- min(loc0$x); x_max <- max(loc0$x); if (x_max > x_min) loc0$x <- (loc0$x - x_min) / (x_max - x_min)
  y_min <- min(loc0$y); y_max <- max(loc0$y); if (y_max > y_min) loc0$y <- (loc0$y - y_min) / (y_max - y_min)
  meta <- readr::read_csv(meta_csv_path, show_col_types = FALSE); barcode_col_csv <- "barcode"
  if (!all(c(barcode_col_csv, domain_col) %in% colnames(meta))) stop("meta.csv is missing required columns: barcode or domain column.")
  meta$barcode <- meta[[barcode_col_csv]]; meta$domain <- meta[[domain_col]]
  meta$spot_id_plain <- as.character(meta$barcode); meta$spot_id_prefix <- paste0(sample_id, "_", as.character(meta$barcode))
  cn <- colnames(count_matrix); use_prefixed <- FALSE
  if (all(grepl("^\\d+_", cn))) use_prefixed <- TRUE
  meta_spot_id <- if (use_prefixed) meta$spot_id_prefix else meta$spot_id_plain
  domain_map <- stats::setNames(as.character(meta$domain), meta_spot_id)
  if (use_prefixed && !all(grepl("^\\d+_", rownames(loc0)))) rownames(loc0) <- paste0(sample_id, "_", rownames(loc0))
  loc0$domain <- domain_map[rownames(loc0)]
  common <- intersect(colnames(count_matrix), rownames(loc0))
  if (length(common) == 0) { warning(paste("No matching spots for sample", sample_id)); return(NULL) }
  count_matrix <- count_matrix[, common, drop = FALSE]; loc_file <- loc0[common, , drop = FALSE]
  keep <- !is.na(loc_file$domain); count_matrix <- count_matrix[, keep, drop = FALSE]; loc_file <- loc_file[keep, , drop = FALSE]
  loc_file$domain <- factor(loc_file$domain)
  if(!identical(colnames(count_matrix), rownames(loc_file))) stop("Mismatch between count matrix columns and coordinate rows.")
  log_matrix <- log1p(count_matrix)
  list(counts = count_matrix, logcounts = log_matrix, coords = loc_file, group = group_label, sample_id = sample_id)
}

get_domain_markers <- function(object, top_n = 2) {
  rep_sample <- object@pilot_data[[1]]; count_mat <- rep_sample$counts; coords <- rep_sample$coords
  log_counts <- log1p(as.matrix(count_mat)); domains <- sort(unique(coords$domain)); marker_list <- list()
  message(">>> Identifying top ", top_n, " markers per domain...")
  for(dom in domains) {
    in_idx  <- which(coords$domain == dom); mean_in <- rowMeans(log_counts[, in_idx, drop=FALSE]); mean_out <- rowMeans(log_counts[, -in_idx, drop=FALSE])
    lfc <- mean_in - mean_out; df <- data.frame(gene = rownames(log_counts), lfc = lfc)
    df <- df[order(df$lfc, decreasing = TRUE), ]
    top_genes <- utils::head(df$gene, top_n)
    marker_list[[dom]] <- top_genes
    cat(sprintf("   Domain %s: %s\n", dom, paste(top_genes, collapse=", ")))
  }
  unique(unlist(marker_list))
}

build_z_from_pilot <- function(object, genes, target_domain="WM", reference_domain="Layer6") {
  pilot <- object@pilot_data; rows <- list(); k <- 0L
  for (i in seq_along(pilot)) {
    samp <- pilot[[i]]; coords <- samp$coords; X <- samp$logcounts; grp <- as.numeric(samp$group)
    idx_T <- which(coords$domain == target_domain); idx_R <- which(coords$domain == reference_domain)
    if (length(idx_T) < 2 || length(idx_R) < 2) next
    for (g in genes) {
      if (!g %in% rownames(X)) next
      z <- mean(as.numeric(X[g, idx_T])) - mean(as.numeric(X[g, idx_R]))
      k <- k + 1L; rows[[k]] <- data.frame(gene=g, sample=i, group=grp, z_hat=z)
    }
  }
  do.call(rbind, rows)
}

make_TR_pool <- function(object, target_domain="WM", reference_domain="Layer6", sample_index=1, mean_cut=0.5, pct_cut=0.05) {
  samp <- object@pilot_data[[sample_index]]; X <- samp$logcounts; coords <- samp$coords
  idx_T <- which(coords$domain == target_domain); idx_R <- which(coords$domain == reference_domain)
  stopifnot(length(idx_T) > 0, length(idx_R) > 0)
  mean_T <- rowMeans(X[, idx_T, drop=FALSE]); mean_R <- rowMeans(X[, idx_R, drop=FALSE])
  pct_T  <- rowMeans(X[, idx_T, drop=FALSE] > 0); pct_R  <- rowMeans(X[, idx_R, drop=FALSE] > 0)
  df <- data.frame(gene = rownames(X), mean_T = mean_T, mean_R = mean_R, pct_T = pct_T, pct_R = pct_R, TR = mean_T - mean_R, stringsAsFactors = FALSE)
  pool <- df$gene[df$mean_T >= mean_cut & df$mean_R >= mean_cut & df$pct_T >= pct_cut & df$pct_R >= pct_cut]
  list(pool=pool, stats=df)
}

make_gene_sets_pilot_only <- function(object, target_domain="WM", reference_domain="Layer6", mean_cut=0.5, pct_cut=0.05, p_exclude = 0.2, de_tr_q = c(0.4, 0.8), n_de = 3, n_test = 10, verbose=TRUE) {
  TR <- make_TR_pool(object, target_domain, reference_domain, sample_index=1, mean_cut=mean_cut, pct_cut=pct_cut); pool <- TR$pool; st <- TR$stats
  z_long <- build_z_from_pilot(object, pool, target_domain, reference_domain)
  p_by_gene <- tapply(seq_len(nrow(z_long)), z_long$gene, function(ii) {
    df <- z_long[ii, ]; if (length(unique(df$group)) < 2) return(c(p_t=NA_real_, p_w=NA_real_))
    p_t <- tryCatch(stats::t.test(z_hat ~ factor(group), data=df)$p.value, error=function(e) NA_real_)
    p_w <- tryCatch(stats::wilcox.test(z_hat ~ factor(group), data=df, exact=FALSE)$p.value, error=function(e) NA_real_)
    c(p_t=p_t, p_w=p_w)
  })
  p_mat <- do.call(rbind, p_by_gene)
  keep <- rownames(p_mat)[is.finite(p_mat[,"p_t"]) & is.finite(p_mat[,"p_w"]) & p_mat[,"p_t"] >= p_exclude & p_mat[,"p_w"] >= p_exclude]
  pool2 <- intersect(pool, keep)
  if (verbose) { message("TR pool size: ", length(pool)); message("After excluding suspicious group-diff genes (p<", p_exclude, "): ", length(pool2)) }
  if (length(pool2) < n_test) warning("pool2 is small; relax mean_cut/pct_cut or p_exclude.")
  st2 <- st[st$gene %in% pool2, , drop=FALSE]; absTR <- abs(st2$TR)
  qlo <- stats::quantile(absTR, probs=de_tr_q[1], na.rm=TRUE); qhi <- stats::quantile(absTR, probs=de_tr_q[2], na.rm=TRUE)
  de_pool <- st2$gene[absTR >= qlo & absTR <= qhi]; de_pool <- de_pool[order(-abs(st2$TR[match(de_pool, st2$gene)]))]
  G_DE <- utils::head(de_pool, n_de); remain <- setdiff(pool2, G_DE)
  remain <- remain[order(abs(st2$TR[match(remain, st2$gene)]))]
  G_null <- utils::head(remain, max(0, n_test - length(G_DE))); G_test <- c(G_DE, G_null)
  list(G_TR_pool = pool, G_test = G_test, G_DE = G_DE, TR_stats = st2, p_groupdiff = p_by_gene)
}

.safe_unit_vec <- function(v, fallback = c(1, 0)) {
  fallback <- as.numeric(fallback); v <- as.numeric(v)
  if (length(v) != length(fallback) || any(!is.finite(v))) return(fallback)
  nv <- sqrt(sum(v^2)); if (!is.finite(nv) || nv <= 1e-12) return(fallback)
  v / nv
}

.safe_mvrnorm_vec <- function(mu, Sigma, jitter = 1e-8) {
  mu <- as.numeric(mu); p <- length(mu)
  if (p == 0L || any(!is.finite(mu))) { mu <- rep(0, max(p, 2L)); p <- length(mu) }
  Sigma <- as.matrix(Sigma)
  if (!all(dim(Sigma) == c(p, p)) || any(!is.finite(Sigma))) Sigma <- diag(jitter, p)
  Sigma <- (Sigma + t(Sigma)) / 2; diag(Sigma) <- diag(Sigma) + jitter
  out <- tryCatch(MASS::mvrnorm(1, mu = mu, Sigma = Sigma), error = function(e) mu)
  as.numeric(out)
}

.safe_rvmf1 <- function(mu, kappa, fallback = c(1, 0)) {
  mu <- .safe_unit_vec(mu, fallback = fallback)
  if (!is.finite(kappa) || kappa <= 1e-8) { ang <- stats::runif(1, 0, 2 * pi); return(c(cos(ang), sin(ang))) }
  out <- tryCatch(as.numeric(Directional::rvmf(1, mu, kappa)), error = function(e) mu)
  .safe_unit_vec(out, fallback = fallback)
}

.fallback_gaussian_density <- function(grid_xy, centroid, scale = 0.08, eps = 1e-12) {
  Z <- sweep(as.matrix(grid_xy), 2, as.numeric(centroid), "-", check.margin = FALSE)
  dens <- exp(-rowSums(Z^2) / max(scale, eps))
  dens[!is.finite(dens)] <- eps
  pmax(dens, eps)
}

.count_labels <- function(labels, doms) {
  tab <- tabulate(match(labels, doms), nbins = length(doms))
  stats::setNames(as.integer(tab), doms)
}

.sample_domain_centroid <- function(geom_p, fallback = c(0.5, 0.5)) {
  if (is.null(geom_p) || is.null(geom_p$placement)) return(as.numeric(fallback))
  mu_G <- geom_p$placement$mu_G; Sigma_G <- geom_p$placement$Sigma_G
  if (is.null(mu_G) || length(mu_G) != 2L || any(!is.finite(mu_G))) mu_G <- as.numeric(fallback)
  if (is.null(Sigma_G) || !all(dim(as.matrix(Sigma_G)) == c(2L, 2L))) Sigma_G <- diag(1e-4, 2L)
  .safe_mvrnorm_vec(mu = mu_G, Sigma = Sigma_G, jitter = 1e-8)
}

sample_one_component_from_global <- function(g, lambda_cloud = 0.7, fallback_phi = c(1, 0), eps = 1e-8) {
  c_anchor <- if (!is.null(g$center_cloud) && length(g$center_cloud) > 0L) as.numeric(g$center_cloud[[sample.int(length(g$center_cloud), 1L)]]) else as.numeric(g$mu_c)
  phi_anchor <- if (!is.null(g$orientation_cloud) && length(g$orientation_cloud) > 0L) .safe_unit_vec(g$orientation_cloud[[sample.int(length(g$orientation_cloud), 1L)]], fallback = fallback_phi) else .safe_unit_vec(g$mu_phi, fallback = fallback_phi)
  c_prior <- .safe_mvrnorm_vec(g$mu_c, g$Sigma_c, jitter = eps)
  r_prior <- abs(stats::rnorm(1, mean = as.numeric(g$mu_r), sd = sqrt(max(as.numeric(g$var_r), eps))))
  phi_prior <- .safe_rvmf1(mu = g$mu_phi, kappa = as.numeric(g$kappa_phi), fallback = fallback_phi)
  tau_prior <- stats::rgamma(1, shape = max(as.numeric(g$a_tau), eps), rate = max(as.numeric(g$b_tau), eps)); tau_prior <- max(tau_prior, 1e-6)
  c_draw <- lambda_cloud * c_anchor + (1 - lambda_cloud) * c_prior
  phi_draw_raw <- lambda_cloud * phi_anchor + (1 - lambda_cloud) * phi_prior
  phi_draw <- .safe_unit_vec(phi_draw_raw, fallback = fallback_phi)
  list(c = as.numeric(c_draw), r = as.numeric(r_prior), phi = as.numeric(phi_draw), tau = as.numeric(tau_prior))
}

sample_shared_sigma_sq_from_pooled <- function(global_theta, eps = 1e-8) {
  a_vec <- vapply(global_theta, function(g) as.numeric(g$a_sigma), numeric(1))
  b_vec <- vapply(global_theta, function(g) as.numeric(g$b_sigma), numeric(1))
  a_bar <- mean(a_vec[is.finite(a_vec) & a_vec > 0], na.rm = TRUE)
  b_bar <- mean(b_vec[is.finite(b_vec) & b_vec > 0], na.rm = TRUE)
  if (!is.finite(a_bar) || a_bar <= 0) a_bar <- 10
  if (!is.finite(b_bar) || b_bar <= 0) b_bar <- 10
  prec <- stats::rgamma(1, shape = a_bar, rate = b_bar); prec <- max(as.numeric(prec), eps)
  1 / prec
}

sample_fgkmm_domain_from_pooled <- function(pooled, lambda_cloud = 0.7, eps = 1e-8) {
  if (is.null(pooled) || is.null(pooled$global_theta) || length(pooled$global_theta) == 0L) stop("pooled must contain non-empty global_theta.")
  M <- length(pooled$global_theta)
  alpha_hat <- pooled$alpha_hat
  if (is.null(alpha_hat) || length(alpha_hat) != M || any(!is.finite(alpha_hat))) alpha_hat <- rep(1, M)
  alpha_hat <- pmax(as.numeric(alpha_hat), eps)
  pi_draw <- as.numeric(MCMCpack::rdirichlet(1, alpha_hat)); pi_draw <- pi_draw / sum(pi_draw)
  comp_list <- lapply(pooled$global_theta, function(g) sample_one_component_from_global(g = g, lambda_cloud = lambda_cloud, fallback_phi = c(1, 0), eps = eps))
  c_hat <- do.call(rbind, lapply(comp_list, `[[`, "c")); phi_hat <- do.call(rbind, lapply(comp_list, `[[`, "phi"))
  r_hat <- vapply(comp_list, `[[`, numeric(1), "r"); tau_hat <- vapply(comp_list, `[[`, numeric(1), "tau")
  sigma_sq <- sample_shared_sigma_sq_from_pooled(pooled$global_theta, eps = eps)
  list(pi = pi_draw, c_hat = c_hat, r_hat = r_hat, phi_hat = phi_hat, tau_hat = tau_hat, sigma_sq = sigma_sq, component_list = comp_list)
}

evaluate_fg_mixture_density <- function(query_xy, centroid, fg_draw, eps = 1e-300) {
  X <- as.matrix(query_xy); n <- nrow(X); if (n == 0L) return(numeric(0))
  centroid <- as.numeric(centroid); if (length(centroid) != ncol(X)) stop("centroid dimension mismatch.")
  X_centered <- sweep(X, 2, centroid, "-", check.margin = FALSE); M <- length(fg_draw$pi)
  ll_mat <- matrix(NA_real_, nrow = n, ncol = M)
  for (m in seq_len(M)) {
    ll_mat[, m] <- log(pmax(fg_draw$pi[m], eps)) + log_fg_one_comp(X = X_centered, c_k = fg_draw$c_hat[m, ], r_k = fg_draw$r_hat[m], sigma_sq = fg_draw$sigma_sq, phi_k = fg_draw$phi_hat[m, ], tau_k = fg_draw$tau_hat[m])
  }
  dens <- exp(matrixStats::rowLogSumExps(ll_mat)); dens[!is.finite(dens)] <- eps
  pmax(dens, eps)
}

# Helper: Inverse warp query for geometric kappa
.inverse_warp_query <- function(query_xy, centroid, kappa_geo) {
  kg <- as.numeric(kappa_geo)
  if (!is.finite(kg) || kg <= 0) kg <- 1

  centroid <- matrix(as.numeric(centroid), nrow = nrow(query_xy), ncol = 2L, byrow = TRUE)
  centroid + (as.matrix(query_xy) - centroid) / sqrt(kg)
}

evaluate_pooled_fg_density <- function(query_xy, centroid, pooled, apply_geo = FALSE, kappa_geo = 1, fg_draw = NULL, lambda_cloud = 0.7, density_eps = 1e-300, param_eps = 1e-8, return_draw = FALSE) {
  if (is.null(fg_draw)) fg_draw <- sample_fgkmm_domain_from_pooled(pooled = pooled, lambda_cloud = lambda_cloud, eps = param_eps)
  Xq <- as.matrix(query_xy); centroid <- as.numeric(centroid)
  if (apply_geo) {
    kg <- as.numeric(kappa_geo); if (!is.finite(kg) || kg <= 0) kg <- 1
    X_eval <- .inverse_warp_query(query_xy = Xq, centroid = centroid, kappa_geo = kg)
    jac <- 1 / kg
  } else { X_eval <- Xq; jac <- 1 }
  dens0 <- evaluate_fg_mixture_density(query_xy = X_eval, centroid = centroid, fg_draw = fg_draw, eps = density_eps)
  dens <- jac * dens0; dens[!is.finite(dens)] <- density_eps; dens <- pmax(dens, density_eps)
  if (return_draw) return(list(density = dens, fg_draw = fg_draw, query_eval = X_eval, jacobian = jac))
  dens
}

.normalize_coords_unit <- function(coords_df) {
  out <- coords_df
  xr <- max(out$x, na.rm = TRUE) - min(out$x, na.rm = TRUE)
  yr <- max(out$y, na.rm = TRUE) - min(out$y, na.rm = TRUE)
  if (!is.finite(xr) || xr <= 0) xr <- 1
  if (!is.finite(yr) || yr <= 0) yr <- 1
  out$x <- (out$x - min(out$x, na.rm = TRUE)) / xr
  out$y <- (out$y - min(out$y, na.rm = TRUE)) / yr
  out
}

.pick_fixed_grid_template <- function(pilot_data, grp_id = NULL) {
  if (length(pilot_data) == 0L) stop("pilot_data is empty.")
  idx_pool <- seq_along(pilot_data)
  if (!is.null(grp_id)) {
    idx_grp <- which(vapply(pilot_data, function(x) as.character(x$group) == as.character(grp_id), logical(1)))
    if (length(idx_grp) > 0L) idx_pool <- idx_grp
  }
  pick <- sample(idx_pool, size = 1L)
  coords <- .normalize_coords_unit(pilot_data[[pick]]$coords)
  list(template_index = pick, grid_df = data.frame(x = coords$x, y = coords$y, stringsAsFactors = FALSE))
}

.extract_geometry_group <- function(params_geometry, grp_chr) {
  if (is.null(params_geometry) || length(params_geometry) == 0L) return(NULL)
  if (!is.null(params_geometry[[grp_chr]])) return(params_geometry[[grp_chr]])
  alt <- if (grp_chr == "0") "1" else "0"
  if (!is.null(params_geometry[[alt]])) {
    warning(sprintf(".extract_geometry_group: key '%s' not found; falling back to '%s'.", grp_chr, alt))
    return(params_geometry[[alt]])
  }
  params_geometry[[1]]
}

build_spatial_cache <- function(sim_coords_df, spatial_mode = c("original", "basis"), graph_k = 20L, basis_rank = 25L, basis_knot_method = c("sample", "grid"), basis_seed = 1L) {
  spatial_mode <- match.arg(spatial_mode); basis_knot_method <- match.arg(basis_knot_method)
  if (!all(c("x", "y", "domain") %in% names(sim_coords_df))) stop("sim_coords_df must contain x, y, domain.")
  coords_mat <- as.matrix(sim_coords_df[, c("x", "y"), drop = FALSE]); n_cells <- nrow(coords_mat); sim_domains <- unique(as.character(sim_coords_df$domain))
  out <- list(spatial_mode = spatial_mode, coords_mat = coords_mat, n_cells = n_cells, sim_domains = sim_domains)
  if (spatial_mode == "original") {
    if (!requireNamespace("FNN", quietly = TRUE)) stop("Need FNN for spatial_mode='original'.")
    k_use <- min(as.integer(graph_k), max(1L, n_cells - 1L))
    if (n_cells <= 1L) {
      out$graph_k <- 0L; out$nn_index <- matrix(integer(0), nrow = n_cells, ncol = 0L); out$nn_dist  <- matrix(numeric(0), nrow = n_cells, ncol = 0L)
    } else {
      nn <- FNN::get.knn(coords_mat, k = k_use)
      out$graph_k  <- k_use; out$nn_index <- nn$nn.index; out$nn_dist  <- nn$nn.dist
    }
  }
  if (spatial_mode == "basis") {
    basis_rank <- min(as.integer(basis_rank), n_cells); basis_rank <- max(1L, basis_rank)
    if (basis_knot_method == "sample") {
      ord <- order(coords_mat[, 1], coords_mat[, 2]); pick <- unique(round(seq(1, n_cells, length.out = basis_rank))); knots <- coords_mat[ord[pick], , drop = FALSE]
    } else {
      gx <- ceiling(sqrt(basis_rank)); gy <- ceiling(basis_rank / gx)
      x_seq <- seq(min(coords_mat[, 1]), max(coords_mat[, 1]), length.out = gx); y_seq <- seq(min(coords_mat[, 2]), max(coords_mat[, 2]), length.out = gy)
      grid_df <- expand.grid(x = x_seq, y = y_seq)
      knots <- as.matrix(grid_df[seq_len(min(nrow(grid_df), basis_rank)), c("x", "y"), drop = FALSE])
    }
    out$basis_rank <- nrow(knots); out$basis_knots <- knots; out$basis_knot_method <- basis_knot_method; out$basis_seed <- basis_seed
  }
  out
}

#' extract_theta_triplet <- function(theta) {
#'   if (is.null(theta)) {
#'     return(c(alpha = 0.2, rho = 10, nugget = 0.05))
#'   }
#'
#'   nm <- names(theta)
#'
#'   alpha_sp <- NA_real_
#'   rho_sp   <- NA_real_
#'   nugget2  <- NA_real_
#'
#'   # case 1: spaCraft legacy storage from estimateExpressionParams()
#'   # theta = c(alpha=..., rho=..., sigma.sq=...) where sigma.sq means nugget
#'   if ("alpha" %in% nm) alpha_sp <- as.numeric(theta[["alpha"]])
#'   if ("rho"   %in% nm) rho_sp   <- as.numeric(theta[["rho"]])
#'
#'   if ("nugget" %in% nm) {
#'     nugget2 <- as.numeric(theta[["nugget"]])
#'   } else if ("alpha" %in% nm && "sigma.sq" %in% nm) {
#'     nugget2 <- as.numeric(theta[["sigma.sq"]])
#'   }
#'
#'   # case 2: raw BRISC Theta passed directly by mistake
#'   if (!is.finite(alpha_sp) && "sigma.sq" %in% nm) alpha_sp <- as.numeric(theta[["sigma.sq"]])
#'   if (!is.finite(rho_sp)   && "phi"      %in% nm) rho_sp   <- as.numeric(theta[["phi"]])
#'   if (!is.finite(nugget2)  && "tau.sq"   %in% nm) nugget2  <- as.numeric(theta[["tau.sq"]])
#'
#'   # sanitize
#'   if (!is.finite(alpha_sp) || alpha_sp < 0) alpha_sp <- 0
#'   if (!is.finite(rho_sp)   || rho_sp <= 0)  rho_sp   <- 10
#'   if (!is.finite(nugget2)  || nugget2 < 0)  nugget2  <- 0
#'
#'   c(alpha = alpha_sp, rho = rho_sp, nugget = nugget2)
#' }
#'
#'
#'
#'
#' # ------------------------------------------------------------------------------
#' # Low-rank pivoted Cholesky GP draw (fallback / legacy helper)
#' # ------------------------------------------------------------------------------
#' low_rank_gp_draw_from_cov <- function(K, rank_max = 20L, jitter = 1e-8) {
#'   n <- nrow(K)
#'   if (n == 0L) return(numeric(0))
#'   if (n == 1L) {
#'     v <- max(K[1, 1], 0)
#'     return(stats::rnorm(1, mean = 0, sd = sqrt(v)))
#'   }
#'
#'   diag(K) <- diag(K) + jitter
#'
#'   out <- tryCatch({
#'     R <- chol(K, pivot = TRUE)
#'
#'     piv <- attr(R, "pivot")
#'     rk  <- attr(R, "rank")
#'     if (is.null(rk) || !is.finite(rk)) rk <- n
#'
#'     r_use <- min(as.integer(rank_max), as.integer(rk), n)
#'     r_use <- max(r_use, 1L)
#'
#'     R_r <- R[seq_len(r_use), , drop = FALSE]
#'
#'     z_r <- stats::rnorm(r_use)
#'     x_piv <- as.numeric(crossprod(R_r, z_r))
#'
#'     x <- numeric(n)
#'     x[piv] <- x_piv
#'     x
#'   }, error = function(e) {
#'     tryCatch({
#'       U <- chol(K, pivot = FALSE)
#'       as.numeric(t(U) %*% stats::rnorm(n))
#'     }, error = function(e2) {
#'       stats::rnorm(n, mean = 0, sd = sqrt(max(mean(diag(K)), 0)))
#'     })
#'   })
#'
#'   out[!is.finite(out)] <- 0
#'   out
#' }
#'
#'
#' # ------------------------------------------------------------------------------
#' # Utility: center + rescale to target variance alpha
#' # ------------------------------------------------------------------------------
#' rescale_spatial_field <- function(x, alpha) {
#'   x <- as.numeric(x)
#'   x[!is.finite(x)] <- 0
#'   x <- x - mean(x)
#'
#'   sx <- stats::sd(x)
#'   if (!is.finite(sx) || sx < 1e-10 || !is.finite(alpha) || alpha <= 0) {
#'     return(rep(0, length(x)))
#'   }
#'
#'   sqrt(alpha) * (x / sx)
#' }
#'
#'
#' # ------------------------------------------------------------------------------
#' # Whole-tissue original graph-based spatial draw
#' #   - uses whole-tissue KNN graph
#' #   - rho controls graph weights
#' #   - alpha controls marginal scale after rescaling
#' # ------------------------------------------------------------------------------
#' draw_graph_gp_whole_tissue <- function(spatial_cache, alpha, rho, include_self = TRUE) {
#'   n <- spatial_cache$n_cells
#'   if (n <= 0L || alpha <= 0) return(numeric(n))
#'   if (!is.finite(rho) || rho <= 0) rho <- 10
#'
#'   nn_index <- spatial_cache$nn_index
#'   nn_dist  <- spatial_cache$nn_dist
#'   k_use <- ncol(nn_index)
#'
#'   z <- stats::rnorm(n)
#'
#'   if (k_use == 0L) {
#'     return(rescale_spatial_field(z, alpha))
#'   }
#'
#'   W <- exp(-nn_dist / rho)
#'
#'   if (include_self) {
#'     z_nb <- cbind(z, matrix(z[nn_index], nrow = n))
#'     W_all <- cbind(rep(1, n), W)
#'   } else {
#'     z_nb <- matrix(z[nn_index], nrow = n)
#'     W_all <- W
#'   }
#'
#'   rs <- rowSums(W_all)
#'   rs[!is.finite(rs) | rs <= 0] <- 1
#'   W_all <- W_all / rs
#'
#'   eta <- rowSums(W_all * z_nb)
#'   rescale_spatial_field(eta, alpha)
#' }
#'
#'
#' # ------------------------------------------------------------------------------
#' # Basis matrix builder for whole-tissue low-rank basis mode
#' # ------------------------------------------------------------------------------
#' build_basis_matrix <- function(coords_mat, knots, rho,
#'                                basis_kernel = c("gaussian", "exponential"),
#'                                normalize_columns = TRUE) {
#'   basis_kernel <- match.arg(basis_kernel)
#'
#'   n <- nrow(coords_mat)
#'   m <- nrow(knots)
#'   if (n == 0L || m == 0L) return(matrix(0, n, m))
#'
#'   if (!is.finite(rho) || rho <= 0) rho <- 10
#'
#'   dx <- outer(coords_mat[, 1], knots[, 1], "-")
#'   dy <- outer(coords_mat[, 2], knots[, 2], "-")
#'   d  <- sqrt(dx^2 + dy^2)
#'
#'   Phi <- switch(
#'     basis_kernel,
#'     gaussian    = exp(-0.5 * (d / rho)^2),
#'     exponential = exp(-d / rho)
#'   )
#'
#'   if (normalize_columns && ncol(Phi) > 0L) {
#'     cn <- sqrt(colSums(Phi^2))
#'     cn[!is.finite(cn) | cn <= 0] <- 1
#'     Phi <- sweep(Phi, 2, cn, "/", check.margin = FALSE)
#'   }
#'
#'   Phi
#' }
#'
#' # -------------------------------------------------------------------------
#' # Whole-tissue low-rank basis GP draw
#' # -------------------------------------------------------------------------
#' draw_basis_gp_whole_tissue <- function(spatial_cache,
#'                                        alpha,
#'                                        rho,
#'                                        basis_kernel = c("gaussian", "exponential")) {
#'   basis_kernel <- match.arg(basis_kernel)
#'
#'   n <- spatial_cache$n_cells
#'   if (n <= 0L || !is.finite(alpha) || alpha <= 0) return(numeric(n))
#'
#'   coords_mat <- spatial_cache$coords_mat
#'   knots <- spatial_cache$basis_knots
#'
#'   if (is.null(knots) || nrow(knots) == 0L) {
#'     return(rescale_spatial_field(stats::rnorm(n), alpha))
#'   }
#'
#'   Phi <- build_basis_matrix(
#'     coords_mat = coords_mat,
#'     knots = knots,
#'     rho = rho,
#'     basis_kernel = basis_kernel,
#'     normalize_columns = TRUE
#'   )
#'
#'   m <- ncol(Phi)
#'   if (m == 0L) {
#'     return(rescale_spatial_field(stats::rnorm(n), alpha))
#'   }
#'
#'   z <- stats::rnorm(m)
#'   eta <- as.numeric(Phi %*% z)
#'   eta[!is.finite(eta)] <- 0
#'
#'   rescale_spatial_field(eta, alpha)
#' }
#'
#'
#'
#' # current work well for LFC, BUT NOT DOAMIN COMPOSITION#
#' gene_worker <- function(g,
#'                            sim_coords_df,
#'                            sim_coords_mat,
#'                            sim_domains,
#'                            pilot_bundle    = NULL,
#'                            p_info_override = NULL,
#'                            params_expr_grp = NULL,
#'                            is_case,
#'                            target_domain,
#'                            de_genes,
#'                            scenario,
#'                            k_nn          = 20L,
#'                            spatial_cache = NULL,
#'                            spatial_mode  = c("original", "basis"),
#'                            gp_rank       = 40L,
#'                            gp_jitter     = 1e-8,
#'                            basis_kernel  = c("gaussian", "exponential"),
#'                            # legacy compatibility
#'                            pilot_sample     = NULL,
#'                            pilot_coords_mat = NULL) {
#'
#'   if (!requireNamespace("FNN", quietly = TRUE)) stop("Package 'FNN' required.")
#'
#'   spatial_mode <- match.arg(spatial_mode)
#'   basis_kernel <- match.arg(basis_kernel)
#'
#'   sim_coords_mat <- as.matrix(sim_coords_mat)
#'   n_cells        <- nrow(sim_coords_mat)
#'   if (n_cells == 0L) return(numeric(0))
#'
#'   # ------------------------------------------------------------------
#'   # 0) Synthetic domain vector (generator-level latent labels)
#'   # ------------------------------------------------------------------
#'   domain_vec <- if ("tilde_d" %in% names(sim_coords_df)) {
#'     as.character(sim_coords_df$tilde_d)
#'   } else {
#'     as.character(sim_coords_df$domain)
#'   }
#'
#'   # ------------------------------------------------------------------
#'   # 1) Parameter extraction (unchanged from v1)
#'   # ------------------------------------------------------------------
#'   if (!is.null(p_info_override)) {
#'     mu_grand  <- p_info_override$mu_grand
#'     theta     <- p_info_override$theta
#'     sigma_bio <- p_info_override$sigma_bio
#'   } else {
#'     p_ex <- if (!is.null(params_expr_grp)) params_expr_grp[[g]] else NULL
#'     if (is.null(p_ex)) {
#'       mu_grand  <- stats::setNames(rep(0.1, length(sim_domains)), sim_domains)
#'       theta     <- c(alpha = 0.2, rho = 10, sigma.sq = 0.05)
#'       sigma_bio <- 0
#'     } else {
#'       mu_grand  <- p_ex$mu_grand
#'       theta     <- p_ex$theta
#'       sigma_bio <- if (!is.null(p_ex$sigma_bio)) p_ex$sigma_bio else 0
#'     }
#'   }
#'
#'   mu_grand <- as.numeric(mu_grand)
#'   if (!is.null(p_info_override) && !is.null(names(p_info_override$mu_grand))) {
#'     names(mu_grand) <- names(p_info_override$mu_grand)
#'   } else if (!is.null(params_expr_grp) &&
#'              !is.null(params_expr_grp[[g]]) &&
#'              !is.null(names(params_expr_grp[[g]]$mu_grand))) {
#'     names(mu_grand) <- names(params_expr_grp[[g]]$mu_grand)
#'   } else {
#'     names(mu_grand) <- sim_domains
#'   }
#'
#'   th       <- extract_theta_triplet(theta)
#'   alpha_sp <- as.numeric(th["alpha"])
#'   rho_sp   <- as.numeric(th["rho"])
#'   nugget2  <- as.numeric(th["nugget"])
#'
#'   if (!is.finite(alpha_sp) || alpha_sp < 0) alpha_sp <- 0
#'   if (!is.finite(rho_sp)   || rho_sp <= 0)  rho_sp   <- 10
#'   if (!is.finite(nugget2)  || nugget2 < 0)  nugget2  <- 0
#'
#'   # Case-specific range perturbation (delta_rho)
#'   rho_eff <- rho_sp
#'   if (is_case) {
#'     delta_rho <- as.numeric(scenario$delta_rho)
#'     if (is.finite(delta_rho) && delta_rho > 0) rho_eff <- rho_sp * delta_rho
#'   }
#'
#'   # ------------------------------------------------------------------
#'   # 2) Pilot sample selection (unchanged from v1)
#'   # ------------------------------------------------------------------
#'   current_pilot_data   <- NULL
#'   current_pilot_coords <- NULL
#'
#'   if (!is.null(pilot_bundle) &&
#'       !is.null(pilot_bundle$pilot_samples) &&
#'       length(pilot_bundle$pilot_samples) > 0L) {
#'
#'     w <- pilot_bundle$pilot_mix_w
#'     if (is.null(w) || length(w) != length(pilot_bundle$pilot_samples) ||
#'         any(!is.finite(w)) || sum(w) <= 0) {
#'       w <- rep(1, length(pilot_bundle$pilot_samples))
#'     }
#'     w   <- w / sum(w)
#'     sel <- sample(seq_along(pilot_bundle$pilot_samples), size = 1L, prob = w)
#'
#'     current_pilot_data <- pilot_bundle$pilot_samples[[sel]]
#'     if (!is.null(pilot_bundle$pilot_coords_mat_list) &&
#'         length(pilot_bundle$pilot_coords_mat_list) >= sel) {
#'       current_pilot_coords <- as.matrix(pilot_bundle$pilot_coords_mat_list[[sel]])
#'     }
#'   } else {
#'     current_pilot_data   <- pilot_sample
#'     current_pilot_coords <- if (!is.null(pilot_coords_mat)) as.matrix(pilot_coords_mat) else NULL
#'   }
#'
#'   # ------------------------------------------------------------------
#'   # 3) Mean surface + DE shift (unchanged from v1)
#'   # ------------------------------------------------------------------
#'   mu_log      <- numeric(n_cells)
#'   mu_fallback <- mean(mu_grand[is.finite(mu_grand)], na.rm = TRUE)
#'   if (!is.finite(mu_fallback)) mu_fallback <- 0.1
#'
#'   de_shift <- as.numeric(scenario$DE_lfc)
#'   if (!is.finite(de_shift)) de_shift <- 0
#'
#'   for (dom in unique(domain_vec)) {
#'     base <- if (!is.null(names(mu_grand)) && dom %in% names(mu_grand)) {
#'       mu_grand[[dom]]
#'     } else {
#'       mu_fallback
#'     }
#'     if (!is.finite(base)) base <- 0.1
#'
#'     if (is_case &&
#'         identical(as.character(dom), as.character(target_domain)) &&
#'         (g %in% de_genes)) {
#'       base <- base + de_shift
#'     }
#'     mu_log[domain_vec == dom] <- base
#'   }
#'
#'   # ------------------------------------------------------------------
#'   # 4) Between-sample biological shift (unchanged from v1)
#'   # ------------------------------------------------------------------
#'   kappa_bio     <- as.numeric(scenario$kappa_bio %||% 1)
#'   if (!is.finite(kappa_bio) || kappa_bio < 0) kappa_bio <- 1
#'   sigma_bio_new <- as.numeric(sigma_bio) * kappa_bio
#'   if (!is.finite(sigma_bio_new) || sigma_bio_new < 0) sigma_bio_new <- 0
#'   B_shift <- stats::rnorm(1, mean = 0, sd = sqrt(sigma_bio_new))
#'
#'   # ------------------------------------------------------------------
#'   # 5) Spatial field: eta = lambda * eta_cond + (1-lambda) * eta_para
#'   # ------------------------------------------------------------------
#'   lambda <- as.numeric(scenario$lambda_cond %||% 0)
#'   if (!is.finite(lambda)) lambda <- 0
#'   lambda <- max(0, min(1, lambda))
#'
#'   # =================================================================
#'   # 5A) eta_cond: DOMAIN-AWARE pilot-guided conditional texture
#'   #
#'   # KEY CHANGE FROM v1:
#'   #   v1: KNN from synthetic spot to ALL pilot spots (domain-blind)
#'   #       -> synthetic HIPPO spot can get CTX texture -> classification fail
#'   #
#'   #   v2: KNN from synthetic spot to pilot spots OF THE SAME DOMAIN
#'   #       -> synthetic HIPPO spot gets HIPPO texture -> classification success
#'   #
#'   # For each unique synthetic domain d:
#'   #   1. Select pilot spots with pilot_domain == d
#'   #   2. Compute pilot residuals for those spots
#'   #   3. KNN: synthetic spots of domain d -> pilot spots of domain d
#'   #   4. Kernel-smooth the domain-matched residuals
#'   # =================================================================
#'   eta_cond <- rep(0, n_cells)
#'
#'   if (lambda > 1e-8 &&
#'       !is.null(current_pilot_data) &&
#'       !is.null(current_pilot_coords) &&
#'       nrow(current_pilot_coords) > 0L) {
#'
#'     gene_in_logcounts <- !is.null(current_pilot_data$logcounts) &&
#'       (g %in% rownames(current_pilot_data$logcounts))
#'     gene_in_counts    <- !is.null(current_pilot_data$counts) &&
#'       (g %in% rownames(current_pilot_data$counts))
#'
#'     if (gene_in_logcounts || gene_in_counts) {
#'
#'       # Full pilot expression vector
#'       p_log <- if (gene_in_logcounts) {
#'         as.numeric(current_pilot_data$logcounts[g, ])
#'       } else {
#'         log1p(as.numeric(current_pilot_data$counts[g, ]))
#'       }
#'
#'       # Pilot domain labels
#'       pilot_dom <- if (!is.null(current_pilot_data$coords) &&
#'                        "domain" %in% names(current_pilot_data$coords)) {
#'         as.character(current_pilot_data$coords$domain)
#'       } else {
#'         rep(NA_character_, length(p_log))
#'       }
#'
#'       # Pilot residuals: r_jg = log1p(X_jg) - mu_{d^pilot_j, c}
#'       p_mu <- vapply(pilot_dom, function(d) {
#'         vv <- if (!is.null(names(mu_grand)) && !is.na(d) && d %in% names(mu_grand)) {
#'           mu_grand[[d]]
#'         } else {
#'           mu_fallback
#'         }
#'         if (!is.finite(vv)) vv <- 0.1
#'         vv
#'       }, numeric(1))
#'
#'       p_res <- p_log - p_mu
#'       p_res[!is.finite(p_res)] <- 0
#'       p_res <- p_res - mean(p_res)   # centre globally
#'
#'       # --- Domain-aware KNN: per-domain processing ---
#'       for (dom in unique(domain_vec)) {
#'         syn_idx   <- which(domain_vec == dom)
#'         if (length(syn_idx) == 0L) next
#'
#'         # Pilot spots of the SAME domain
#'         pilot_idx <- which(pilot_dom == dom)
#'
#'         # Fallback: if no pilot spots of this domain, use all pilot spots
#'         if (length(pilot_idx) < 3L) {
#'           pilot_idx <- seq_along(pilot_dom)
#'         }
#'
#'         k_use <- min(as.integer(k_nn), length(pilot_idx))
#'         k_use <- max(k_use, 1L)
#'
#'         nn <- FNN::get.knnx(
#'           data  = current_pilot_coords[pilot_idx, , drop = FALSE],
#'           query = sim_coords_mat[syn_idx, , drop = FALSE],
#'           k     = k_use
#'         )
#'
#'         d2 <- nn$nn.dist^2
#'         h  <- stats::median(d2[is.finite(d2)], na.rm = TRUE)
#'         if (!is.finite(h) || h <= 0) h <- mean(d2[is.finite(d2)], na.rm = TRUE)
#'         if (!is.finite(h) || h <= 0) h <- 1
#'
#'         W  <- exp(-d2 / (h + 1e-8))
#'         rs <- rowSums(W)
#'         rs[!is.finite(rs) | rs <= 0] <- 1
#'         W  <- W / rs
#'
#'         # Pilot residuals restricted to domain-matched spots
#'         p_res_dom <- p_res[pilot_idx]
#'         eta_cond[syn_idx] <- rowSums(
#'           W * matrix(p_res_dom[nn$nn.index], nrow = length(syn_idx))
#'         )
#'       }
#'
#'       eta_cond[!is.finite(eta_cond)] <- 0
#'       eta_cond <- eta_cond - mean(eta_cond)   # re-centre after domain-wise fill
#'     }
#'   }
#'
#'   # -- 5B) eta_para: parametric whole-tissue spatial draw (unchanged) ---
#'   eta_para <- rep(0, n_cells)
#'
#'   if (alpha_sp > 0) {
#'     used <- FALSE
#'
#'     if (!is.null(spatial_cache)) {
#'       if (spatial_mode == "original" &&
#'           exists("draw_graph_gp_whole_tissue", mode = "function")) {
#'         eta_para <- draw_graph_gp_whole_tissue(spatial_cache, alpha_sp, rho_eff,
#'                                                include_self = TRUE)
#'         used <- TRUE
#'       }
#'       if (spatial_mode == "basis" &&
#'           exists("draw_basis_gp_whole_tissue", mode = "function")) {
#'         eta_para <- draw_basis_gp_whole_tissue(spatial_cache, alpha_sp, rho_eff,
#'                                                basis_kernel)
#'         used <- TRUE
#'       }
#'     }
#'
#'     if (!used) {
#'       if (n_cells <= 1L) {
#'         eta_para <- stats::rnorm(n_cells, 0, sqrt(alpha_sp))
#'       } else {
#'         D_all    <- as.matrix(stats::dist(sim_coords_mat))
#'         K_all    <- alpha_sp * exp(-D_all / rho_eff)
#'         diag(K_all) <- diag(K_all) + gp_jitter
#'         eta_para <- low_rank_gp_draw_from_cov(
#'           K_all,
#'           rank_max = min(as.integer(gp_rank), n_cells),
#'           jitter   = gp_jitter
#'         )
#'       }
#'       eta_para <- rescale_spatial_field(eta_para, alpha_sp)
#'     }
#'   }
#'   eta_para[!is.finite(eta_para)] <- 0
#'
#'   # -- 5C) Blended spatial residual ---------------------------------
#'   eta <- lambda * eta_cond + (1 - lambda) * eta_para
#'   eta[!is.finite(eta)] <- 0
#'
#'   # ------------------------------------------------------------------
#'   # 6) Nugget noise
#'   # ------------------------------------------------------------------
#'   eps_noise <- if (nugget2 > 0) {
#'     stats::rnorm(n_cells, 0, sqrt(nugget2))
#'   } else {
#'     rep(0, n_cells)
#'   }
#'
#'   # ------------------------------------------------------------------
#'   # 7) Final log1p(count)
#'   # ------------------------------------------------------------------
#'   log1p_y <- mu_log + B_shift + eta + eps_noise
#'   log1p_y[!is.finite(log1p_y)] <- 0
#'   log1p_y[log1p_y < 0]         <- 0
#'
#'   return(log1p_y)
#' }
#'
#'
#'
#'
#'
#' get_conditional_parameters <- function(
#'     mu_full, cov_full, observed_indices, observed_values) {
#'
#'   # Partition mean vector
#'   mu_obs <- mu_full[observed_indices]
#'   mu_unobs <- mu_full[-observed_indices]
#'
#'   # Partition covariance matrix
#'   Sigma_obs_obs <- cov_full[observed_indices, observed_indices]
#'   Sigma_unobs_unobs <- cov_full[-observed_indices, -observed_indices]
#'   Sigma_unobs_obs <- cov_full[-observed_indices, observed_indices] # Sigma_unobs_obs is Sigma_YX
#'
#'   # Inverse of observed block (add jitter for numerical stability if needed)
#'   diag(Sigma_obs_obs) <- diag(Sigma_obs_obs) + 1e-8 # Add jitter
#'
#'   # Use tryCatch for robust inversion
#'   Sigma_obs_obs_inv <- tryCatch(solve(Sigma_obs_obs), error = function(e) {
#'     warning(paste("Singular matrix encountered for Sigma_obs_obs:", e$message, "Returning NULL."))
#'     return(NULL)
#'   })
#'
#'   if (is.null(Sigma_obs_obs_inv)) {
#'     return(list(mu_conditional = NULL, Sigma_conditional_chol = NULL)) # Return NULL for chol too
#'   }
#'
#'   # Pre-calculate common matrix products
#'   Sigma_uo_Sigma_oo_inv <- Sigma_unobs_obs %*% Sigma_obs_obs_inv
#'
#'   # Conditional Mean
#'   mu_conditional <- mu_unobs + Sigma_uo_Sigma_oo_inv %*% (observed_values - mu_obs)
#'
#'   # Conditional Covariance
#'   Sigma_conditional <- Sigma_unobs_unobs - Sigma_uo_Sigma_oo_inv %*% t(Sigma_unobs_obs)
#'
#'   # Add jitter to conditional covariance for Cholesky decomposition
#'   diag(Sigma_conditional) <- diag(Sigma_conditional) + 1e-8
#'
#'   # Pre-calculate Cholesky of conditional covariance
#'   Sigma_conditional_chol <- tryCatch(chol(Sigma_conditional), error = function(e) {
#'     warning(paste("Non-positive definite conditional covariance:", e$message, "Returning NULL."))
#'     return(NULL)
#'   })
#'
#'   return(list(mu_conditional = mu_conditional, Sigma_conditional_chol = Sigma_conditional_chol))
#' }
#'
#'
#'
#' generate_fitted_data <- function(
#'     mu_est, alpha_spatial_est, rho_est, sigma2_B_est, sigma2_noise_est,
#'     coords, K_syn, Y_list_original, subratio=0.5, seed = NULL) { # subratio default 0.7 as per directions
#'
#'   if (!is.null(seed)) set.seed(seed)
#'
#'   # --- Initial Setup ---
#'   n <- nrow(coords)
#'   p <- length(mu_est)
#'   D <- as.matrix(stats::dist(coords))
#'
#'   # Ensure parameters are vectors of the correct length
#'   if (length(mu_est) == 1) mu_est <- rep(mu_est, p)
#'   if (length(rho_est) == 1) rho_est <- rep(rho_est, p)
#'   if (length(alpha_spatial_est) == 1) alpha_spatial_est <- rep(alpha_spatial_est, p)
#'   if (length(sigma2_B_est) == 1) sigma2_B_est <- rep(sigma2_B_est, p)
#'   if (length(sigma2_noise_est) == 1) sigma2_noise_est <- rep(sigma2_noise_est, p)
#'
#'   # MODIFIED: Moved outside loop
#'   n_obs <- round(subratio * n)
#'   observed_indices <- sample(1:n, n_obs)
#'   unobserved_indices <- setdiff(1:n, observed_indices)
#'
#'   Y_original_array <- simplify2array(Y_list_original) # n x p x K
#'   Y_bar_original <- apply(Y_original_array, c(1, 2), mean, na.rm = TRUE) # n x p
#'   Y_bar_sub <- Y_bar_original[observed_indices, ] # n_obs x p
#'
#'   # MODIFIED START: Pre-calculate conditional parameters for ALL genes (g) outside K_syn loop
#'   precomputed_cond_params <- lapply(1:p, function(g) {
#'     # Full spatial covariance for gene g
#'     Sigma_spatial_full_g <- alpha_spatial_est[g] * exp(-D / rho_est[g])
#'
#'     # Total covariance for Y_g vector for current gene g
#'     inflation_factor <- 1
#'     Sigma_Y_full_g <- Sigma_spatial_full_g + diag(sigma2_B_est[g] + inflation_factor*sigma2_noise_est[g], n)
#'
#'     # Observed Y_bar for this gene
#'     observed_y_bar_g <- Y_bar_sub[, g]
#'
#'     # Mean vector for the full Y_g
#'     mu_full_g <- rep(mu_est[g], n)
#'
#'     # Get conditional parameters (now includes pre-computed Cholesky)
#'     get_conditional_parameters(
#'       mu_full = mu_full_g,
#'       cov_full = Sigma_Y_full_g,
#'       observed_indices = observed_indices,
#'       observed_values = observed_y_bar_g
#'     )
#'   })
#'   # MODIFIED END
#'
#'   Y_fitted_list <- lapply(1:K_syn, function(k_syn) {
#'     Y_k_syn <- matrix(NA, nrow = n, ncol = p)
#'
#'     # MODIFIED START: Vectorize random number generation for all genes
#'     # Generate all standard normal random vectors for unobserved parts at once
#'     set.seed(seed+k_syn)
#'     Z_unobs_all_genes <- matrix(rnorm(length(unobserved_indices) * p),
#'                                 nrow = length(unobserved_indices), ncol = p)
#'
#'     for (g in 1:p) {
#'       current_cond_params <- precomputed_cond_params[[g]]
#'
#'       if (is.null(current_cond_params$mu_conditional) || is.null(current_cond_params$Sigma_conditional_chol)) {
#'         # If conditional params could not be computed for this gene, fill with NA
#'         Y_k_syn[, g] <- NA
#'         warning(paste("Skipping gene", g, "due to failed conditional parameter computation."))
#'       } else {
#'         # Sample from the conditional MVN using precomputed Cholesky and vectorized Z
#'         unobserved_values_g_sampled <- current_cond_params$mu_conditional +
#'           t(current_cond_params$Sigma_conditional_chol) %*% Z_unobs_all_genes[, g]
#'
#'         temp_Y_g <- numeric(n)
#'         temp_Y_g[observed_indices] <- Y_bar_sub[, g] # Use the observed Y_bar_sub for fixed values
#'         temp_Y_g[unobserved_indices] <- unobserved_values_g_sampled
#'         Y_k_syn[, g] <- temp_Y_g
#'       }
#'     }
#'     # MODIFIED END
#'
#'     return(Y_k_syn)
#'   })
#'
#'   return(Y_fitted_list)
#' }
#'
#'
#'
#'
#' # MODIFIED: estimate_parameters_by_group (wrapper)
#' # Added new arguments: effective_gene_indices, non_effective_gene_indices
#' estimate_parameters_by_group <- function(
#'     Y_list_by_group, groups, n, p, coords, p_B, sigma_B_sq_vec=NULL,
#'     effective_gene_indices, non_effective_gene_indices # NEW ARGUMENTS
#' ) {
#'   # Determine K (number of pilot samples) for a representative group
#'   if (length(groups) == 0 || length(Y_list_by_group[[as.character(groups[1])]]) == 0) {
#'     stop("Y_list_by_group is empty or malformed. Cannot determine K.")
#'   }
#'
#'   K <- length(Y_list_by_group[[as.character(groups[1])]])
#'
#'   if (K > 1) {
#'     # MODIFIED: Pass new arguments to multiple_K version
#'     result <- estimate_parameters_by_group_multiple_K(
#'       Y_list_by_group, groups, n, p, coords, p_B, sigma_B_sq_vec,
#'       effective_gene_indices, non_effective_gene_indices
#'     )
#'   } else if (K == 1) {
#'     # MODIFIED: Pass new arguments to K1 version
#'     result <- estimate_parameters_by_group_K1(
#'       Y_list_by_group, groups, n, p, coords, p_B,
#'       effective_gene_indices, non_effective_gene_indices
#'     )
#'   } else {
#'     stop("Invalid value for K (number of pilot samples). K must be >= 1.")
#'   }
#'
#'   # Attach the original Y_list_by_group to the result.
#'   result$Y_list_original_for_conditional_sampling <- Y_list_by_group
#'
#'   return(result)
#' }
#'
#' estimate_parameters_by_group_multiple_K <- function(
#'     Y_list_by_group, groups, n, p, coords, p_B, sigma_B_sq_vec,
#'     effective_gene_indices, non_effective_gene_indices
#' ) {
#'   results_by_group <- list()
#'
#'   for (c in seq_along(groups)) {
#'     group <- groups[c]
#'     Y_list      <- Y_list_by_group[[group]]
#'     p_B_group   <- p_B[c]
#'     sigma_B_sq  <- sigma_B_sq_vec[group] # True value for comparison
#'     K <- length(Y_list)
#'
#'     # --- Step 1: Fit BRISC on ORIGINAL data for each sample and gene ---
#'     brisc_param_array <- array(NA, dim = c(p, 4, K),
#'                                dimnames = list(paste0("gene", 1:p),
#'                                                c("mu", "alpha_spatial", "rho", "sigma2_noise"),
#'                                                paste0("sample", 1:K)))
#'     n_neighbors = 10
#'     task_grid <- expand.grid(g = 1:p, k = 1:K)
#'
#'     results_list <- pbmclapply(1:nrow(task_grid), function(i) {
#'       g <- task_grid$g[i]
#'       k <- task_grid$k[i]
#'       Y_gk <- Y_list[[k]][, g]
#'
#'       if (any(!is.finite(Y_gk))) return(rep(NA, 4))
#'       if (var(Y_gk, na.rm = TRUE) == 0) Y_gk <- Y_gk + rnorm(length(Y_gk), 0, 1e-6)
#'
#'       brisc_fit <- suppressMessages(try(BRISC::BRISC_estimation(
#'         coords = coords, y = Y_gk, x = matrix(1, nrow = n),
#'         cov.model = "exponential", n.neighbors = n_neighbors, order = "AMMD",
#'         verbose = FALSE, nugget_status = 1
#'       ), silent = TRUE))
#'
#'       if (!inherits(brisc_fit, "try-error")) {
#'         return(c(mu = as.numeric(brisc_fit$Beta),
#'                  alpha_spatial = brisc_fit$Theta["sigma.sq"],
#'                  rho = brisc_fit$Theta["phi"],
#'                  sigma2_noise = brisc_fit$Theta["tau.sq"]))
#'       } else {
#'         return(rep(NA, 4))
#'       }
#'     }, mc.cores = N_CORES)
#'
#'     for (i in 1:nrow(task_grid)) {
#'       g <- task_grid$g[i]
#'       k <- task_grid$k[i]
#'       brisc_param_array[g, , k] <- results_list[[i]]
#'     }
#'
#'     # --- Step 2: Aggregate per-sample estimates to get final parameters ---
#'     brisc_param_avg <- apply(brisc_param_array, c(1, 2), mean, na.rm = TRUE)
#'     sigma2_B_data <- apply(brisc_param_array[, "mu", ], 1, var, na.rm = TRUE)
#'
#'     # Calculate hybrid sigma2_B using the new, more accurate sigma2_B_data
#'     Y_concat       <- do.call(rbind, Y_list)
#'     total_var_vec  <- apply(Y_concat, 2, var, na.rm = TRUE)
#'     sigma2_B_prior <- p_B_group * total_var_vec
#'     w_Kp           <- min(1, (K - 1) / 10)
#'     sigma2_B_hybrid <- (1 - w_Kp) * sigma2_B_prior + w_Kp * sigma2_B_data
#'
#'     # --- Save all outputs for this group ---
#'     results_by_group[[group]] <- list(
#'       brisc_param_avg  = brisc_param_avg,
#'       sigma2_B_data    = sigma2_B_data,
#'       sigma2_B_prior   = sigma2_B_prior,
#'       sigma2_B_hybrid  = sigma2_B_hybrid,
#'       sigma2_B_true    = rep(sigma_B_sq, p)
#'     )
#'   }
#'   return(results_by_group)
#' }
#'
#'
#'
#' estimate_parameters_by_group_K1 <- function(
#'     Y_list_by_group, groups, n, p, coords, p_B,
#'     effective_gene_indices, non_effective_gene_indices
#' ) {
#'   brisc_results_by_group <- list()
#'
#'   # --- Step 1: Fit BRISC on the single sample for each group ---
#'   for (group in groups) {
#'     Y_mat <- Y_list_by_group[[group]][[1]]
#'
#'     brisc_param_matrix <- matrix(NA, nrow = p, ncol = 4,
#'                                  dimnames = list(paste0("gene", 1:p),
#'                                                  c("mu", "alpha_spatial", "rho", "sigma2_noise")))
#'     n_neighbors = 15
#'     for (g in 1:p) {
#'       Y_g <- Y_mat[, g]
#'       if (var(Y_g, na.rm = TRUE) == 0) Y_g <- Y_g + rnorm(length(Y_g), 0, 1e-6)
#'
#'       brisc_fit <- suppressMessages(try(BRISC::BRISC_estimation(
#'         coords = coords, y = Y_g, x = matrix(1, nrow = n),
#'         cov.model = "exponential", n.neighbors = n_neighbors, order = "AMMD",
#'         verbose = FALSE, nugget_status = 1
#'       ), silent = TRUE))
#'
#'       if (!inherits(brisc_fit, "try-error")) {
#'         brisc_param_matrix[g, ] <- c(as.numeric(brisc_fit$Beta),
#'                                      brisc_fit$Theta["sigma.sq"],
#'                                      brisc_fit$Theta["phi"],
#'                                      brisc_fit$Theta["tau.sq"])
#'       }
#'     }
#'     brisc_results_by_group[[group]] <- brisc_param_matrix
#'   }
#'
#'   # --- Step 2: Post-hoc pooling of spatial parameters ---
#'   alpha_hat_pooled <- (brisc_results_by_group$case[, "alpha_spatial"] + brisc_results_by_group$control[, "alpha_spatial"]) / 2
#'   rho_hat_pooled   <- (brisc_results_by_group$case[, "rho"] + brisc_results_by_group$control[, "rho"]) / 2
#'
#'   final_results <- list()
#'   sigma2_B_hybrid_by_group <- list()
#'   for(c_idx in seq_along(groups)){
#'     group <- groups[c_idx]
#'     p_B_group <- p_B[c_idx]
#'     estimated_nuggets_g <- brisc_results_by_group[[group]][, "sigma2_noise"]
#'     # Average across all genes to get a single value, as per Algorithm 1, line 14
#'     sigma2_B_hybrid_by_group[[group]] <- mean(p_B_group * estimated_nuggets_g, na.rm = TRUE)
#'   }
#'
#'   # Now, build the final brisc_param_avg for each group
#'   for (c_idx in seq_along(groups)) {
#'     group <- groups[c_idx]
#'
#'     # Start with the original BRISC estimates for this group
#'     brisc_param_avg <- brisc_results_by_group[[group]]
#'
#'     # Overwrite with pooled spatial parameters
#'     brisc_param_avg[, "alpha_spatial"] <- alpha_hat_pooled
#'     brisc_param_avg[, "rho"] <- rho_hat_pooled
#'
#'     # The `sigma2_B_hybrid` is now the single, group-level value calculated above
#'     sigma2_B_hybrid_current_group <- sigma2_B_hybrid_by_group[[group]]
#'
#'     final_results[[group]] <- list(
#'       brisc_param_avg = brisc_param_avg,
#'       sigma2_B_hybrid = sigma2_B_hybrid_current_group # This is now a single value
#'     )
#'   }
#'
#'   return(final_results)
#' }
#'
#'
#'
#'
#' get_param_table <- function(results_by_group, groups, p) {
#'   final_param_estimates <- do.call(rbind, lapply(groups, function(group) {
#'     brisc_param_avg <- results_by_group[[group]]$brisc_param_avg
#'     sigma2_B_hybrid <- results_by_group[[group]]$sigma2_B_hybrid
#'     data.frame(
#'       Gene = paste0("Gene", 1:p),
#'       mu_est = brisc_param_avg[, "mu"],
#'       alpha_spatial_est = brisc_param_avg[, "alpha_spatial"],
#'       rho_est = brisc_param_avg[, "rho"],
#'       sigma2_B_hybrid = sigma2_B_hybrid,
#'       sigma2_noise_est = brisc_param_avg[, "sigma2_noise"],
#'       Group = group
#'     )
#'   }))
#'   rownames(final_param_estimates) <- NULL
#'   return(final_param_estimates)
#' }
#'
#'
#' get_sample_means <- function(Y_list) {
#'   # Y_list: list of K matrices, each n x p
#'   K <- length(Y_list)
#'   p <- ncol(Y_list[[1]])
#'   means <- matrix(NA, nrow = p, ncol = K)
#'   for (k in 1:K) {
#'     means[, k] <- colMeans(Y_list[[k]])
#'   }
#'   rownames(means) <- paste0("Gene", 1:p)
#'   colnames(means) <- paste0("Sample", 1:K)
#'   means
#' }
#'
#'
#' get_pooled_variances <- function(Y_list) {
#'   p <- ncol(Y_list[[1]])
#'   s2_list <- lapply(Y_list, function(sample_matrix) {
#'     apply(sample_matrix, 2, var, na.rm = TRUE)
#'   })
#'   s2_matrix <- do.call(rbind, s2_list) # K x p matrix of variances
#'   pooled_s2 <- colMeans(s2_matrix, na.rm = TRUE) # A single p-length vector
#'   return(pooled_s2)
#' }
#'
#'
#'
#'
#' make_alternative_from_true_beta <- function(beta_mat, tol = 0) {
#'   # beta_mat: 2열(case, control)을 가진 행렬/데이터프레임 (행=유전자)
#'   # tol: 두 평균 차이가 tol 이하이면 two.sided로 처리 (기본 0: 부호만 보고 결정)
#'   if (!all(c("case","control") %in% colnames(beta_mat)))
#'     stop("beta_mat must have columns named 'case' and 'control'.")
#'
#'   diff <- beta_mat[, "case"] - beta_mat[, "control"]
#'   alt  <- character(length(diff))
#'   alt[diff >  tol] <- "greater"
#'   alt[diff < -tol] <- "less"
#'   alt[abs(diff) <= tol] <- "two.sided"
#'   alt
#' }
#'
#'
#' moderated_ttest <- function(
#'     Y_case, Y_control,
#'     gene_effective = NULL,
#'     alpha = 0.05,
#'     alternative = "auto",   # "auto" | "greater" | "less" | "two.sided" | length-p vector of those
#'     direction = "greater"   # used only when alternative == "auto": one-sided direction for effective genes
#' ) {
#'   # --- 0) helpers ---
#'   .match_alt <- function(x) {
#'     x <- tolower(x)
#'     if (!x %in% c("greater", "less", "two.sided"))
#'       stop("alternative must be one of: 'greater','less','two.sided'")
#'     x
#'   }
#'
#'   # --- 1. Compute per-sample means ---
#'   case_means    <- get_sample_means(Y_case)
#'   control_means <- get_sample_means(Y_control)
#'
#'   p <- nrow(case_means)
#'   K_case    <- ncol(case_means)
#'   K_control <- ncol(control_means)
#'
#'   # --- 2. Calculate robust, pooled sample variances ---
#'   s2_case_all    <- get_pooled_variances(Y_case)
#'   s2_control_all <- get_pooled_variances(Y_control)
#'
#'   # --- 3. Moderation (Smyth 2004-style moment matching) ---
#'   bar_s2_case    <- mean(s2_case_all, na.rm = TRUE)
#'   bar_s2_control <- mean(s2_control_all, na.rm = TRUE)
#'   v_case         <- var(s2_case_all, na.rm = TRUE)
#'   v_control      <- var(s2_control_all, na.rm = TRUE)
#'   d0_case        <- if (!is.na(v_case) && v_case > 0) 2 * bar_s2_case^2 / v_case else 10
#'   d0_control     <- if (!is.na(v_control) && v_control > 0) 2 * bar_s2_control^2 / v_control else 10
#'
#'   # --- 4. Build per-gene 'alternative' vector ---
#'   if (length(alternative) == 1L) {
#'     alt_mode <- tolower(alternative)
#'     if (alt_mode == "auto") {
#'       # effective genes: one-sided with 'direction', others: two-sided
#'       if (is.null(gene_effective)) {
#'         # no prior -> default to two-sided everywhere
#'         alt_vec <- rep("two.sided", p)
#'       } else {
#'         alt_vec <- rep("two.sided", p)
#'         direction <- .match_alt(direction)
#'         if (length(gene_effective)) {
#'           idx <- as.integer(gene_effective)
#'           idx <- idx[idx >= 1 & idx <= p]
#'           alt_vec[idx] <- direction
#'         }
#'       }
#'     } else {
#'       # single mode for all genes
#'       alt_vec <- rep(.match_alt(alt_mode), p)
#'     }
#'   } else {
#'     # user supplied per-gene vector
#'     if (length(alternative) != p)
#'       stop("When providing per-gene 'alternative', its length must equal p (number of genes).")
#'     alt_vec <- tolower(alternative)
#'     ok <- alt_vec %in% c("greater","less","two.sided")
#'     if (!all(ok)) stop("Invalid entries in 'alternative'. Use only 'greater','less','two.sided'.")
#'   }
#'
#'   # --- 5. Loop through genes to calculate moderated t-statistics & p-values ---
#'   results <- lapply(seq_len(p), function(g) {
#'     # Per-gene variances
#'     var_case    <- s2_case_all[g]
#'     var_control <- s2_control_all[g]
#'
#'     # Posterior shrinkage variance for this gene
#'     s2_tilde_case    <- (d0_case * bar_s2_case    + (K_case    - 1) * var_case   ) / (d0_case    + K_case    - 1)
#'     s2_tilde_control <- (d0_control * bar_s2_control + (K_control - 1) * var_control) / (d0_control + K_control - 1)
#'
#'     # Moderated t-statistic
#'     mean_case_g     <- mean(case_means[g, ])
#'     mean_control_g  <- mean(control_means[g, ])
#'     mean_diff       <- mean_case_g - mean_control_g
#'     denom           <- sqrt(s2_tilde_case / K_case + s2_tilde_control / K_control)
#'
#'     # guard against zero/NA denom
#'     if (!is.finite(denom) || denom <= 0) {
#'       t_mod <- NA_real_
#'       df_mod <- NA_real_
#'       p_mod <- NA_real_
#'     } else {
#'       t_mod <- mean_diff / denom
#'
#'       # Welch-Satterthwaite d.f. under moderation
#'       num   <- (s2_tilde_case / K_case + s2_tilde_control / K_control)^2
#'       denom_df <- ((s2_tilde_case / K_case)^2    / (d0_case    + K_case    - 1)) +
#'         ((s2_tilde_control / K_control)^2 / (d0_control + K_control - 1))
#'       df_mod <- num / denom_df
#'
#'       # Per-gene alternative
#'       alt_g <- alt_vec[g]
#'       if (alt_g == "greater") {
#'         p_mod <- stats::pt(t_mod, df = df_mod, lower.tail = FALSE)
#'       } else if (alt_g == "less") {
#'         p_mod <- stats::pt(t_mod, df = df_mod, lower.tail = TRUE)
#'       } else {
#'         # two-sided
#'         p_mod <- 2 * stats::pt(-abs(t_mod), df = df_mod)
#'       }
#'     }
#'
#'     c(mean_case = mean_case_g,
#'       mean_control = mean_control_g,
#'       mean_diff = mean_diff,
#'       sd = denom,
#'       t_moderated = t_mod,
#'       p_moderated = p_mod,
#'       df_moderated = df_mod)
#'   })
#'
#'   # --- 6. Combine results and add adjusted p-values ---
#'   t_test_results <- as.data.frame(do.call(rbind, results))
#'   t_test_results <- cbind(Gene = paste0("Gene", seq_len(p)),
#'                           Alternative = alt_vec,
#'                           t_test_results)
#'
#'   t_test_results$q_value    <- p.adjust(t_test_results$p_moderated, method = "BH")
#'   t_test_results$p_value_bf <- p.adjust(t_test_results$p_moderated, method = "bonferroni")
#'
#'   t_test_results$significant_fdr  <- t_test_results$q_value    < alpha
#'   t_test_results$significant_fwer <- t_test_results$p_value_bf < alpha
#'
#'   t_test_results
#' }
#'
#'
#'
#'
#' run_wilcoxon_test_per_gene <- function(case_means, control_means) {
#'   p <- nrow(case_means)
#'   p_values <- numeric(p)
#'
#'   for (g in 1:p) {
#'     x_case    <- as.numeric(case_means[g, ])
#'     x_control <- as.numeric(control_means[g, ])
#'
#'     # Skip if data has no variance (e.g., all identical values)
#'     if (var(x_case, na.rm = TRUE) == 0 && var(x_control, na.rm = TRUE) == 0) {
#'       p_values[g] <- 1.0
#'       next
#'     }
#'
#'     # Set the alternative hypothesis based on the gene index
#'     alt_hypothesis <- "two.sided"
#'     if (g %in% gene_effective) { # H1: case > control for effective genes
#'       alt_hypothesis <- "greater"
#'     }
#'
#'     # Perform the Wilcoxon Rank-Sum test for the current gene
#'     test_result <- wilcox.test(x_case, x_control, alternative = alt_hypothesis)
#'     p_values[g] <- test_result$p.value
#'   }
#'
#'   return(p_values)
#' }
#'
#'
#' compute_directional_neighbors <- function(loc_df, radius = 0.1) {
#'   coords <- as.matrix(loc_df[, c("x", "y")])
#'   domain <- as.character(loc_df$domain)
#'   unique_domains <- sort(unique(domain))
#'
#'   # 결과를 저장할 리스트를 초기화합니다.
#'   domain_neighbors <- setNames(vector("list", length(unique_domains)), unique_domains)
#'
#'   # 각 도메인의 x 좌표 중앙값을 미리 계산해 둡니다.
#'   domain_x_medians <- tapply(coords[, "x"], domain, median)
#'
#'   for (i in seq_along(unique_domains)) {
#'     d1 <- unique_domains[i]
#'     idx1 <- which(domain == d1)
#'
#'     # 다른 모든 도메인에 대해 반복합니다.
#'     for (j in seq_along(unique_domains)) {
#'       d2 <- unique_domains[j]
#'       if (d1 == d2) next
#'
#'       idx2 <- which(domain == d2)
#'
#'       # 두 도메인 간의 최소 거리를 계산하여 '인접' 여부를 판단합니다.
#'       dist_mat <- as.matrix(stats::dist(rbind(coords[idx1, , drop=FALSE], coords[idx2, , drop=FALSE])))
#'       dist_sub <- dist_mat[1:length(idx1), (length(idx1)+1):(length(idx1)+length(idx2))]
#'
#'       if (any(dist_sub <= radius)) {
#'         # 인접한 경우, 좌우 구분 없이 리스트에 추가합니다.
#'         domain_neighbors[[d1]] <- c(domain_neighbors[[d1]], d2)
#'       }
#'     }
#'   }
#'
#'   # 최종적으로 각 리스트의 중복을 제거하고 정렬합니다.
#'   for (d in unique_domains) {
#'     domain_neighbors[[d]] <- sort(unique(domain_neighbors[[d]]))
#'   }
#'
#'   return(domain_neighbors)
#' }
#'
#'
#' compute_domain_neighbors <- function(loc_df, radius = 0.1) {
#'   coords <- as.matrix(loc_df[, c("x", "y")])
#'   domain <- as.character(loc_df$domain)
#'   unique_domains <- sort(unique(domain))
#'
#'   domain_neighbors <- setNames(vector("list", length(unique_domains)), unique_domains)
#'
#'   for (i in seq_along(unique_domains)) {
#'     d1 <- unique_domains[i]
#'     idx1 <- which(domain == d1)
#'     pts1 <- coords[idx1, , drop = FALSE]
#'
#'     neighbor_set <- c()
#'
#'     for (j in seq_along(unique_domains)) {
#'       d2 <- unique_domains[j]
#'       if (d1 == d2) next
#'       idx2 <- which(domain == d2)
#'       pts2 <- coords[idx2, , drop = FALSE]
#'
#'       dist_mat <- as.matrix(stats::dist(rbind(pts1, pts2)))
#'       dist_sub <- dist_mat[1:nrow(pts1), (nrow(pts1)+1):(nrow(pts1)+nrow(pts2))]
#'
#'       if (any(dist_sub <= radius)) {
#'         neighbor_set <- c(neighbor_set, d2)
#'       }
#'     }
#'
#'     domain_neighbors[[d1]] <- sort(unique(neighbor_set))
#'   }
#'
#'   return(domain_neighbors)
#' }
#'
#' generate_pattern <- function(location_df,
#'                              perturb_domains = NULL,
#'                              perturb_probs = NULL,
#'                              seed = 123,
#'                              boundary_radius = 0.05,
#'                              boundary_threshold = 0.5) {
#'   set.seed(seed)
#'
#'   loc_df <- location_df
#'   coords <- as.matrix(loc_df[, c("x", "y")])
#'   domain <- as.character(loc_df$domain)
#'   original_domains <- sort(unique(domain))
#'
#'   domain_neighbors <- compute_directional_neighbors(loc_df, radius = boundary_radius)
#'   domain_centroids <- aggregate(cbind(x, y) ~ domain, data = loc_df, FUN = median)
#'
#'   for (d_char in sort(original_domains, decreasing = TRUE)) {
#'     # Skip domains not in the perturbation list, if specified.
#'     if (!is.null(perturb_domains) && !(d_char %in% perturb_domains)) next
#'
#'     idx_d <- which(domain == d_char)
#'     n_d <- length(idx_d)
#'
#'     # Probability to perturb the current domain
#'     prob_d <- if (!is.null(perturb_probs)) {
#'       perturb_probs[[as.character(d_char)]]
#'     } else {
#'       0.05
#'     }
#'     if (is.na(prob_d) || prob_d <= 0) next
#'     n_change <- floor(n_d * prob_d)
#'
#'     # *** Updated Logic for Multiple Target Domains ***
#'     neighbors_of_d <- domain_neighbors[[d_char]]
#'     if (length(neighbors_of_d) == 0) next
#'
#'     # Identify all boundary points of the current domain
#'     is_boundary_point <- function(i) {
#'       pt <- matrix(coords[i, ], nrow = 1)
#'       dists <- sqrt(rowSums(sweep(coords, 2, pt)^2))
#'       nearby_idx <- which(dists > 0 & dists <= boundary_radius)
#'       if (length(nearby_idx) == 0) return(FALSE)
#'
#'       # A boundary point is now one that is near ANY neighboring domain.
#'       mean(domain[nearby_idx] %in% neighbors_of_d) >= boundary_threshold
#'     }
#'
#'     boundary_idx <- idx_d[sapply(idx_d, is_boundary_point)]
#'     if (length(boundary_idx) == 0) next
#'
#'     # Calculate distances from boundary points to ALL neighboring domains
#'     boundary_dists_matrix <- sapply(neighbors_of_d, function(target_d) {
#'       sapply(boundary_idx, function(i) {
#'         pt <- matrix(coords[i, ], nrow = 1)
#'         target_coords <- coords[domain == target_d, , drop = FALSE]
#'         min(sqrt(rowSums(sweep(target_coords, 2, pt)^2)), na.rm = TRUE)
#'       })
#'     })
#'
#'     # Sort points based on their minimum distance to *any* neighbor
#'     min_dists <- apply(boundary_dists_matrix, 1, min)
#'     sorted_idx <- boundary_idx[order(min_dists)]
#'
#'     # Select indices to be changed
#'     change_idx <- head(sorted_idx, n_change)
#'
#'     if (length(change_idx) > 0) {
#'       # Get the centroids of all neighboring domains
#'       neighbor_centroids <- domain_centroids[domain_centroids$domain %in% neighbors_of_d, ]
#'
#'       # For each point to be changed, find the closest neighbor centroid
#'       for (i in change_idx) {
#'         point_coords <- coords[i, ]
#'         distances_to_neighbors <- sqrt(rowSums(sweep(neighbor_centroids[, c("x", "y")], 2, point_coords)^2))
#'         closest_neighbor_index <- which.min(distances_to_neighbors)
#'         new_domain <- as.character(neighbor_centroids$domain[closest_neighbor_index])
#'
#'         domain[i] <- new_domain
#'       }
#'     }
#'   }
#'   loc_df$domain <- factor(domain, levels = original_domains)
#'   loc_df <- loc_df[order(loc_df$domain), ]
#'   rownames(loc_df) <- NULL
#'   return(loc_df)
#' }
#'
#'
#' library(matrixStats)
#'
#' # --- 1. Global Composition Test (Dirichlet-Multinomial LRT) ---
#'
#' #' Estimate Dirichlet Alpha using Method of Moments
#' estimate_dirichlet_alpha_mom <- function(counts, epsilon=1e-2) {
#'   # counts: K x D matrix
#'   # Smooth proportions
#'   props <- (counts + epsilon) / rowSums(counts + epsilon * ncol(counts))
#'
#'   pi_bar <- colMeans(props)
#'   v_bar <- colVars(props)
#'   v_bar[v_bar < 1e-6] <- 1e-6 # Avoid division by zero
#'
#'   # Estimate precision alpha0
#'   # alpha0 = (E[p]*(1-E[p]) / Var(p)) - 1
#'   alpha0_ests <- (pi_bar * (1 - pi_bar)) / v_bar - 1
#'   alpha0 <- mean(alpha0_ests, na.rm = TRUE)
#'
#'   # 1. Stability Floor (하한선): 0 이하가 되지 않도록 방지
#'   alpha0 <- max(alpha0, 1e-2)
#'
#'   # 2. Variance Protection Cap (상한선): 분산이 과소추정되어 alpha0가 폭발하는 것 방지
#'   # alpha0가 100을 넘어가면 분산이 거의 없다고 가정하게 되어 Type I Error가 급증합니다.
#'   alpha0 <- min(alpha0, 100)
#'
#'   return(alpha0 * pi_bar)
#' }
#'
#' #' Dirichlet-Multinomial Log-Likelihood
#' dm_loglik <- function(alpha, counts) {
#'   # alpha: vector of length D
#'   # counts: matrix K x D
#'   K <- nrow(counts)
#'   D <- ncol(counts)
#'   alpha0 <- sum(alpha)
#'
#'   # Term 1: K * logGamma(alpha0)
#'   term1 <- K * lgamma(alpha0)
#'
#'   # Term 2: sum_k logGamma(N_k + alpha0)
#'   term2 <- sum(lgamma(rowSums(counts) + alpha0))
#'
#'   # Term 3: sum_k sum_d logGamma(n_kd + alpha_d)
#'   # Use sweep to add alpha to each row of counts correctly
#'   term3 <- sum(lgamma(sweep(counts, 2, alpha, "+")))
#'
#'   # Term 4: K * sum_d logGamma(alpha_d)
#'   term4 <- K * sum(lgamma(alpha))
#'
#'   ll <- term1 - term2 + term3 - term4
#'   return(ll)
#' }
#'
#' #' Global Test for Compositional Difference
#' run_global_composition_test <- function(counts_case, counts_ctrl) {
#'   D <- ncol(counts_case)
#'
#'   # Null Model: Pool all samples
#'   counts_pooled <- rbind(counts_case, counts_ctrl)
#'   alpha_null <- estimate_dirichlet_alpha_mom(counts_pooled)
#'   ll_null <- dm_loglik(alpha_null, counts_pooled)
#'
#'   # Alternative Model: Separate groups
#'   alpha_case <- estimate_dirichlet_alpha_mom(counts_case)
#'   alpha_ctrl <- estimate_dirichlet_alpha_mom(counts_ctrl)
#'   ll_alt <- dm_loglik(alpha_case, counts_case) + dm_loglik(alpha_ctrl, counts_ctrl)
#'
#'   # LRT Statistic
#'   lrt_stat <- 2 * (ll_alt - ll_null)
#'
#'   # P-value
#'   # df = (2D - D) = D (assuming alpha0 is not shared/constrained between groups)
#'   p_val <- pchisq(lrt_stat, df = D, lower.tail = FALSE)
#'
#'   return(list(statistic = lrt_stat, p_value = p_val, df = D))
#' }
#'
#' # --- 2. Morphological Analysis (Feature Extraction & Joint Test) ---
#'
#' #' Extract Geometry Features (R, A) from synthetic data list
#' extract_geometry_features <- function(sample_list, D) {
#'   K <- length(sample_list)
#'   R_mat <- matrix(NA, K, D)
#'   A_mat <- matrix(NA, K, D)
#'
#'   for (k in 1:K) {
#'     df <- sample_list[[k]]
#'     if(!"domain" %in% colnames(df)) next
#'
#'     for (d in 1:D) {
#'       coords <- df[df$domain == d, c("x", "y")]
#'       if (nrow(coords) < 5) next
#'
#'       Sigma <- cov(coords)
#'       eig <- eigen(Sigma)$values
#'       lambda1 <- max(eig)
#'       lambda2 <- min(eig)
#'       if (lambda2 <= 0) lambda2 <- 1e-6
#'
#'       R_mat[k, d] <- 0.5 * log(lambda1 + lambda2)
#'       A_mat[k, d] <- log(lambda1 / lambda2)
#'     }
#'   }
#'   return(list(R = R_mat, A = A_mat))
#' }
#'
#'
#'
#' # Updated function to estimate prior hyperparameters via MLE with improved stability
#' # Argument name changed to 'df1' to match the call in run_moderated_test
#' estimate_prior_params_mle <- function(sample_variances, df1) {
#'
#'   # Remove NAs to prevent errors in optimization
#'   sample_variances <- sample_variances[!is.na(sample_variances)]
#'
#'   # The function to be minimized (Negative Log-Likelihood)
#'   neg_log_likelihood <- function(params) {
#'     s0_sq <- params[1]
#'     d0 <- params[2]
#'
#'     # Constraint check (d0 must be > 4 for finite variance, s0_sq > 0)
#'     if (s0_sq <= 0 || d0 <= 4.0001) {
#'       return(1e100)
#'     }
#'
#'     # Calculate the log-likelihood for a scaled F-distribution
#'     # Model: s^2 ~ s0^2 * F(df1, d0)
#'     # Variable transformation: y = s^2/s0^2 ~ F(df1, d0)
#'     # Log PDF = log(df(y, df1, df2)) - log(s0^2) [Jacobian]
#'
#'     vals <- sample_variances / s0_sq
#'
#'     # df1 is passed from the outer function scope
#'     log_f <- df(vals, df1 = df1, df2 = d0, log = TRUE)
#'     log_lik_terms <- log_f - log(s0_sq)
#'
#'     # Check for non-finite values (NA, NaN, Inf)
#'     if (any(!is.finite(log_lik_terms))) {
#'       return(1e100)
#'     }
#'
#'     # Return sum of negative log-likelihoods
#'     return(-sum(log_lik_terms))
#'   }
#'
#'   # --- Initialization via Method of Moments (MoM) ---
#'   initial_s0_sq <- mean(sample_variances)
#'   v_s2 <- var(sample_variances)
#'
#'   # Approximate d0 based on variance properties
#'   if (is.na(v_s2) || v_s2 <= 1e-8) {
#'     initial_d0 <- 20
#'   } else {
#'     # Simple heuristic for initial d0
#'     initial_d0 <- 2 * initial_s0_sq^2 / v_s2 + 4
#'     initial_d0 <- min(max(initial_d0, 5), 50)
#'   }
#'
#'   if (initial_s0_sq <= 0) initial_s0_sq <- 1e-3
#'
#'   initial_params <- c(initial_s0_sq, initial_d0)
#'
#'   # --- Optimization ---
#'   opt_result <- tryCatch({
#'     optim(
#'       par = initial_params,
#'       fn = neg_log_likelihood,
#'       method = "L-BFGS-B",
#'       lower = c(1e-6, 4.1), # d0 > 4 constraint
#'       upper = c(Inf, 1e5)
#'     )
#'   }, error = function(e) {
#'     list(par = initial_params)
#'   })
#'
#'   return(list(s0_sq = opt_result$par[1], d0 = opt_result$par[2]))
#' }
#'
#'
#' #' Joint Test for Geometric Differences (R & A)
#' run_joint_geometry_test <- function(feat_case, feat_ctrl, alpha = 0.05) {
#'   n_domains <- ncol(feat_case$R)
#'   results_df <- data.frame(domain = 1:n_domains, p_value = NA, T_R = NA, T_A = NA, Q = NA)
#'
#'   # Helper to compute simple stats
#'   calc_stats <- function(mat) {
#'     list(
#'       mean = colMeans(mat, na.rm=TRUE),
#'       var = colVars(mat, na.rm=TRUE),
#'       n = colSums(!is.na(mat))
#'     )
#'   }
#'
#'   stats_R_case <- calc_stats(feat_case$R)
#'   stats_R_ctrl <- calc_stats(feat_ctrl$R)
#'   stats_A_case <- calc_stats(feat_case$A)
#'   stats_A_ctrl <- calc_stats(feat_ctrl$A)
#'
#'   # Prior Estimation (Pooling variances across ALL domains to stabilize)
#'   # We pool case variances and ctrl variances separately
#'   prior_R_case <- estimate_prior_params_mle(stats_R_case$var, median(stats_R_case$n)-1)
#'   prior_R_ctrl <- estimate_prior_params_mle(stats_R_ctrl$var, median(stats_R_ctrl$n)-1)
#'   prior_A_case <- estimate_prior_params_mle(stats_A_case$var, median(stats_A_case$n)-1)
#'   prior_A_ctrl <- estimate_prior_params_mle(stats_A_ctrl$var, median(stats_A_ctrl$n)-1)
#'
#'   for (d in 1:n_domains) {
#'     # R Statistic
#'     df1 <- stats_R_case$n[d] - 1
#'     df2 <- stats_R_ctrl$n[d] - 1
#'     if(df1 < 1 || df2 < 1) next
#'
#'     s2_R_case_mod <- (prior_R_case$d0 * prior_R_case$s0_sq + df1 * stats_R_case$var[d]) / (prior_R_case$d0 + df1)
#'     s2_R_ctrl_mod <- (prior_R_ctrl$d0 * prior_R_ctrl$s0_sq + df2 * stats_R_ctrl$var[d]) / (prior_R_ctrl$d0 + df2)
#'
#'     se_R <- sqrt(s2_R_case_mod/stats_R_case$n[d] + s2_R_ctrl_mod/stats_R_ctrl$n[d])
#'     t_R <- (stats_R_case$mean[d] - stats_R_ctrl$mean[d]) / se_R
#'
#'     # A Statistic
#'     df1_A <- stats_A_case$n[d] - 1
#'     df2_A <- stats_A_ctrl$n[d] - 1
#'
#'     s2_A_case_mod <- (prior_A_case$d0 * prior_A_case$s0_sq + df1_A * stats_A_case$var[d]) / (prior_A_case$d0 + df1_A)
#'     s2_A_ctrl_mod <- (prior_A_ctrl$d0 * prior_A_ctrl$s0_sq + df2_A * stats_A_ctrl$var[d]) / (prior_A_ctrl$d0 + df2_A)
#'
#'     se_A <- sqrt(s2_A_case_mod/stats_A_case$n[d] + s2_A_ctrl_mod/stats_A_ctrl$n[d])
#'     t_A <- (stats_A_case$mean[d] - stats_A_ctrl$mean[d]) / se_A
#'
#'     # Joint Q
#'     Q <- t_R^2 + t_A^2
#'     p_val <- pchisq(Q, df = 2, lower.tail = FALSE)
#'
#'     results_df$p_value[d] <- p_val
#'     results_df$T_R[d] <- t_R
#'     results_df$T_A[d] <- t_A
#'     results_df$Q[d] <- Q
#'   }
#'
#'   results_df$q_value <- p.adjust(results_df$p_value, method = "BH")
#'   results_df$rejected <- results_df$q_value < alpha
#'
#'   return(results_df)
#' }
#'
#' # ==============================================================================
#' # 1. Helper Function: Process a Single Visium Sample Folder (Final Optimized)
#' # ==============================================================================
#' #' Process a single Visium sample folder containing .RData and meta.csv
#' #'
#' #' Updates:
#' #' 1. Calculates log-transformed counts (log1p).
#' #' 2. Normalizes spatial coordinates to [0, 1] range.
#' #'
#' #' @param folder_path Path to the sample folder
#' #' @param group_label Binary group label (0 for Control, 1 for Case)
#' #' @param domain_col Column name in meta.csv for domain/cluster labels
#' #' @return A list containing counts, logcounts, coords, group, and sample_id
#' process_single_visium <- function(folder_path, group_label, domain_col="benmarklabel") {
#'
#'   # Extract sample ID from the folder name
#'   sample_id <- basename(folder_path)
#'   message(sprintf("Processing Sample: %s (Group: %s)...", sample_id, group_label))
#'
#'   # 1-1) Define file paths
#'   rdata_path <- file.path(folder_path, "brain_processed.RData")
#'   meta_csv_path <- file.path(folder_path, "meta.csv")
#'
#'   # Check if files exist
#'   if(!file.exists(rdata_path) || !file.exists(meta_csv_path)) {
#'     warning(paste("Missing files in:", folder_path))
#'     return(NULL)
#'   }
#'
#'   # 1-2) Load .RData safely
#'   e <- new.env()
#'   load(rdata_path, envir = e)
#'
#'   obj_name <- ls(e)[1]
#'   brain <- e[[obj_name]]
#'
#'   if(!inherits(brain, "Seurat")) {
#'     warning(paste("Object in", rdata_path, "is not a Seurat object."))
#'     return(NULL)
#'   }
#'
#'   # 1-3) Extract Count Matrix
#'   assay_use <- "Spatial"
#'   if(!"Spatial" %in% names(brain@assays)) assay_use <- DefaultAssay(brain)
#'
#'   count_matrix <- GetAssayData(brain, assay = assay_use, layer = "counts")
#'
#'   # 1-4) Extract Coordinates
#'   if(length(brain@images) == 0) {
#'     warning(paste("No images found in Seurat object for sample", sample_id))
#'     return(NULL)
#'   }
#'   image_name <- names(brain@images)[1]
#'   coords <- brain@images[[image_name]]@coordinates
#'
#'   loc0 <- data.frame(
#'     x = as.numeric(coords$col),
#'     y = as.numeric(coords$row),
#'     stringsAsFactors = FALSE
#'   )
#'   rownames(loc0) <- rownames(coords)
#'
#'   # =========================================================
#'   # [NEW] Coordinate Normalization to [0, 1]
#'   # =========================================================
#'   # Normalize X
#'   x_min <- min(loc0$x); x_max <- max(loc0$x)
#'   if (x_max > x_min) loc0$x <- (loc0$x - x_min) / (x_max - x_min)
#'
#'   # Normalize Y
#'   y_min <- min(loc0$y); y_max <- max(loc0$y)
#'   if (y_max > y_min) loc0$y <- (loc0$y - y_min) / (y_max - y_min)
#'   # =========================================================
#'
#'   # 1-5) Load Meta CSV and Map Domains
#'   meta <- read_csv(meta_csv_path, show_col_types = FALSE)
#'   barcode_col_csv <- "barcode"
#'
#'   if (!all(c(barcode_col_csv, domain_col) %in% colnames(meta))) {
#'     stop("meta.csv is missing required columns: barcode or domain column.")
#'   }
#'
#'   # Handle Barcode Prefixes
#'   meta <- meta %>%
#'     mutate(
#'       barcode = .data[[barcode_col_csv]],
#'       domain  = .data[[domain_col]],
#'       spot_id_plain  = as.character(barcode),
#'       spot_id_prefix = paste0(sample_id, "_", as.character(barcode))
#'     )
#'
#'   cn <- colnames(count_matrix)
#'   use_prefixed <- FALSE
#'   if (all(grepl("^\\d+_", cn))) {
#'     use_prefixed <- TRUE
#'   }
#'
#'   meta_spot_id <- if (use_prefixed) meta$spot_id_prefix else meta$spot_id_plain
#'   domain_map <- setNames(as.character(meta$domain), meta_spot_id)
#'
#'   if (use_prefixed && !all(grepl("^\\d+_", rownames(loc0)))) {
#'     rownames(loc0) <- paste0(sample_id, "_", rownames(loc0))
#'   }
#'
#'   loc0$domain <- domain_map[rownames(loc0)]
#'
#'   # 1-6) Intersect and Clean Data
#'   common <- intersect(colnames(count_matrix), rownames(loc0))
#'
#'   if (length(common) == 0) {
#'     warning(paste("No matching spots for sample", sample_id))
#'     return(NULL)
#'   }
#'
#'   count_matrix <- count_matrix[, common, drop = FALSE]
#'   loc_file <- loc0[common, , drop = FALSE]
#'
#'   keep <- !is.na(loc_file$domain)
#'   count_matrix <- count_matrix[, keep, drop = FALSE]
#'   loc_file <- loc_file[keep, , drop = FALSE]
#'   loc_file$domain <- factor(loc_file$domain)
#'
#'   if(!identical(colnames(count_matrix), rownames(loc_file))) {
#'     stop("Mismatch between count matrix columns and coordinate rows.")
#'   }
#'
#'   # =========================================================
#'   # [NEW] Calculate Log-Transformed Counts (log1p)
#'   # =========================================================
#'   # Safe log transformation (handles sparse matrices via S4 dispatch)
#'   log_matrix <- log1p(count_matrix)
#'
#'   # Return structured list including logcounts and normalized coords
#'   return(list(
#'     counts = count_matrix,    # Raw counts
#'     logcounts = log_matrix,   # Log-transformed counts
#'     coords = loc_file,        # [0,1] Normalized coords
#'     group = group_label,
#'     sample_id = sample_id
#'   ))
#' }
#'
#' get_domain_markers <- function(object, top_n = 2) {
#'   # 1. Use the first pilot sample as reference
#'   rep_sample <- object@pilot_data[[1]]
#'   count_mat  <- rep_sample$counts
#'   coords     <- rep_sample$coords
#'
#'   # 2. Log-transform
#'   log_counts <- log1p(as.matrix(count_mat))
#'   domains    <- sort(unique(coords$domain))
#'
#'   marker_list <- list()
#'
#'   # 3. Iterate over domains to find specific markers
#'   message(">>> Identifying top ", top_n, " markers per domain...")
#'   for(dom in domains) {
#'     # Identify cells inside vs outside the domain
#'     in_idx  <- which(coords$domain == dom)
#'
#'     # Calculate Mean In vs Mean Out for all genes
#'     # (RowMeans is fast enough for typical gene sets)
#'     mean_in  <- rowMeans(log_counts[, in_idx, drop=FALSE])
#'     mean_out <- rowMeans(log_counts[, -in_idx, drop=FALSE])
#'
#'     # Log Fold Change
#'     lfc <- mean_in - mean_out
#'
#'     # Select genes with highest positive LFC (specific to this domain)
#'     # We create a dataframe to sort easily
#'     df <- data.frame(gene = rownames(log_counts), lfc = lfc)
#'     top_genes <- df %>%
#'       arrange(desc(lfc)) %>%
#'       head(top_n) %>%
#'       pull(gene)
#'
#'     marker_list[[dom]] <- top_genes
#'     cat(sprintf("   Domain %s: %s\n", dom, paste(top_genes, collapse=", ")))
#'   }
#'
#'   # 4. Return unique union of these top markers
#'   return(unique(unlist(marker_list)))
#' }
#'
#'
#' # ------------------------------------------------------------
#' # 0) 샘플별 z_hat = mean(T) - mean(R) 만들기 (pilot only)
#' # ------------------------------------------------------------
#' build_z_from_pilot <- function(object, genes, target_domain="WM", reference_domain="Layer6") {
#'   pilot <- object@pilot_data
#'   rows <- list(); k <- 0L
#'
#'   for (i in seq_along(pilot)) {
#'     samp <- pilot[[i]]
#'     coords <- samp$coords
#'     X <- samp$logcounts
#'     grp <- as.numeric(samp$group)
#'
#'     idx_T <- which(coords$domain == target_domain)
#'     idx_R <- which(coords$domain == reference_domain)
#'     if (length(idx_T) < 2 || length(idx_R) < 2) next
#'
#'     for (g in genes) {
#'       if (!g %in% rownames(X)) next
#'       z <- mean(as.numeric(X[g, idx_T])) - mean(as.numeric(X[g, idx_R]))
#'       k <- k + 1L
#'       rows[[k]] <- data.frame(gene=g, sample=i, group=grp, z_hat=z)
#'     }
#'   }
#'   do.call(rbind, rows)
#' }
#'
#' # ------------------------------------------------------------
#' # 1) Target/Ref에서 “검정 가능한” gene pool 만들기 (발현/검출률 필터)
#' # ------------------------------------------------------------
#' make_TR_pool <- function(object, target_domain="WM", reference_domain="Layer6",
#'                          sample_index=1, mean_cut=0.5, pct_cut=0.05) {
#'   samp <- object@pilot_data[[sample_index]]
#'   X <- samp$logcounts
#'   coords <- samp$coords
#'   idx_T <- which(coords$domain == target_domain)
#'   idx_R <- which(coords$domain == reference_domain)
#'   stopifnot(length(idx_T) > 0, length(idx_R) > 0)
#'
#'   mean_T <- rowMeans(X[, idx_T, drop=FALSE])
#'   mean_R <- rowMeans(X[, idx_R, drop=FALSE])
#'   pct_T  <- rowMeans(X[, idx_T, drop=FALSE] > 0)
#'   pct_R  <- rowMeans(X[, idx_R, drop=FALSE] > 0)
#'
#'   df <- data.frame(
#'     gene = rownames(X),
#'     mean_T = mean_T, mean_R = mean_R,
#'     pct_T = pct_T, pct_R = pct_R,
#'     TR = mean_T - mean_R,
#'     stringsAsFactors = FALSE
#'   )
#'
#'   pool <- df$gene[df$mean_T >= mean_cut & df$mean_R >= mean_cut &
#'                     df$pct_T  >= pct_cut  & df$pct_R  >= pct_cut]
#'   list(pool=pool, stats=df)
#' }
#'
#' # ------------------------------------------------------------
#' # 2) pilot만으로 SaLFC용 G_test, G_DE 구성
#' #    - group 차이 의심 유전자 제외
#' #    - G_DE는 중간 난이도에서 선택
#' # ------------------------------------------------------------
#' make_gene_sets_pilot_only <- function(object,
#'                                       target_domain="WM", reference_domain="Layer6",
#'                                       # 발현 필터
#'                                       mean_cut=0.5, pct_cut=0.05,
#'                                       # "pilot에서 group 차이 의심" 제외 기준
#'                                       p_exclude = 0.2,
#'                                       # G_DE 선택 난이도 (TR 기준 분위수)
#'                                       de_tr_q = c(0.4, 0.8),
#'                                       n_de = 3,
#'                                       n_test = 10,
#'                                       verbose=TRUE) {
#'
#'   # (A) TR pool
#'   TR <- make_TR_pool(object, target_domain, reference_domain,
#'                      sample_index=1, mean_cut=mean_cut, pct_cut=pct_cut)
#'   pool <- TR$pool
#'   st <- TR$stats
#'
#'   # (B) pilot에서 z_hat 만들어 그룹 차이 검사
#'   z_long <- build_z_from_pilot(object, pool, target_domain, reference_domain)
#'
#'   p_by_gene <- tapply(seq_len(nrow(z_long)), z_long$gene, function(ii) {
#'     df <- z_long[ii, ]
#'     if (length(unique(df$group)) < 2) return(c(p_t=NA_real_, p_w=NA_real_))
#'
#'     p_t <- tryCatch(t.test(z_hat ~ factor(group), data=df)$p.value, error=function(e) NA_real_)
#'     p_w <- tryCatch(wilcox.test(z_hat ~ factor(group), data=df, exact=FALSE)$p.value, error=function(e) NA_real_)
#'
#'     c(p_t=p_t, p_w=p_w)
#'   })
#'
#'   p_mat <- do.call(rbind, p_by_gene)
#'
#'   # 보수적 제외: 둘 중 하나라도 작으면 제외
#'   keep <- rownames(p_mat)[is.finite(p_mat[,"p_t"]) & is.finite(p_mat[,"p_w"]) &
#'                             p_mat[,"p_t"] >= p_exclude & p_mat[,"p_w"] >= p_exclude]
#'
#'   pool2 <- intersect(pool, keep)
#'
#'   if (verbose) {
#'     message("TR pool size: ", length(pool))
#'     message("After excluding suspicious group-diff genes (p<", p_exclude, "): ", length(pool2))
#'   }
#'   if (length(pool2) < n_test) warning("pool2 is small; relax mean_cut/pct_cut or p_exclude.")
#'
#'   # (C) G_DE: TR의 '중간~상대적으로 큰' 구간에서 선택 (너무 쉬운/어려운 것 피함)
#'   st2 <- st[st$gene %in% pool2, , drop=FALSE]
#'   absTR <- abs(st2$TR)
#'   qlo <- stats::quantile(absTR, probs=de_tr_q[1], na.rm=TRUE)
#'   qhi <- stats::quantile(absTR, probs=de_tr_q[2], na.rm=TRUE)
#'
#'   de_pool <- st2$gene[absTR >= qlo & absTR <= qhi]
#'   # 발현이 너무 낮은 것 추가 배제(선택사항)
#'   de_pool <- de_pool[order(-abs(st2$TR[match(de_pool, st2$gene)]))]
#'
#'   G_DE <- head(de_pool, n_de)
#'
#'   # (D) G_test: G_DE 포함 + 나머지는 "TR이 너무 크지 않은" 쪽에서 채움(=null로 남기기 쉬움)
#'   remain <- setdiff(pool2, G_DE)
#'   # null 후보는 abs(TR) 작은 쪽 우선
#'   remain <- remain[order(abs(st2$TR[match(remain, st2$gene)]))]
#'   G_null <- head(remain, max(0, n_test - length(G_DE)))
#'
#'   G_test <- c(G_DE, G_null)
#'
#'   list(
#'     G_TR_pool = pool,
#'     G_test = G_test,
#'     G_DE = G_DE,
#'     TR_stats = st2,
#'     p_groupdiff = p_by_gene
#'   )
#' }
#'
#'
#'
#' # -------------------------------------------------------------------------
#' # Helper: safe normalization
#' # -------------------------------------------------------------------------
#' .safe_unit_vec <- function(v, fallback = c(1, 0)) {
#'   fallback <- as.numeric(fallback)
#'   v <- as.numeric(v)
#'   if (length(v) != length(fallback) || any(!is.finite(v))) return(fallback)
#'   nv <- sqrt(sum(v^2))
#'   if (!is.finite(nv) || nv <= 1e-12) return(fallback)
#'   v / nv
#' }
#'
#' # -------------------------------------------------------------------------
#' # Helper: safe 2D Gaussian draw
#' # -------------------------------------------------------------------------
#' .safe_mvrnorm_vec <- function(mu, Sigma, jitter = 1e-8) {
#'   mu <- as.numeric(mu)
#'   p <- length(mu)
#'
#'   if (p == 0L || any(!is.finite(mu))) {
#'     mu <- rep(0, max(p, 2L))
#'     p <- length(mu)
#'   }
#'
#'   Sigma <- as.matrix(Sigma)
#'   if (!all(dim(Sigma) == c(p, p)) || any(!is.finite(Sigma))) {
#'     Sigma <- diag(jitter, p)
#'   }
#'
#'   Sigma <- (Sigma + t(Sigma)) / 2
#'   diag(Sigma) <- diag(Sigma) + jitter
#'
#'   out <- tryCatch(
#'     MASS::mvrnorm(1, mu = mu, Sigma = Sigma),
#'     error = function(e) mu
#'   )
#'   as.numeric(out)
#' }
#'
#' # -------------------------------------------------------------------------
#' # Helper: vMF draw with fallback
#' #   - for very small kappa, random unit direction is more natural than fixed mu
#' # -------------------------------------------------------------------------
#' .safe_rvmf1 <- function(mu, kappa, fallback = c(1, 0)) {
#'   mu <- .safe_unit_vec(mu, fallback = fallback)
#'
#'   if (!is.finite(kappa) || kappa <= 1e-8) {
#'     ang <- stats::runif(1, 0, 2 * pi)
#'     return(c(cos(ang), sin(ang)))
#'   }
#'
#'   out <- tryCatch(
#'     as.numeric(Directional::rvmf(1, mu, kappa)),
#'     error = function(e) mu
#'   )
#'   .safe_unit_vec(out, fallback = fallback)
#' }
#'
#' # -------------------------------------------------------------------------
#' # Helper: fallback Gaussian compatibility density on grid
#' # -------------------------------------------------------------------------
#' .fallback_gaussian_density <- function(grid_xy, centroid, scale = 0.08, eps = 1e-12) {
#'   Z <- sweep(as.matrix(grid_xy), 2, as.numeric(centroid), "-", check.margin = FALSE)
#'   dens <- exp(-rowSums(Z^2) / max(scale, eps))
#'   dens[!is.finite(dens)] <- eps
#'   pmax(dens, eps)
#' }
#'
#'
#' # -------------------------------------------------------------------------
#' # Keep names
#' # -------------------------------------------------------------------------
#' .count_labels <- function(labels, doms) {
#'   tab <- tabulate(match(labels, doms), nbins = length(doms))
#'   stats::setNames(as.integer(tab), doms)
#' }
#'
#' # -------------------------------------------------------------------------
#' # Helper: domain centroid sampler
#' # -------------------------------------------------------------------------
#' .sample_domain_centroid <- function(geom_p, fallback = c(0.5, 0.5)) {
#'   if (is.null(geom_p) || is.null(geom_p$placement)) return(as.numeric(fallback))
#'
#'   mu_G <- geom_p$placement$mu_G
#'   Sigma_G <- geom_p$placement$Sigma_G
#'
#'   if (is.null(mu_G) || length(mu_G) != 2L || any(!is.finite(mu_G))) {
#'     mu_G <- as.numeric(fallback)
#'   }
#'   if (is.null(Sigma_G) || !all(dim(as.matrix(Sigma_G)) == c(2L, 2L))) {
#'     Sigma_G <- diag(1e-4, 2L)
#'   }
#'
#'   .safe_mvrnorm_vec(mu = mu_G, Sigma = Sigma_G, jitter = 1e-8)
#' }
#'
#'
#' # -------------------------------------------------------------------------
#' # Helper: draw one component from pooled global_theta
#' #   - pooled prior draw + cloud-anchored hybrid
#' # -------------------------------------------------------------------------
#' sample_one_component_from_global <- function(g,
#'                                              lambda_cloud = 0.7,
#'                                              fallback_phi = c(1, 0),
#'                                              eps = 1e-8) {
#'   # --- 1. anchor picks from empirical cloud ---
#'   c_anchor <- if (!is.null(g$center_cloud) && length(g$center_cloud) > 0L) {
#'     as.numeric(g$center_cloud[[sample.int(length(g$center_cloud), 1L)]])
#'   } else {
#'     as.numeric(g$mu_c)
#'   }
#'
#'   phi_anchor <- if (!is.null(g$orientation_cloud) && length(g$orientation_cloud) > 0L) {
#'     .safe_unit_vec(
#'       g$orientation_cloud[[sample.int(length(g$orientation_cloud), 1L)]],
#'       fallback = fallback_phi
#'     )
#'   } else {
#'     .safe_unit_vec(g$mu_phi, fallback = fallback_phi)
#'   }
#'
#'   # --- 2. prior draws ---
#'   c_prior <- .safe_mvrnorm_vec(g$mu_c, g$Sigma_c, jitter = eps)
#'
#'   r_prior <- abs(stats::rnorm(
#'     1,
#'     mean = as.numeric(g$mu_r),
#'     sd   = sqrt(max(as.numeric(g$var_r), eps))
#'   ))
#'
#'   phi_prior <- .safe_rvmf1(
#'     mu = g$mu_phi,
#'     kappa = as.numeric(g$kappa_phi),
#'     fallback = fallback_phi
#'   )
#'
#'   tau_prior <- stats::rgamma(
#'     1,
#'     shape = max(as.numeric(g$a_tau), eps),
#'     rate  = max(as.numeric(g$b_tau), eps)
#'   )
#'   tau_prior <- max(tau_prior, 1e-6)
#'
#'   # --- 3. hybridize cloud anchor + pooled prior draw ---
#'   c_draw <- lambda_cloud * c_anchor + (1 - lambda_cloud) * c_prior
#'
#'   phi_draw_raw <- lambda_cloud * phi_anchor + (1 - lambda_cloud) * phi_prior
#'   phi_draw <- .safe_unit_vec(phi_draw_raw, fallback = fallback_phi)
#'
#'   list(
#'     c   = as.numeric(c_draw),
#'     r   = as.numeric(r_prior),
#'     phi = as.numeric(phi_draw),
#'     tau = as.numeric(tau_prior)
#'   )
#' }
#'
#' # -------------------------------------------------------------------------
#' # Helper: sample ONE domain-level shared sigma_sq from pooled prior
#' # FG-EM uses a common sigma_sq across components within a domain
#' # -------------------------------------------------------------------------
#' sample_shared_sigma_sq_from_pooled <- function(global_theta, eps = 1e-8) {
#'   a_vec <- vapply(global_theta, function(g) as.numeric(g$a_sigma), numeric(1))
#'   b_vec <- vapply(global_theta, function(g) as.numeric(g$b_sigma), numeric(1))
#'
#'   a_bar <- mean(a_vec[is.finite(a_vec) & a_vec > 0], na.rm = TRUE)
#'   b_bar <- mean(b_vec[is.finite(b_vec) & b_vec > 0], na.rm = TRUE)
#'
#'   if (!is.finite(a_bar) || a_bar <= 0) a_bar <- 10
#'   if (!is.finite(b_bar) || b_bar <= 0) b_bar <- 10
#'
#'   prec <- stats::rgamma(1, shape = a_bar, rate = b_bar)
#'   prec <- max(as.numeric(prec), eps)
#'   1 / prec
#' }
#'
#' # -------------------------------------------------------------------------
#' # Main: sample sample-specific FG mixture from pooled prior
#' # output matches evaluate stage:
#' #   pi, c_hat, r_hat, phi_hat, tau_hat, sigma_sq
#' # -------------------------------------------------------------------------
#' sample_fgkmm_domain_from_pooled <- function(pooled,
#'                                             lambda_cloud = 0.7,
#'                                             eps = 1e-8) {
#'   if (is.null(pooled) || is.null(pooled$global_theta) || length(pooled$global_theta) == 0L) {
#'     stop("pooled must contain non-empty global_theta.")
#'   }
#'
#'   M <- length(pooled$global_theta)
#'
#'   alpha_hat <- pooled$alpha_hat
#'   if (is.null(alpha_hat) || length(alpha_hat) != M || any(!is.finite(alpha_hat))) {
#'     alpha_hat <- rep(1, M)
#'   }
#'   alpha_hat <- pmax(as.numeric(alpha_hat), eps)
#'
#'   pi_draw <- as.numeric(MCMCpack::rdirichlet(1, alpha_hat))
#'   pi_draw <- pi_draw / sum(pi_draw)
#'
#'   comp_list <- lapply(pooled$global_theta, function(g) {
#'     sample_one_component_from_global(
#'       g = g,
#'       lambda_cloud = lambda_cloud,
#'       fallback_phi = c(1, 0),
#'       eps = eps
#'     )
#'   })
#'
#'   c_hat <- do.call(rbind, lapply(comp_list, `[[`, "c"))
#'   phi_hat <- do.call(rbind, lapply(comp_list, `[[`, "phi"))
#'   r_hat <- vapply(comp_list, `[[`, numeric(1), "r")
#'   tau_hat <- vapply(comp_list, `[[`, numeric(1), "tau")
#'
#'   sigma_sq <- sample_shared_sigma_sq_from_pooled(pooled$global_theta, eps = eps)
#'
#'   list(
#'     pi = pi_draw,
#'     c_hat = c_hat,
#'     r_hat = r_hat,
#'     phi_hat = phi_hat,
#'     tau_hat = tau_hat,
#'     sigma_sq = sigma_sq,
#'     component_list = comp_list
#'   )
#' }
#'
#' # -------------------------------------------------------------------------
#' # Evaluate q_d(s) from one sample-specific FG mixture draw
#' #   q_d(s) = sum_m pi_m f_FG(s - centroid | theta_m)
#' # query_xy: absolute coordinates (n x 2)
#' # centroid: absolute domain centroid (length 2)
#' # fg_draw: output of sample_fgkmm_domain_from_pooled()
#' # -------------------------------------------------------------------------
#' evaluate_fg_mixture_density <- function(query_xy,
#'                                         centroid,
#'                                         fg_draw,
#'                                         eps = 1e-300) {
#'   X <- as.matrix(query_xy)
#'   n <- nrow(X)
#'   if (n == 0L) return(numeric(0))
#'
#'   centroid <- as.numeric(centroid)
#'   if (length(centroid) != ncol(X)) stop("centroid dimension mismatch.")
#'
#'   X_centered <- sweep(X, 2, centroid, "-", check.margin = FALSE)
#'   M <- length(fg_draw$pi)
#'
#'   ll_mat <- matrix(NA_real_, nrow = n, ncol = M)
#'
#'   for (m in seq_len(M)) {
#'     ll_mat[, m] <-
#'       log(pmax(fg_draw$pi[m], eps)) +
#'       log_fg_one_comp(
#'         X = X_centered,
#'         c_k = fg_draw$c_hat[m, ],
#'         r_k = fg_draw$r_hat[m],
#'         sigma_sq = fg_draw$sigma_sq,
#'         phi_k = fg_draw$phi_hat[m, ],
#'         tau_k = fg_draw$tau_hat[m]
#'       )
#'   }
#'
#'   dens <- exp(matrixStats::rowLogSumExps(ll_mat))
#'   dens[!is.finite(dens)] <- eps
#'   pmax(dens, eps)
#' }
#'
#' # -------------------------------------------------------------------------
#' # Wrapper:
#' #   if apply_geo = TRUE,
#' #     q_geo(s) = kappa^{-1} q(T^{-1}(s))
#' #   else
#' #     q(s)
#' #
#' # If fg_draw is NULL, it samples one sample-specific mixture once.
#' # -------------------------------------------------------------------------
#' evaluate_pooled_fg_density <- function(query_xy,
#'                                        centroid,
#'                                        pooled,
#'                                        apply_geo = FALSE,
#'                                        kappa_geo = 1,
#'                                        fg_draw = NULL,
#'                                        lambda_cloud = 0.7,
#'                                        density_eps = 1e-300,
#'                                        param_eps = 1e-8,
#'                                        return_draw = FALSE) {
#'   if (is.null(fg_draw)) {
#'     fg_draw <- sample_fgkmm_domain_from_pooled(
#'       pooled = pooled,
#'       lambda_cloud = lambda_cloud,
#'       eps = param_eps
#'     )
#'   }
#'
#'   Xq <- as.matrix(query_xy)
#'   centroid <- as.numeric(centroid)
#'
#'   if (apply_geo) {
#'     kg <- as.numeric(kappa_geo)
#'     if (!is.finite(kg) || kg <= 0) kg <- 1
#'
#'     X_eval <- .inverse_warp_query(
#'       query_xy = Xq,
#'       centroid = centroid,
#'       kappa_geo = kg
#'     )
#'     jac <- 1 / kg
#'   } else {
#'     X_eval <- Xq
#'     jac <- 1
#'   }
#'
#'   dens0 <- evaluate_fg_mixture_density(
#'     query_xy = X_eval,
#'     centroid = centroid,
#'     fg_draw = fg_draw,
#'     eps = density_eps
#'   )
#'
#'   dens <- jac * dens0
#'   dens[!is.finite(dens)] <- density_eps
#'   dens <- pmax(dens, density_eps)
#'
#'   if (return_draw) {
#'     return(list(
#'       density = dens,
#'       fg_draw = fg_draw,
#'       query_eval = X_eval,
#'       jacobian = jac
#'     ))
#'   }
#'
#'   dens
#' }
#'
#'
#'
#'
#' ####simulate helpers###
#' # ==============================================================================
#' # File: simulatespaCraft.R
#' # Replacement fixed-grid exact-count relabeling implementation
#' # ==============================================================================
#' .normalize_coords_unit <- function(coords_df) {
#'   out <- coords_df
#'   xr <- max(out$x, na.rm = TRUE) - min(out$x, na.rm = TRUE)
#'   yr <- max(out$y, na.rm = TRUE) - min(out$y, na.rm = TRUE)
#'   if (!is.finite(xr) || xr <= 0) xr <- 1
#'   if (!is.finite(yr) || yr <= 0) yr <- 1
#'   out$x <- (out$x - min(out$x, na.rm = TRUE)) / xr
#'   out$y <- (out$y - min(out$y, na.rm = TRUE)) / yr
#'   out
#' }
#'
#' .pick_fixed_grid_template <- function(pilot_data, grp_id = NULL) {
#'   if (length(pilot_data) == 0L) stop("pilot_data is empty.")
#'
#'   idx_pool <- seq_along(pilot_data)
#'   if (!is.null(grp_id)) {
#'     idx_grp <- which(vapply(
#'       pilot_data,
#'       function(x) as.character(x$group) == as.character(grp_id),
#'       logical(1)
#'     ))
#'     if (length(idx_grp) > 0L) idx_pool <- idx_grp
#'   }
#'
#'   pick <- sample(idx_pool, size = 1L)
#'   coords <- .normalize_coords_unit(pilot_data[[pick]]$coords)
#'
#'   list(
#'     template_index = pick,
#'     grid_df = data.frame(
#'       x = coords$x,
#'       y = coords$y,
#'       stringsAsFactors = FALSE
#'     )
#'   )
#' }
#'
#' .extract_geometry_group <- function(params_geometry, grp_chr) {
#'   if (is.null(params_geometry) || length(params_geometry) == 0L) return(NULL)
#'
#'   # 1. exact group key
#'   if (!is.null(params_geometry[[grp_chr]])) return(params_geometry[[grp_chr]])
#'
#'   # 2. try the other numeric string key ("0" or "1")
#'   alt <- if (grp_chr == "0") "1" else "0"
#'   if (!is.null(params_geometry[[alt]])) {
#'     warning(sprintf(
#'       ".extract_geometry_group: key '%s' not found; falling back to '%s'.",
#'       grp_chr, alt
#'     ))
#'     return(params_geometry[[alt]])
#'   }
#'
#'   # 3. last resort: first element
#'   params_geometry[[1]]
#' }
#'
#'
#'
#'
#' # ------------------------------------------------------------------------------
#' # Whole-tissue spatial cache
#' #   spatial_mode = "original" : whole-tissue KNN graph
#' #   spatial_mode = "basis"    : whole-tissue low-rank basis knots
#' # ------------------------------------------------------------------------------
#' build_spatial_cache <- function(sim_coords_df,
#'                                 spatial_mode = c("original", "basis"),
#'                                 graph_k = 20L,
#'                                 basis_rank = 25L,
#'                                 basis_knot_method = c("sample", "grid"),
#'                                 basis_seed = 1L) {
#'   spatial_mode <- match.arg(spatial_mode)
#'   basis_knot_method <- match.arg(basis_knot_method)
#'
#'   if (!all(c("x", "y", "domain") %in% names(sim_coords_df))) {
#'     stop("sim_coords_df must contain x, y, domain.")
#'   }
#'
#'   coords_mat <- as.matrix(sim_coords_df[, c("x", "y"), drop = FALSE])
#'   n_cells <- nrow(coords_mat)
#'   sim_domains <- unique(as.character(sim_coords_df$domain))
#'
#'   out <- list(
#'     spatial_mode = spatial_mode,
#'     coords_mat = coords_mat,
#'     n_cells = n_cells,
#'     sim_domains = sim_domains
#'   )
#'
#'   # whole-tissue KNN graph cache
#'   if (spatial_mode == "original") {
#'     if (!requireNamespace("FNN", quietly = TRUE)) stop("Need FNN for spatial_mode='original'.")
#'
#'     k_use <- min(as.integer(graph_k), max(1L, n_cells - 1L))
#'
#'     if (n_cells <= 1L) {
#'       out$graph_k <- 0L
#'       out$nn_index <- matrix(integer(0), nrow = n_cells, ncol = 0L)
#'       out$nn_dist  <- matrix(numeric(0), nrow = n_cells, ncol = 0L)
#'     } else {
#'       nn <- FNN::get.knn(coords_mat, k = k_use)
#'       out$graph_k  <- k_use
#'       out$nn_index <- nn$nn.index
#'       out$nn_dist  <- nn$nn.dist
#'     }
#'   }
#'
#'   # low-rank basis cache
#'   if (spatial_mode == "basis") {
#'     basis_rank <- min(as.integer(basis_rank), n_cells)
#'     basis_rank <- max(1L, basis_rank)
#'
#'     if (basis_knot_method == "sample") {
#'       ord <- order(coords_mat[, 1], coords_mat[, 2])
#'       pick <- unique(round(seq(1, n_cells, length.out = basis_rank)))
#'       knots <- coords_mat[ord[pick], , drop = FALSE]
#'     } else {
#'       gx <- ceiling(sqrt(basis_rank))
#'       gy <- ceiling(basis_rank / gx)
#'
#'       x_seq <- seq(min(coords_mat[, 1]), max(coords_mat[, 1]), length.out = gx)
#'       y_seq <- seq(min(coords_mat[, 2]), max(coords_mat[, 2]), length.out = gy)
#'       grid_df <- expand.grid(x = x_seq, y = y_seq)
#'       knots <- as.matrix(grid_df[seq_len(min(nrow(grid_df), basis_rank)), c("x", "y"), drop = FALSE])
#'     }
#'
#'     out$basis_rank <- nrow(knots)
#'     out$basis_knots <- knots
#'     out$basis_knot_method <- basis_knot_method
#'     out$basis_seed <- basis_seed
#'   }
#'
#'   out
#' }
#'
#'
