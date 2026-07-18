# ||row_i||^2 for the rows of a matrix.
#' @noRd
row_norms2 <- function(A) rowSums(A * A)


# Log normalizing constant of the von Mises-Fisher distribution (exponentially
# scaled Bessel for numerical stability; the +tau scale is added back by callers).
#' @noRd
log_Cd <- function(tau, d) {
  if (tau < 0) stop("tau must be non-negative")
  if (tau == 0) return(-log(2) - (d/2) * log(pi) + lgamma(d/2))
  res <- (d/2 - 1) * log(tau) -
    log(besselI(tau, nu = (d/2 - 1), expon.scaled = TRUE)) - (d/2) * log(2 * pi)
  return(res)
}


# Row-wise log density of one Fisher-Gaussian component, log f_k(x_i | theta_k);
# returns a numeric vector of length nrow(X).
#' @noRd
log_fg_one_comp <- function(X, c_k, r_k, sigma_sq, phi_k, tau_k) {
  n <- nrow(X); d <- ncol(X)
  if (sigma_sq <= 0) stop("sigma_sq must be positive")
  if (tau_k <= 0) stop("tau_k must be positive")

  # z_i = x_i - c_k
  Z <- sweep(X, 2, c_k, FUN = "-")
  z_norm2 <- row_norms2(Z)

  # v_i = tau_k * phi_k + (r_k/sigma_sq) * z_i
  V <- Z * (r_k / sigma_sq)
  V <- sweep(V, 2, tau_k * phi_k, FUN = "+")
  v_norm <- sqrt(row_norms2(V))

  # log C_d(.) is exponentially scaled, so add back +v_norm terms below
  logCd_tau  <- log_Cd(tau_k, d)
  logCd_v    <- vapply(v_norm, log_Cd, numeric(1), d = d)

  const <- logCd_tau - (d * log(2 * pi * sigma_sq) / 2) - (r_k^2) / (2 * sigma_sq) - tau_k
  const - logCd_v - (z_norm2 / (2 * sigma_sq)) + v_norm
}


# Incomplete (observed-data) log-likelihood of the FG mixture.
#' @noRd
log_likelihood <- function(X, pi_hat, c_hat, r_hat, phi_hat, tau_hat, sigma_sq) {
  n <- nrow(X); M <- length(pi_hat)
  ll_mat <- matrix(NA_real_, n, M)

  for (k in 1:M) {
    ll_mat[, k] <- log(pi_hat[k]) +
      log_fg_one_comp(X, c_hat[k, ], r_hat[k], sigma_sq, phi_hat[k, ], tau_hat[k])
  }
  # log-sum-exp across components per row, then sum over rows
  sum(matrixStats::rowLogSumExps(ll_mat))
}


# k-means initialization of the FG-EM parameters (centers, radii, directions,
# concentrations, shared variance).
#' @noRd
initial_kmeans <- function(X, M) {
  # ---- 0) Hard numeric coercion ----
  if (!is.matrix(X)) X <- as.matrix(X)
  storage.mode(X) <- "double"
  if (!all(is.finite(X))) {
    X <- X[apply(X, 1, function(r) all(is.finite(r))), , drop = FALSE]
  }

  n <- nrow(X); d <- ncol(X)
  if (n < M) stop(sprintf("Need at least M=%d points for initialization; got n=%d", M, n))

  pi_hat <- rep(1/M, M)

  # ---- 1) k-means init (base kmeans for portability) ----
  # Fixed seed for reproducible initialization; restore the caller's RNG on exit.
  if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
    .old_seed <- get(".Random.seed", envir = globalenv(), inherits = FALSE)
    on.exit(assign(".Random.seed", .old_seed, envir = globalenv()), add = TRUE)
  }
  set.seed(1L)
  km <- stats::kmeans(X, centers = M, nstart = 10, iter.max = 100, algorithm = "Lloyd")
  c_hat <- km$centers
  cl    <- km$cluster

  storage.mode(c_hat) <- "double"

  # ---- 2) r_hat: mean distance to center per cluster ----
  r_hat <- numeric(M)
  for (k in seq_len(M)) {
    idx <- which(cl == k)
    if (length(idx) == 0L) { r_hat[k] <- 0; next }
    Zk <- sweep(X[idx, , drop = FALSE], 2, c_hat[k, ], "-")
    storage.mode(Zk) <- "double"
    r_hat[k] <- mean(sqrt(rowSums(Zk * Zk)))
  }

  # ---- 3) sigma_sq from squared residuals ----
  Z <- X - c_hat[cl, , drop = FALSE]
  storage.mode(Z) <- "double"
  sigma_sq <- mean(rowSums(Z * Z))

  # ---- 4) Recover Y = unit direction vectors ----
  eps <- 1e-12
  r_by_row <- pmax(r_hat[cl], eps)
  Y <- Z / r_by_row
  storage.mode(Y) <- "double"

  rn <- sqrt(rowSums(Y * Y))
  rn[rn < eps] <- 1.0
  Y <- Y / rn
  storage.mode(Y) <- "double"

  # ---- 5) phi_hat and tau_hat per cluster (fully guarded) ----
  phi_hat <- matrix(0, M, d)
  storage.mode(phi_hat) <- "double"
  tau_hat <- numeric(M)

  for (k in seq_len(M)) {
    idx <- which(cl == k)

    if (length(idx) <= 1L) {
      phi_hat[k, ] <- rep(1 / sqrt(d), d)
      storage.mode(phi_hat[k, ]) <- "double"
      tau_hat[k] <- 1
      next
    }

    ybar <- colSums(Y[idx, , drop = FALSE])
    storage.mode(ybar) <- "double"
    nrm <- sqrt(sum(ybar * ybar))

    if (!is.finite(nrm) || nrm < 1e-12) {
      phi_hat[k, ] <- rep(1 / sqrt(d), d)
    } else {
      phi_hat[k, ] <- ybar / nrm
    }
    storage.mode(phi_hat[k, ]) <- "double"

    # robust tau via angular spread
    Yk <- Y[idx, , drop = FALSE]
    storage.mode(Yk) <- "double"

    phi_col <- matrix(phi_hat[k, ], ncol = 1)
    storage.mode(phi_col) <- "double"

    dots <- Yk %*% phi_col
    dots <- as.numeric(dots)
    dots[!is.finite(dots)] <- 0
    dots <- pmin(pmax(dots, -1), 1)

    ang <- acos(dots)
    ang <- ang[is.finite(ang)]
    if (length(ang) == 0L) {
      tau_hat[k] <- 1
    } else {
      med_ang <- stats::median(ang)
      if (!is.finite(med_ang) || med_ang <= 0) med_ang <- 1
      tau_hat[k] <- max(1 / med_ang, 1e-3)
    }
  }

  list(
    pi_hat   = pi_hat,
    c_hat    = c_hat,
    r_hat    = r_hat,
    sigma_sq = sigma_sq,
    phi_hat  = phi_hat,
    tau_hat  = tau_hat
  )
}


# Fit the Fisher-Gaussian kernel mixture by EM. The iteration runs in the
# compiled engine FG_EM_cpp; this wrapper only builds the initialization.
#' @useDynLib spaCraft, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @noRd
FG_EM <- function(x, M = 5, iter_max = 1000, tol = 1e-2, initial = NULL, n_cores = 1) {
  x <- as.matrix(x)
  storage.mode(x) <- "double"

  if (is.null(initial)) {
    initialPars <- initial_kmeans(X = x, M = M)
  } else {
    initialPars <- initial
  }

  res <- FG_EM_cpp(
    X = x,
    initial = initialPars,
    iter_max = as.integer(iter_max),
    tol = as.double(tol),
    n_cores = as.integer(n_cores)
  )

  return(res)
}


# Fit one FG mixture of order M and attach AIC / BIC.
#' @noRd
FG_EM_with_criteria <- function(x, M, iter_max, tol) {
  fit <- FG_EM(x, M, iter_max, tol)

  # Number of free parameters:
  # pi_k (M-1), c_k (M*d), r_k (M), phi_k (M*(d-1) under unit norm), tau_k (M), sigma^2 (1)
  n <- nrow(x)
  d <- ncol(x)
  p <- (M - 1) + M * d + M + M * (d - 1) + M + 1

  final_loglik <- utils::tail(fit$log_likelihood, 1)

  BIC <- -2 * final_loglik + p * log(n)
  AIC <- -2 * final_loglik + 2 * p

  return(list(
    model = fit,
    AIC = AIC,
    BIC = BIC,
    n_parameters = p,
    converged_iter = fit$n_iter
  ))
}


# Select the mixture order M by AIC/BIC over a candidate grid.
#' @noRd
select_best_M <- function(x, M_candidates = 2:5, iter_max = 1000, tol = 1e-2) {
  results <- list()
  criteria_df <- data.frame(
    M = M_candidates,
    AIC = NA_real_,
    BIC = NA_real_,
    converged = NA,
    n_iter = NA
  )

  for (i in seq_along(M_candidates)) {
    M <- M_candidates[i]
    message("\nFitting M = ", M)
    result <- FG_EM_with_criteria(x, M, iter_max, tol)

    results[[paste0("M_", M)]] <- result

    criteria_df$AIC[i] <- result$AIC
    criteria_df$BIC[i] <- result$BIC
    criteria_df$converged[i] <- (result$model$n_iter < iter_max)
    criteria_df$n_iter[i] <- result$model$n_iter
  }

  best_AIC_idx <- which.min(criteria_df$AIC)
  best_BIC_idx <- which.min(criteria_df$BIC)

  return(list(
    all_models = results,
    criteria_table = criteria_df,
    best_M_AIC = M_candidates[best_AIC_idx],
    best_M_BIC = M_candidates[best_BIC_idx],
    best_model_AIC = results[[best_AIC_idx]]$model,
    best_model_BIC = results[[best_BIC_idx]]$model
  ))
}


# Pool per-domain (or per-sample) FGKMM fits into a global prior. Supports both
# the multi-sample (K > 1) empirical pooling and the single-sample (K = 1)
# shrinkage pooling (when p_B and gamma are supplied).
#' @noRd
pool_FGKMM_parameters <- function(theta_list, omega_list, eps = 1e-6,
                                  p_B = NULL, gamma = NULL) {

  is_domain_pooling <- !is.null(p_B) && !is.null(gamma)

  #== Case 1: Single-sample (K = 1) shrinkage pooling =========================
  if (is_domain_pooling) {

    successful_indices <- !sapply(theta_list, is.null)
    theta_list_d <- theta_list[successful_indices]
    omega_list_d <- omega_list[successful_indices]

    D <- length(theta_list_d)
    if (D < 1) stop("K=1 case requires at least 1 valid domain fit.")
    M <- length(theta_list_d[[1]])

    omega_mat <- do.call(rbind, omega_list_d)
    omega_bar <- colMeans(omega_mat)
    omega_var <- pmax(matrixStats::colVars(omega_mat), eps)
    alpha0_hat <- mean((omega_bar * (1 - omega_bar)) / omega_var - 1, na.rm = TRUE)
    alpha0_hat <- max(alpha0_hat, eps)
    alpha_hat <- alpha0_hat * omega_bar
    alpha_hat <- pmax(alpha_hat, eps)

    global_theta <- lapply(1:M, function(m) {
      c_list     <- lapply(theta_list_d, function(s) s[[m]]$c)
      r_list     <- sapply(theta_list_d, function(s) s[[m]]$r)
      tau_list   <- sapply(theta_list_d, function(s) s[[m]]$tau)
      sigma2_list<- sapply(theta_list_d, function(s) s[[m]]$sigma2)
      phi_list   <- lapply(theta_list_d, function(s) s[[m]]$phi)

      # --- Center ---
      mu_c <- colMeans(do.call(rbind, c_list))
      emp_cov_c <- stats::cov(do.call(rbind, c_list))
      if (anyNA(emp_cov_c)) emp_cov_c <- diag(eps, length(mu_c))
      Sigma_c_prior <- diag((p_B * abs(mu_c))^2)
      Sigma_c <- gamma * emp_cov_c + (1 - gamma) * Sigma_c_prior

      # --- Orientation ---
      phi_mat <- do.call(rbind, phi_list)
      R_vec <- colMeans(phi_mat)
      R_bar <- sqrt(sum(R_vec^2))
      mu_phi <- if (R_bar > 0) R_vec / R_bar else rep(0, ncol(phi_mat))

      kappa_emp   <- estimate_kappa(R_bar)
      kappa_prior <- 2 * (1 - p_B)
      kappa_phi   <- gamma * kappa_emp + (1 - gamma) * kappa_prior

      # --- Other parameters ---
      mu_r <- mean(r_list)
      mu_tau <- mean(tau_list)
      mean_inv_sigma <- mean(1 / sigma2_list)

      var_r <- (p_B * abs(mu_r))^2
      var_tau <- (p_B * abs(mu_tau))^2
      var_inv_sigma <- (p_B * abs(mean_inv_sigma))^2

      a_tau <- mu_tau^2 / var_tau; b_tau <- mu_tau / var_tau
      a_sigma <- mean_inv_sigma^2 / var_inv_sigma
      b_sigma <- mean_inv_sigma / var_inv_sigma

      list(
        mu_c = mu_c, Sigma_c = Sigma_c,
        mu_r = mu_r, var_r = var_r,
        mu_phi = mu_phi, kappa_phi = kappa_phi,
        a_tau = a_tau, b_tau = b_tau,
        a_sigma = a_sigma, b_sigma = b_sigma,
        center_cloud = c_list,
        orientation_cloud = phi_list
      )
    })

    return(list(global_theta = global_theta, alpha_hat = alpha_hat))

    #== Case 2: Multi-sample (K > 1) empirical pooling ========================
  } else {

    successful_indices <- !sapply(theta_list, is.null)
    theta_list_clean <- theta_list[successful_indices]
    omega_list_clean <- omega_list[successful_indices]

    K <- length(theta_list_clean)
    if (K <= 1) stop("Pooling by samples requires K > 1, or supply p_B and gamma for K=1.")
    M <- length(theta_list_clean[[1]])

    omega_mat <- do.call(rbind, omega_list_clean)
    omega_bar <- colMeans(omega_mat)
    omega_var <- pmax(matrixStats::colVars(omega_mat), eps)
    alpha0_hat <- mean((omega_bar * (1 - omega_bar)) / omega_var - 1, na.rm = TRUE)
    alpha0_hat <- max(alpha0_hat, eps)
    alpha_hat <- alpha0_hat * omega_bar
    alpha_hat <- pmax(alpha_hat, eps)

    # Global parameters per FG component (pooled over pilot samples K).
    global_theta <- lapply(1:M, function(m) {
      c_mat <- do.call(rbind, lapply(theta_list_clean, function(s) s[[m]]$c))
      phi_mat <- do.call(rbind, lapply(theta_list_clean, function(s) s[[m]]$phi))
      r_vec <- sapply(theta_list_clean, function(s) s[[m]]$r)
      sigma2_vec <- sapply(theta_list_clean, function(s) s[[m]]$sigma2)
      tau_vec <- sapply(theta_list_clean, function(s) s[[m]]$tau)

      R_vec <- colMeans(phi_mat)
      R_bar <- sqrt(sum(R_vec^2))
      mu_phi <- if (R_bar > 0) R_vec / R_bar else rep(0, ncol(phi_mat))
      kappa_phi <- estimate_kappa(R_bar)

      mean_tau <- mean(tau_vec); var_tau <- pmax(stats::var(tau_vec, na.rm = TRUE), eps)
      a_tau <- mean_tau^2 / var_tau; b_tau <- mean_tau / var_tau

      inv_sigma_vec <- 1 / sigma2_vec
      mean_inv_sigma <- mean(inv_sigma_vec)
      var_inv_sigma <- pmax(stats::var(inv_sigma_vec, na.rm = TRUE), eps)
      a_sigma <- mean_inv_sigma^2 / var_inv_sigma; b_sigma <- mean_inv_sigma / var_inv_sigma

      list(
        mu_c = colMeans(c_mat), Sigma_c = stats::cov(c_mat),
        mu_r = mean(r_vec), var_r = pmax(stats::var(r_vec, na.rm = TRUE), eps),
        mu_phi = mu_phi, kappa_phi = kappa_phi,
        a_tau = a_tau, b_tau = b_tau,
        a_sigma = a_sigma, b_sigma = b_sigma,
        center_cloud = lapply(theta_list_clean, function(s) s[[m]]$c),
        orientation_cloud = lapply(theta_list_clean, function(s) s[[m]]$phi)
      )
    })

    return(list(global_theta = global_theta, alpha_hat = alpha_hat))
  }
}


# Invert the mean resultant length R_bar to the vMF concentration kappa.
#' @noRd
estimate_kappa <- function(R_bar, d = 2) {
  if (R_bar < 1e-6) return(1e-3)  # Avoid instability

  if (d == 2) {
    kappa <- (2 * R_bar - R_bar^3) / (1 - R_bar^2)
  } else if (d == 3) {
    kappa <- R_bar * (d - R_bar^2) / (1 - R_bar^2)
  } else {
    kappa <- R_bar * (d - 1 - R_bar^2) / (1 - R_bar^2)
  }

  return(max(kappa, 1e-3))
}


# Per-component quantile bounds for (r, sigma, tau) from the pilot FGKMM fits,
# used to bound the posterior shape sampler.
#' @noRd
get_r_sigma_bounds <- function(theta_list, M, lower_q = 0.25, upper_q = 0.75) {

  # 1. Keep only fits that match the expected M
  valid_thetas <- list()
  if (!is.null(theta_list)) {
    for (i in seq_along(theta_list)) {
      if (!is.null(theta_list[[i]]) && length(theta_list[[i]]) == M) {
        valid_thetas[[length(valid_thetas) + 1L]] <- theta_list[[i]]
      }
    }
  }

  K <- length(valid_thetas)

  # 2. No valid pilot data -> broad default bounds
  if (K == 0L) {
    return(list(
      r_min     = rep(0, M),      r_max     = rep(10, M),
      sigma_min = rep(0.001, M),  sigma_max = rep(1.0, M),
      tau_min   = rep(0.1, M),    tau_max   = rep(500, M)
    ))
  }

  # 3. Assemble per-(sample, component) matrices
  r_mat     <- matrix(NA_real_, nrow = K, ncol = M)
  sigma_mat <- matrix(NA_real_, nrow = K, ncol = M)
  tau_mat   <- matrix(NA_real_, nrow = K, ncol = M)

  for (k in seq_len(K)) {
    for (m in seq_len(M)) {
      r_mat[k, m] <- valid_thetas[[k]][[m]]$r

      # Handle sigma_sq vs sigma2 naming differences robustly
      s2 <- if (!is.null(valid_thetas[[k]][[m]]$sigma_sq)) {
        valid_thetas[[k]][[m]]$sigma_sq
      } else {
        valid_thetas[[k]][[m]]$sigma2
      }
      sigma_mat[k, m] <- if (!is.null(s2)) sqrt(s2) else 0.1

      tau_mat[k, m] <- valid_thetas[[k]][[m]]$tau
    }
  }

  # 4. Quantile bounds (ignoring NAs)
  other_lower_q <- max(0.05, lower_q)
  other_upper_q <- min(upper_q, 0.95)

  list(
    r_min     = apply(r_mat, 2, stats::quantile, probs = other_lower_q, na.rm = TRUE),
    r_max     = apply(r_mat, 2, stats::quantile, probs = other_upper_q, na.rm = TRUE),
    sigma_min = apply(sigma_mat, 2, stats::quantile, probs = lower_q, na.rm = TRUE),
    sigma_max = apply(sigma_mat, 2, stats::quantile, probs = upper_q, na.rm = TRUE),
    tau_min   = apply(tau_mat, 2, stats::quantile, probs = lower_q, na.rm = TRUE),
    tau_max   = apply(tau_mat, 2, stats::quantile, probs = upper_q, na.rm = TRUE)
  )
}


# Draw N synthetic spot coordinates (centered at the origin) from the pooled
# FG posterior: a hybrid of cloud anchors and the pooled prior. The caller adds
# the absolute domain centroid (see simulatespaCraft).
#' @noRd
simulate_from_FG_post <- function(N,
                                  global_theta,
                                  alpha_hat,
                                  r_min,
                                  r_max,
                                  sigma_min,
                                  sigma_max,
                                  tau_min,
                                  tau_max,
                                  lambda = 0.7) {

  M      <- length(global_theta)
  Pi     <- MCMCpack::rdirichlet(1, alpha_hat)[1, ]
  labels <- sample(1:M, size = N, replace = TRUE, prob = Pi)
  coords <- matrix(NA_real_, nrow = N, ncol = 2)

  for (m in 1:M) {
    idx    <- which(labels == m)
    if (length(idx) == 0L) next

    g          <- global_theta[[m]]
    cloud_len  <- length(g$center_cloud)
    n_gen      <- length(idx)

    # cloud index sampling (j ~ Uniform{1,...,K})
    k_idx <- if (n_gen <= cloud_len) {
      sample(cloud_len, n_gen, replace = FALSE)
    } else {
      sample(cloud_len, n_gen, replace = TRUE)
    }

    c_cloud_mat   <- do.call(rbind, g$center_cloud[k_idx])
    phi_cloud_mat <- do.call(rbind, g$orientation_cloud[k_idx])
    mu_c          <- g$mu_c
    mu_phi        <- g$mu_phi

    # hybrid center : c_new = lambda * c_cloud + (1 - lambda) * mu_c
    # hybrid phi    : phi_new = normalize(lambda * phi_cloud + (1 - lambda) * mu_phi)
    c_k_mat  <- matrix(NA_real_, nrow = n_gen, ncol = 2)
    phi_mat  <- matrix(NA_real_, nrow = n_gen, ncol = 2)
    for (j in seq_len(n_gen)) {
      c_k_mat[j, ]  <- lambda * c_cloud_mat[j, ] + (1 - lambda) * mu_c
      phi_j         <- lambda * phi_cloud_mat[j, ] + (1 - lambda) * mu_phi
      phi_norm      <- sqrt(sum(phi_j^2))
      phi_mat[j, ]  <- if (phi_norm > 1e-12) phi_j / phi_norm else mu_phi
    }

    # scalar parameters: moment-based + quantile-bounded sampling
    r_k     <- pmin(r_max[m],     pmax(r_min[m],
                                       abs(stats::rnorm(n_gen, g$mu_r, sqrt(g$var_r)))))
    tau_k   <- pmin(tau_max[m],   pmax(tau_min[m],
                                       stats::rgamma(n_gen, g$a_tau, g$b_tau)))
    sigma_k <- pmin(sigma_max[m], pmax(sigma_min[m],
                                       1 / sqrt(stats::rgamma(n_gen, g$a_sigma, g$b_sigma))))

    # spot coordinates: y_j ~ vMF(phi_j, tau_j), s_j = c_j + r_j * y_j + noise
    for (j in seq_along(idx)) {
      y_j    <- Directional::rvmf(1, phi_mat[j, ], tau_k[j])
      z_j    <- c_k_mat[j, ] + r_k[j] * y_j
      noise  <- stats::rnorm(2, 0, sigma_k[j])
      coords[idx[j], ] <- z_j + noise
    }
  }

  data.frame(x = coords[, 1], y = coords[, 2])
}
