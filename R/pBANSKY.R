#' @title Pilot-Guided Spatial Clustering (pBANKSY, Optimized)
#' @description
#' Pilot-centroid-guided spatial clustering on a robust three-block BANKSY
#' embedding (own expression \eqn{C}, neighbourhood mean \eqn{M}, and
#' neighbourhood gradient \eqn{G}). Latent truth may be supplied in
#' \code{coords$tilde_d} or \code{coords$domain}; recovered labels are written
#' to \code{coords$hat_d}.
#'
#' The recovered labels are statistically equivalent to a reference loop-based
#' implementation (up to floating-point reassociation and PCA sign conventions,
#' neither of which affects the labels), but substantially faster. The main
#' algorithmic gains include:
#' \enumerate{
#'   \item The per-spot neighbour-mean (\eqn{M}) and neighbour-gradient
#'     (\eqn{G}) loops become sparse products with the row-normalized
#'     neighbour-weight matrix \eqn{W}, via
#'     \eqn{M = C W}{M = C W} and
#'     \eqn{G^{2} = (C \odot C)\,W - 2\,C \odot (C W) + C \odot C}{G^2 = (C^2) W - 2 C (C W) + C^2}
#'     (then an elementwise square root and re-scaling).
#'   \item PCA uses a truncated SVD (\code{RSpectra::svds}) for the top
#'     \code{n_pcs} components, falling back to \code{stats::prcomp} when
#'     \pkg{RSpectra} is unavailable.
#'   \item The zero-variance row filter is fully vectorized.
#'   \item Samples are processed in parallel across \code{n_cores}; the shared
#'     pilot embedding is read-only.
#' }
#'
#' @param sim_data List of simulated samples (each with \code{counts} or
#'   \code{logcounts} and a \code{coords} data frame); typically the output of
#'   \code{\link{simulatespaCraft}}.
#' @param pilot_data List of pilot samples. The sample with the most domains
#'   (then the most spots) seeds the cluster centroids.
#' @param lambda Numeric in \eqn{[0,1]}. BANKSY spatial-mixing weight
#'   (0 = expression only, 1 = neighbourhood only).
#' @param k_neighbors Integer. Number of spatial neighbours in the BANKSY graph.
#'   Default 15.
#' @param verbose Logical. Print progress. Default \code{TRUE}.
#' @param use_G Logical. Include the neighbourhood-gradient (\eqn{G}) block.
#'   Default \code{TRUE}.
#' @param balance_mode Character. How the spatial budget is split between the
#'   \eqn{M} and \eqn{G} blocks: one of \code{"energy_proportional"} (default),
#'   \code{"equal_split"}, or \code{"paper_mu"}.
#' @param mu Numeric. Gradient down-weighting, used only when
#'   \code{balance_mode = "paper_mu"}. Default 1.5.
#' @param scale_pilot_coords Logical. Min-max scale the pilot coordinates to the
#'   unit square before graph construction. Default \code{TRUE}.
#' @param scale_sim_coords Logical. Min-max scale the simulated coordinates.
#'   Default \code{TRUE}.
#' @param hard_guidance Logical. If \code{TRUE}, assign each spot to its nearest
#'   pilot centroid instead of running k-means. Default \code{FALSE}.
#' @param n_pcs Integer. Number of principal components of the joint embedding.
#'   Default 20.
#' @param do_hungarian Logical. If \code{TRUE} and latent truth is available,
#'   align recovered labels to truth via the Hungarian algorithm. Default
#'   \code{FALSE}.
#' @param n_cores Integer. Samples processed in parallel (clamped to 1 on
#'   Windows). Default 1.
#'
#' @return A list of processed samples. Each element carries the input
#'   \code{counts}, a \code{coords} data frame with the recovered domain in
#'   \code{hat_d} (and a per-spot \code{is_correct} flag when truth is
#'   available), plus clustering diagnostics (\code{global_accuracy},
#'   \code{n_domains}).
#'
#' @importFrom stats prcomp kmeans setNames dist
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @export
pBANSKY <- function(sim_data,
                    pilot_data,
                    lambda,
                    k_neighbors = 15,
                    verbose = TRUE,
                    use_G = TRUE,
                    balance_mode = c("energy_proportional", "equal_split", "paper_mu"),
                    mu = 1.5,
                    scale_pilot_coords = TRUE,
                    scale_sim_coords = TRUE,
                    hard_guidance = FALSE,
                    n_pcs = 20,
                    do_hungarian = FALSE,
                    n_cores = 1L) {
  balance_mode <- match.arg(balance_mode)

  if (!requireNamespace("FNN", quietly = TRUE))    stop("Package 'FNN' required.")
  if (!requireNamespace("Matrix", quietly = TRUE)) stop("Package 'Matrix' required.")
  if (do_hungarian && !requireNamespace("clue", quietly = TRUE))
    stop("Package 'clue' required when do_hungarian = TRUE.")

  if (!is.numeric(lambda) || length(lambda) != 1L || !is.finite(lambda) || lambda < 0 || lambda > 1)
    stop("`lambda` must be a single numeric value in [0,1].")
  if (!is.numeric(k_neighbors) || length(k_neighbors) != 1L || k_neighbors < 1)
    stop("`k_neighbors` must be a positive integer.")
  if (length(pilot_data) == 0L) stop("`pilot_data` is empty.")
  if (length(sim_data) == 0L)   stop("`sim_data` is empty.")

  if (.Platform$OS.type == "windows" && n_cores > 1L) {
    if (verbose) message(">>> Parallelisation limited to 1 core on Windows.")
    n_cores <- 1L
  }
  n_cores <- max(1L, as.integer(n_cores))

  # --------------------------------------------------------------------------
  # Helpers
  # --------------------------------------------------------------------------
  normalize_coords_01 <- function(coords_df) {
    xmin <- min(coords_df$x, na.rm = TRUE); xmax <- max(coords_df$x, na.rm = TRUE)
    ymin <- min(coords_df$y, na.rm = TRUE); ymax <- max(coords_df$y, na.rm = TRUE)
    max_range <- max(xmax - xmin, ymax - ymin)
    if (!is.finite(max_range) || max_range <= 0) return(coords_df)
    out <- coords_df
    out$x <- (out$x - xmin) / max_range
    out$y <- (out$y - ymin) / max_range
    out
  }

  # Row-wise standardization. Uses the proven transpose idiom, which coerces its
  # input to a matrix and is robust to degenerate single-row / single-column input
  # (a bare vector reaching rowMeans() is what triggered the earlier crash).
  safe_scale_rows <- function(mat) {
    m <- t(scale(t(as.matrix(mat)), center = TRUE, scale = TRUE))
    m[!is.finite(m)] <- 0
    m
  }

  .get_expr_sparse <- function(obj) {
    if (!is.null(obj$logcounts)) return(obj$logcounts)
    if (!is.null(obj$counts))    return(log1p(obj$counts))
    stop("Each sample must contain either `logcounts` or `counts`.")
  }
  .get_raw_counts_sparse <- function(obj) {
    if (!is.null(obj$counts)) return(obj$counts)
    stop("Each sample must contain `counts`.")
  }

  .pick_reference_pilot <- function(pilot_data) {
    score_mat <- vapply(seq_along(pilot_data), function(i) {
      coords <- pilot_data[[i]]$coords
      n_dom  <- if (!is.null(coords) && "domain" %in% names(coords))
        length(unique(as.character(coords$domain))) else 0
      n_spot <- if (!is.null(coords)) nrow(coords) else 0
      c(n_dom = n_dom, n_spot = n_spot)
    }, numeric(2))
    ord <- order(score_mat["n_dom", ], score_mat["n_spot", ], decreasing = TRUE)
    pilot_data[[ord[1]]]
  }

  # Top-k scores (U %*% diag(d)); RSpectra truncated SVD when available.
  .fast_scores <- function(Xc, k) {
    max_k <- min(dim(Xc))
    k <- max(1L, min(as.integer(k), max_k))
    if (k < (max_k - 1L) && requireNamespace("RSpectra", quietly = TRUE)) {
      sv <- tryCatch(RSpectra::svds(Xc, k = k), error = function(e) NULL)
      if (!is.null(sv)) return(sv$u %*% diag(sv$d, nrow = k, ncol = k))
    }
    pr <- stats::prcomp(Xc, center = FALSE, scale. = FALSE)
    pr$x[, seq_len(min(k, ncol(pr$x))), drop = FALSE]
  }

  # --------------------------------------------------------------------------
  # BANKSY embedding builder (vectorized via sparse neighbour-weight matrix)
  # --------------------------------------------------------------------------
  build_banksy <- function(expr_mat, pos_xy, lambda, k_neighbors,
                           expr_is_log = TRUE, use_G = TRUE,
                           balance_mode = c("energy_proportional", "equal_split", "paper_mu"),
                           mu = 1.5, eps = 1e-12) {
    balance_mode <- match.arg(balance_mode)
    block_energy <- function(A) if (is.null(A)) 0 else mean(A^2, na.rm = TRUE)

    logm <- if (isTRUE(expr_is_log)) as.matrix(expr_mat) else log1p(as.matrix(expr_mat))
    C <- safe_scale_rows(logm)
    n_genes <- nrow(C); n_spots <- ncol(C)
    if (n_spots <= 0L) stop("Expression matrix has zero spots.")
    if (nrow(pos_xy) != n_spots)
      stop(sprintf(paste0("Spot-count mismatch: expression has %d spots (columns) ",
                          "but coordinates have %d rows."), n_spots, nrow(pos_xy)))

    k_use <- min(as.integer(k_neighbors), max(1L, n_spots - 1L))

    if (n_spots == 1L) {
      M <- matrix(0, n_genes, n_spots)
      G <- if (use_G) matrix(0, n_genes, n_spots) else NULL
    } else {
      knn  <- FNN::get.knn(pos_xy, k = k_use)
      dist <- knn$nn.dist
      idx  <- knn$nn.index

      midk <- ceiling(k_use / 2)
      sig  <- dist[, midk]
      sig0 <- (!is.finite(sig) | sig <= 0)
      if (any(sig0)) {
        fallback <- mean(dist[is.finite(dist) & dist > 0], na.rm = TRUE)
        if (!is.finite(fallback) || fallback <= 0) fallback <- 1
        sig[sig0] <- fallback
      }

      w  <- exp(-sweep(dist^2, 1, sig^2, FUN = "/"))
      rs <- rowSums(w); rs[!is.finite(rs) | rs <= 0] <- 1
      w  <- w / rs                               # rows sum to 1 (used in G identity)

      # Sparse neighbour-weight matrix: Wsp[a, i] = sum_{j: idx[i,j]=a} w[i,j]
      Wsp <- Matrix::sparseMatrix(
        i    = as.integer(idx),                  # column-major -> neighbour (row) index
        j    = rep.int(seq_len(n_spots), k_use), # spot i (column)
        x    = as.numeric(w),
        dims = c(n_spots, n_spots)
      )

      M_raw <- as.matrix(C %*% Wsp)              # neighbour mean (one sparse multiply)
      M <- safe_scale_rows(M_raw)

      G <- NULL
      if (use_G) {
        Csq <- C * C
        G2  <- as.matrix(Csq %*% Wsp) - 2 * C * M_raw + Csq
        G   <- safe_scale_rows(sqrt(pmax(G2, 0)))   # G2 first: preserve matrix dims
      }
    }

    eC <- max(block_energy(C), eps)
    eM <- max(block_energy(M), eps)

    if (is.null(G) || !use_G) {
      aC <- sqrt((1 - lambda) / eC); aM <- sqrt(lambda / eM)
      B  <- rbind(aC * C, aM * M)
      return(list(B = B, balance = list(mode = "2-block", eC = eC, eM = eM, eG = 0,
                                        aC = aC, aM = aM, aG = 0, pi_M = 1, pi_G = 0)))
    }

    eG <- max(block_energy(G), eps)

    if (balance_mode == "paper_mu") {
      aC <- sqrt(1 - lambda); aM <- sqrt(lambda / mu); aG <- sqrt(lambda / (2 * mu))
      B  <- rbind(aC * C, aM * M, aG * G)
      return(list(B = B, balance = list(mode = "paper_mu", mu = mu, eC = eC, eM = eM,
                                        eG = eG, aC = aC, aM = aM, aG = aG,
                                        pi_M = NA_real_, pi_G = NA_real_)))
    }

    if (balance_mode == "equal_split") {
      pi_M <- 0.5; pi_G <- 0.5
    } else {
      denom <- max(eM + eG, eps); pi_M <- eM / denom; pi_G <- eG / denom
    }

    aC <- sqrt((1 - lambda) / eC)
    aM <- sqrt((lambda * pi_M) / eM)
    aG <- sqrt((lambda * pi_G) / eG)
    B  <- rbind(aC * C, aM * M, aG * G)
    list(B = B, balance = list(mode = balance_mode, eC = eC, eM = eM, eG = eG,
                               aC = aC, aM = aM, aG = aG, pi_M = pi_M, pi_G = pi_G))
  }

  # --------------------------------------------------------------------------
  # Pilot pre-processing (once, shared across samples)
  # --------------------------------------------------------------------------
  if (verbose) message(">>> Pre-processing pilot data ...")
  ref_pilot       <- .pick_reference_pilot(pilot_data)
  p_counts_sparse <- .get_raw_counts_sparse(ref_pilot)
  p_expr_sparse   <- .get_expr_sparse(ref_pilot)
  p_coordsraw     <- ref_pilot$coords

  if (!all(c("x", "y") %in% colnames(p_coordsraw)))
    stop("Reference pilot coords must have columns x and y.")
  if (!("domain" %in% colnames(p_coordsraw)))
    stop("Reference pilot coords must have column 'domain'.")
  if (is.null(rownames(p_counts_sparse)) || is.null(rownames(p_expr_sparse)))
    stop("Reference pilot expression/count matrices must have rownames.")

  p_coords <- p_coordsraw
  if (scale_pilot_coords) p_coords <- normalize_coords_01(p_coords)

  sim_counts1 <- sim_data[[1]]$counts
  if (is.null(dim(sim_counts1)))
    stop("sim_data[[1]]$counts is not a matrix.")
  sim_genes <- rownames(sim_counts1)
  if (is.null(sim_genes)) stop("sim_data[[1]]$counts must have rownames.")

  common_genes <- intersect(rownames(p_expr_sparse), sim_genes)
  if (length(common_genes) == 0L) stop("No matching genes between pilot and simulated data.")

  p_expr <- as.matrix(p_expr_sparse[common_genes, , drop = FALSE])
  p_pos  <- as.matrix(p_coords[, c("x", "y"), drop = FALSE])

  p_banksy <- build_banksy(p_expr, p_pos, lambda, k_neighbors,
                           expr_is_log = TRUE, use_G = use_G,
                           balance_mode = balance_mode, mu = mu)
  B_pilot <- p_banksy$B

  p_domains   <- as.character(p_coords$domain)
  unique_doms <- sort(unique(p_domains))
  n_k         <- length(unique_doms)

  # --------------------------------------------------------------------------
  # Per-sample worker
  # --------------------------------------------------------------------------
  .process_one <- function(i) {
    samp            <- sim_data[[i]]
    s_counts_sparse <- .get_raw_counts_sparse(samp)
    s_expr_sparse   <- .get_expr_sparse(samp)
    s_coords0       <- samp$coords

    latent_truth <- if ("tilde_d" %in% names(s_coords0)) as.character(s_coords0$tilde_d)
    else if ("domain" %in% names(s_coords0)) as.character(s_coords0$domain)
    else rep(NA_character_, nrow(s_coords0))

    use_genes <- intersect(common_genes, rownames(s_expr_sparse))
    if (length(use_genes) == 0L) return(NULL)

    s_counts <- as.matrix(s_counts_sparse[use_genes, , drop = FALSE])
    s_expr   <- as.matrix(s_expr_sparse[use_genes, , drop = FALSE])
    s_coords <- s_coords0
    if (scale_sim_coords) s_coords <- normalize_coords_01(s_coords)
    s_pos <- as.matrix(s_coords[, c("x", "y"), drop = FALSE])

    s_banksy <- build_banksy(s_expr, s_pos, lambda, k_neighbors,
                             expr_is_log = TRUE, use_G = use_G,
                             balance_mode = balance_mode, mu = mu)
    B_sim <- s_banksy$B

    if (!identical(use_genes, common_genes)) {
      p_expr_i  <- p_expr[use_genes, , drop = FALSE]
      B_pilot_i <- build_banksy(p_expr_i, p_pos, lambda, k_neighbors,
                                expr_is_log = TRUE, use_G = use_G,
                                balance_mode = balance_mode, mu = mu)$B
    } else {
      B_pilot_i <- B_pilot
    }

    B_joint <- cbind(B_pilot_i, B_sim)

    # Vectorized zero-variance filter (replaces apply(., 1, var))
    rmean <- rowMeans(B_joint)
    keep  <- (rowSums((B_joint - rmean)^2) > 0)
    if (!any(keep)) return(NULL)

    pca_in      <- t(B_joint[keep, , drop = FALSE])
    npc_use     <- min(n_pcs, ncol(pca_in), nrow(pca_in))
    embed_total <- .fast_scores(pca_in, npc_use)

    n_p         <- ncol(B_pilot_i)
    embed_pilot <- embed_total[seq_len(n_p), , drop = FALSE]
    embed_sim   <- embed_total[(n_p + 1):nrow(embed_total), , drop = FALSE]

    start_centers <- matrix(0, nrow = n_k, ncol = ncol(embed_pilot))
    rownames(start_centers) <- unique_doms
    for (d in unique_doms) {
      idx_d <- which(p_domains == d)
      start_centers[d, ] <- if (length(idx_d) > 1L)
        colMeans(embed_pilot[idx_d, , drop = FALSE]) else embed_pilot[idx_d, , drop = FALSE]
    }
    guided_names <- rownames(start_centers)

    if (hard_guidance) {
      dmat   <- as.matrix(stats::dist(rbind(embed_sim, start_centers)))
      n_sim  <- nrow(embed_sim)
      d_sc   <- dmat[seq_len(n_sim), (n_sim + 1):(n_sim + n_k), drop = FALSE]
      raw_clusters <- max.col(-d_sc)
    } else {
      km_res <- tryCatch(stats::kmeans(embed_sim, centers = start_centers, nstart = 1),
                         error = function(e) stats::kmeans(embed_sim, centers = n_k, nstart = 25))
      raw_clusters <- km_res$cluster
    }
    pred_labels_named <- guided_names[raw_clusters]

    if (do_hungarian && all(!is.na(latent_truth))) {
      all_doms <- unique(c(latent_truth, guided_names))
      profit <- matrix(0, nrow = length(all_doms), ncol = length(all_doms),
                       dimnames = list(all_doms, all_doms))
      for (td in all_doms) for (pd in all_doms)
        profit[td, pd] <- sum(latent_truth == td & pred_labels_named == pd)
      cost   <- max(profit) - profit
      assign <- clue::solve_LSAP(cost)
      pred_to_truth <- stats::setNames(rep(NA_character_, length(all_doms)), all_doms)
      pred_to_truth[all_doms[as.integer(assign)]] <- all_doms
      final_labels <- pred_to_truth[pred_labels_named]
      na_idx <- is.na(final_labels)
      if (any(na_idx)) final_labels[na_idx] <- pred_labels_named[na_idx]
    } else {
      final_labels <- pred_labels_named
    }

    coords_df <- s_coords0
    coords_df$hat_d <- final_labels
    if ("tilde_d" %in% names(coords_df)) {
      coords_df$is_correct <- (as.character(coords_df$hat_d) == as.character(coords_df$tilde_d))
      match_acc <- mean(coords_df$is_correct)
    } else if ("domain" %in% names(coords_df)) {
      coords_df$is_correct <- (as.character(coords_df$hat_d) == as.character(coords_df$domain))
      match_acc <- mean(coords_df$is_correct)
    } else {
      coords_df$is_correct <- NA; match_acc <- NA_real_
    }

    list(counts = s_counts, coords = coords_df,
         group = if (!is.null(samp$meta$group)) samp$meta$group else NA,
         sample_id = if (!is.null(samp$meta$sample_id)) samp$meta$sample_id else i,
         global_accuracy = match_acc, n_domains = n_k,
         hard_guidance = hard_guidance, scale_sim_coords = scale_sim_coords,
         banksy_use_G = use_G, banksy_balance_mode = balance_mode)
  }

  # --------------------------------------------------------------------------
  # Execute (parallel across samples; serial path shows a progress bar)
  # --------------------------------------------------------------------------
  if (n_cores == 1L) {
    pb <- if (verbose) utils::txtProgressBar(min = 0, max = length(sim_data), style = 3) else NULL
    processed_samples <- vector("list", length(sim_data))
    for (i in seq_along(sim_data)) {
      processed_samples[[i]] <- .process_one(i)
      if (!is.null(pb)) utils::setTxtProgressBar(pb, i)
    }
    if (!is.null(pb)) close(pb)
  } else {
    if (verbose) message(sprintf(">>> Recovering domains for %d samples (n_cores=%d)...",
                                 length(sim_data), n_cores))
    processed_samples <- parallel::mclapply(seq_along(sim_data), .process_one,
                                            mc.cores = n_cores, mc.set.seed = TRUE)
  }
  if (verbose) message("\nDone.")

  processed_samples[!vapply(processed_samples, is.null, logical(1))]
}
