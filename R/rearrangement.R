#' @title Rearrange Synthetic Spots onto Pilot Domain Geometry (Per-Domain OT)
#' @description
#' Reshapes blob-like simulated domains so that each domain occupies the same
#' spatial region as the corresponding domain in a reference pilot sample, via a
#' two-stage per-domain optimal transport (exact 2D-Gaussian affine alignment,
#' then bijective nearest-neighbour snapping onto real pilot spots).
#'
#' \strong{Composition-aware mode (\code{boundary_shift = TRUE}).} The default
#' per-domain map sends each synthetic domain onto the SAME domain's pilot
#' footprint, so an increased case-group target proportion is rendered as higher
#' \emph{density} inside the fixed (control) pilot target footprint rather than a
#' larger \emph{area}. When \code{boundary_shift = TRUE} (with
#' \code{target_domain} and \code{reference_domain} supplied), the target and
#' reference pilot spots are pooled and re-partitioned per section so that the
#' target \strong{territory area fraction matches that section's own target
#' proportion} \eqn{f = n_{target}/(n_{target}+n_{reference})}: the target region
#' is grown contiguously into the nearest reference spots until it holds an
#' \eqn{f} fraction of the pool. Synthetic target spots are then placed in this
#' (enlarged, for case) territory and reference spots in the remainder. Labels
#' and per-domain counts are untouched (so the LOR statistic is unchanged); only
#' the spatial footprint reflects the effect. Control sections (baseline \eqn{f})
#' look like the pilot; case sections show the target band expanding into the
#' reference. All other domains map by the standard per-domain OT.
#'
#' @param synthetic_data List of simulated samples; each has \code{coords} with
#'   \code{x}, \code{y}, the \code{match_by} label, and \code{meta$group}/\code{meta$sample_id}.
#' @param pilot_data List of pilot samples; each has \code{coords} with \code{x},
#'   \code{y}, \code{domain}, \code{group}, \code{sample_id}.
#' @param match_by \code{"hat_d"} (default) or \code{"domain"}.
#' @param jitter_frac Numeric. Oversampling jitter as a fraction of local spacing.
#'   Default 0.4.
#' @param k_nn Integer. Candidate neighbours for the greedy assignment. Default 80.
#' @param seed Integer base seed (global RNG restored on exit). Default 1.
#' @param verbose Logical. Default \code{TRUE}.
#' @param boundary_shift Logical. If TRUE, render the target/reference composition
#'   difference as a moving boundary (see Description). Default FALSE (original
#'   per-domain behaviour, unchanged).
#' @param target_domain,reference_domain Character. The domain pair whose boundary
#'   moves when \code{boundary_shift = TRUE} (e.g. "WM" and "Layer6"). Must match
#'   the labels in the \code{match_by} column.
#'
#' @return The input list with \code{coords$x}/\code{coords$y} replaced by the
#'   rearranged coordinates; originals kept as \code{x_raw}/\code{y_raw}, and
#'   \code{matched_pilot_id} records the reference section.
#'
#' @seealso \code{\link{pBANSKY}}, \code{\link{simulatespaCraft}}
#'
#' @importFrom stats cov rnorm median
#' @export
#'
rearrangeSyntheticToPilot <- function(synthetic_data,
                                      pilot_data,
                                      match_by         = c("hat_d", "domain"),
                                      jitter_frac      = 0.4,
                                      k_nn             = 80L,
                                      seed             = 1L,
                                      verbose          = TRUE,
                                      boundary_shift   = FALSE,
                                      target_domain    = NULL,
                                      reference_domain = NULL) {

  match_by <- match.arg(match_by)
  if (!requireNamespace("RANN", quietly = TRUE) &&
      !requireNamespace("FNN", quietly = TRUE))
    stop("Package 'RANN' or 'FNN' is required.")
  if (isTRUE(boundary_shift) && (is.null(target_domain) || is.null(reference_domain)))
    stop("boundary_shift = TRUE requires target_domain and reference_domain.")

  if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
    .old <- get(".Random.seed", envir = globalenv(), inherits = FALSE)
    on.exit(assign(".Random.seed", .old, envir = globalenv()), add = TRUE)
  }
  set.seed(seed)

  out <- vector("list", length(synthetic_data))
  names(out) <- names(synthetic_data)

  for (si in seq_along(synthetic_data)) {
    sid  <- names(synthetic_data)[si] %||% paste0("sim_", si)
    samp <- synthetic_data[[si]]
    S_xy  <- as.matrix(samp$coords[, c("x", "y")])
    S_dom <- as.character(samp$coords[[match_by]])
    grp   <- samp$meta$group

    pool <- which(vapply(pilot_data, function(p) isTRUE(p$group == grp), logical(1)))
    if (length(pool) == 0) pool <- seq_along(pilot_data)
    ref  <- pilot_data[[pool[sample.int(length(pool), 1)]]]
    P_xy  <- as.matrix(ref$coords[, c("x", "y")])
    P_dom <- as.character(ref$coords$domain)

    if (verbose) message(sprintf(
      "  -> Morphing %s (N=%d) onto pilot %s%s ...",
      sid, nrow(S_xy), ref$sample_id %||% paste0("pilot_", pool[1]),
      if (boundary_shift) " [boundary-shift]" else ""))

    new_xy <- if (isTRUE(boundary_shift))
      .morph_sample_DA(S_xy, S_dom, P_xy, P_dom, target_domain, reference_domain,
                       jitter_frac, k_nn)
    else
      .morph_sample(S_xy, S_dom, P_xy, P_dom, jitter_frac, k_nn)

    samp$coords$x_raw <- samp$coords$x
    samp$coords$y_raw <- samp$coords$y
    samp$coords$x <- new_xy[, 1]
    samp$coords$y <- new_xy[, 2]
    samp$matched_pilot_id <- ref$sample_id
    out[[si]] <- samp
  }
  out
}

# ---- internal helpers --------------------------------------------------------

.knnx <- function(data, query, k) {
  k <- min(k, nrow(data))
  if (requireNamespace("RANN", quietly = TRUE)) RANN::nn2(data, query, k = k)$nn.idx
  else FNN::get.knnx(data = data, query = query, k = k)$nn.index
}

.knnx_dist <- function(data, query, k) {
  k <- min(k, nrow(data))
  if (requireNamespace("RANN", quietly = TRUE)) RANN::nn2(data, query, k = k)$nn.dists
  else FNN::get.knnx(data = data, query = query, k = k)$nn.dist
}

.sym_pow2 <- function(M, p) {
  M <- (M + t(M)) / 2
  e <- eigen(M, symmetric = TRUE)
  e$vectors %*% diag(pmax(e$values, 1e-10)^p, 2) %*% t(e$vectors)
}

# Exact OT map between two 2D Gaussians (match mean + covariance).
.gaussian_affine <- function(S, P) {
  if (nrow(S) < 3 || nrow(P) < 3)
    return(sweep(S, 2, colMeans(S)) + matrix(colMeans(P), nrow(S), 2, byrow = TRUE))
  muS <- colMeans(S); muP <- colMeans(P)
  Ss <- stats::cov(S) + diag(1e-9, 2)
  Sp <- stats::cov(P) + diag(1e-9, 2)
  Ss_h  <- .sym_pow2(Ss, 0.5)
  Ss_ih <- .sym_pow2(Ss, -0.5)
  A     <- Ss_ih %*% .sym_pow2(Ss_h %*% Sp %*% Ss_h, 0.5) %*% Ss_ih
  sweep(S, 2, muS) %*% A + matrix(muP, nrow(S), 2, byrow = TRUE)
}

.resample_region <- function(P, n, jitter_frac) {
  m <- nrow(P)
  if (m == 0L) return(matrix(0, n, 2))
  if (n <= m) return(P[sample.int(m, n), , drop = FALSE])
  idx <- c(seq_len(m), sample.int(m, n - m, replace = TRUE))
  out <- P[idx, , drop = FALSE]
  if (m >= 2) {
    nnd <- if (requireNamespace("RANN", quietly = TRUE))
      RANN::nn2(P, P, k = 2)$nn.dists[, 2] else FNN::get.knn(P, k = 1)$nn.dist[, 1]
    s <- jitter_frac * stats::median(nnd[is.finite(nnd) & nnd > 0])
    if (is.finite(s) && s > 0) out <- out + matrix(stats::rnorm(2 * n, 0, s), n, 2)
  }
  out
}

.greedy_assign <- function(query, ref, nn_idx) {
  if (exists("greedy_assign_cpp", mode = "function"))
    return(greedy_assign_cpp(query, ref, nn_idx))
  n <- nrow(query); k <- ncol(nn_idx)
  used <- logical(nrow(ref)); assigned <- integer(n)
  for (i in seq_len(n)) {
    got <- FALSE
    for (j in seq_len(k)) {
      cand <- nn_idx[i, j]
      if (!used[cand]) { assigned[i] <- cand; used[cand] <- TRUE; got <- TRUE; break }
    }
    if (!got) {
      free <- which(!used)
      d2 <- (ref[free, 1] - query[i, 1])^2 + (ref[free, 2] - query[i, 2])^2
      b <- free[which.min(d2)]; assigned[i] <- b; used[b] <- TRUE
    }
  }
  assigned
}

# place a set of synthetic points (S) into a target territory (Pterr)
.place_into <- function(S, Pterr, n, jitter_frac, k_nn) {
  if (n == 0L) return(matrix(nrow = 0, ncol = 2))
  Pr  <- .resample_region(Pterr, n, jitter_frac)
  Sa  <- .gaussian_affine(S, Pr)
  asg <- .greedy_assign(Sa, Pr, .knnx(Pr, Sa, k_nn))
  Pr[asg, , drop = FALSE]
}

# Standard per-domain OT (original behaviour).
.morph_sample <- function(S_xy, S_dom, P_xy, P_dom, jitter_frac, k_nn) {
  N <- nrow(S_xy); new_xy <- matrix(NA_real_, N, 2)
  for (d in intersect(unique(S_dom), unique(P_dom))) {
    iS <- which(S_dom == d)
    P  <- P_xy[P_dom == d, , drop = FALSE]
    new_xy[iS, ] <- .place_into(S_xy[iS, , drop = FALSE], P, length(iS), jitter_frac, k_nn)
  }
  miss <- which(is.na(new_xy[, 1]))
  if (length(miss) > 0) {
    nn <- .knnx(P_xy, S_xy[miss, , drop = FALSE], 1)
    new_xy[miss, ] <- P_xy[nn[, 1], ]
  }
  new_xy
}

# Composition-aware OT: move the target/reference boundary to match this
# section's target proportion (grow target into nearest reference spots).
.morph_sample_DA <- function(S_xy, S_dom, P_xy, P_dom, target, reference,
                             jitter_frac, k_nn) {
  N <- nrow(S_xy); new_xy <- matrix(NA_real_, N, 2)

  haveTR <- (target %in% P_dom) && (reference %in% P_dom) &&
    (target %in% S_dom || reference %in% S_dom)

  if (haveTR) {
    iT_S <- which(S_dom == target); iR_S <- which(S_dom == reference)
    nT <- length(iT_S); nR <- length(iR_S)

    idxT_P <- which(P_dom == target); idxR_P <- which(P_dom == reference)
    P_pool <- P_xy[c(idxT_P, idxR_P), , drop = FALSE]
    core   <- P_xy[idxT_P, , drop = FALSE]                       # target core
    np     <- nrow(P_pool)

    # distance of each pool spot to the nearest target-core spot (0 on the core)
    d2t <- .knnx_dist(core, P_pool, 1)[, 1]

    # target territory = an f-fraction of the pool, grown outward from the core
    f      <- if ((nT + nR) > 0) nT / (nT + nR) else length(idxT_P) / np
    n_terr <- max(1L, min(np - 1L, as.integer(round(f * np))))
    ord    <- order(d2t)
    terrT  <- P_pool[ord[seq_len(n_terr)], , drop = FALSE]
    terrR  <- P_pool[ord[(n_terr + 1L):np], , drop = FALSE]

    new_xy[iT_S, ] <- .place_into(S_xy[iT_S, , drop = FALSE], terrT, nT, jitter_frac, k_nn)
    new_xy[iR_S, ] <- .place_into(S_xy[iR_S, , drop = FALSE], terrR, nR, jitter_frac, k_nn)
  }

  done <- if (haveTR) c(target, reference) else character(0)
  for (d in setdiff(intersect(unique(S_dom), unique(P_dom)), done)) {
    iS <- which(S_dom == d)
    P  <- P_xy[P_dom == d, , drop = FALSE]
    new_xy[iS, ] <- .place_into(S_xy[iS, , drop = FALSE], P, length(iS), jitter_frac, k_nn)
  }

  miss <- which(is.na(new_xy[, 1]))
  if (length(miss) > 0) {
    nn <- .knnx(P_xy, S_xy[miss, , drop = FALSE], 1)
    new_xy[miss, ] <- P_xy[nn[, 1], ]
  }
  new_xy
}


# =============================================================================
# Plot method
# =============================================================================

# =============================================================================
# Plot method
# =============================================================================

# Quiet R CMD check for ggplot2/dplyr non-standard-evaluation column names.
utils::globalVariables(c(".col", "group_label", "sample_id"))

#' @include createspaCraftObject.R
#' @title Plot Spatial Domains for a spaCraft Object
#' @description
#' A unified plotting method for \code{spaCraft} objects. Visualizes the
#' spatial distribution of domains for both pilot data and synthetic data. It
#' automatically extracts unified domain levels so that colours are consistent
#' across all plots.
#'
#' @param x A \code{spaCraft} object.
#' @param y Not used. Required for S4 generic compatibility.
#' @param type Character. Which data to plot: \code{"pilot"} (ground truth),
#'   \code{"synthetic_raw"} (before rearrangement, continuous), or
#'   \code{"synthetic_mapped"} (after rearrangement, KD-tree snapped).
#' @param color_by Character. For synthetic data, colour by the recovered domains
#'   (\code{"hat_d"}, default) or by truth (\code{"domain"}).
#' @param n_per_group Integer. Number of samples to show per group. Default 3.
#' @param ncol Integer. Number of columns in the facet wrap. Default 3.
#' @param sample_ids Optional character vector of specific sample IDs to plot.
#' @param cex Numeric. Point-size expansion factor. Default 1.0.
#' @param alpha Numeric. Point transparency in \eqn{[0,1]}. Default 0.9.
#' @param title Character. Custom plot title; if \code{NULL}, a default is used.
#' @param subtitle Character. Custom subtitle; if \code{NULL}, a default is used.
#' @param colors Character vector of colours. If \code{NULL}, uses the "Set1" palette.
#' @param reverse_y Logical. If \code{TRUE}, reverses the y-axis to match image coordinates. Default \code{FALSE}.
#' @param ... Additional arguments (currently ignored).
#'
#' @return A \code{ggplot} object.
#'
#' @import ggplot2
#' @importFrom dplyr %>% distinct group_by slice_head pull
#' @importFrom graphics plot
#' @export
#'
setMethod("plot", signature(x = "spaCraft", y = "missing"),
          function(x, y,
                   type = c("pilot", "synthetic_raw", "synthetic_mapped"),
                   color_by = c("hat_d", "domain"),
                   n_per_group = 3,
                   ncol = 3,
                   sample_ids = NULL,
                   cex = 1.0,
                   alpha = 0.9,
                   title = NULL,
                   subtitle = NULL,
                   colors = NULL,
                   reverse_y = FALSE,
                   ...) {

            type <- match.arg(type)
            color_by <- match.arg(color_by)

            # Base point size (ggplot size scale), expanded by cex.
            base_size   <- 1.2
            actual_size <- base_size * cex

            # ------------------------------------------------------------------
            # 1. Extract unified domain levels (for perfect colour matching)
            # ------------------------------------------------------------------
            pilot_doms <- unlist(lapply(x@pilot_data, function(s) as.character(s$coords$domain)))
            syn_doms <- if (length(x@synthetic_data) > 0) {
              unlist(lapply(x@synthetic_data, function(s) as.character(s$coords[[color_by]])))
            } else { character(0) }

            all_domains <- sort(unique(c(pilot_doms, syn_doms)))

            if (is.null(colors) && length(all_domains) > 9) {
              warning("More than 9 unique domains found. The default 'Set1' palette supports up to 9 colors. Consider providing a custom 'colors' vector.")
            }

            # ------------------------------------------------------------------
            # 2. Set data source, titles, and axes based on 'type'
            # ------------------------------------------------------------------
            if (type == "pilot") {
              data_list <- x@pilot_data
              if (length(data_list) == 0) stop("No pilot data found in the object.")

              default_title    <- "Pilot Data: Ground Truth Domains"
              default_subtitle <- NULL
              legend_title     <- "Domain"
              col_var <- "domain"
              use_raw <- FALSE

            } else {
              data_list <- x@synthetic_data
              if (length(data_list) == 0) stop("No synthetic data found. Generate or map synthetic data first.")

              use_raw      <- (type == "synthetic_raw")
              title_suffix <- if (use_raw) "(Before Rearrangement)" else "(After Rearrangement)"

              default_title    <- paste("Synthetic Data:", title_suffix)
              default_subtitle <- "Top Row: Control | Bottom Row: Case"

              # Mathematical expression for "hat_d" in the legend.
              legend_title <- if (color_by == "hat_d") {
                expression(paste("Recovered (", hat(d), ")"))
              } else {
                "Truth (domain)"
              }
              col_var <- color_by
            }

            final_title    <- if (!is.null(title)) title else default_title
            final_subtitle <- if (!is.null(subtitle)) subtitle else default_subtitle

            # ------------------------------------------------------------------
            # 3. Assemble data frame (extraction and formatting)
            # ------------------------------------------------------------------
            df_list <- lapply(seq_along(data_list), function(i) {
              s  <- data_list[[i]]
              df <- s$coords

              # Extract meta safely.
              sid <- if (!is.null(s$sample_id)) s$sample_id else if (!is.null(s$meta$sample_id)) s$meta$sample_id else paste0("Sample_", i)
              grp <- if (!is.null(s$group)) s$group else if (!is.null(s$meta$group)) s$meta$group else NA

              # Decide X and Y (raw coordinates for "synthetic_raw").
              plot_x <- if (use_raw && "x_raw" %in% names(df)) df$x_raw else df$x
              plot_y <- if (use_raw && "y_raw" %in% names(df)) df$y_raw else df$y

              data.frame(
                x           = plot_x,
                y           = plot_y,
                .col        = factor(as.character(df[[col_var]]), levels = all_domains),
                sample_id   = sid,
                group_label = factor(ifelse(grp == 0, "Control", "Case"), levels = c("Control", "Case"))
              )
            })
            plot_df <- do.call(rbind, df_list)

            # ------------------------------------------------------------------
            # 4. Filter output samples (by n_per_group or sample_ids)
            # ------------------------------------------------------------------
            if (!is.null(sample_ids)) {
              plot_df <- plot_df[plot_df$sample_id %in% sample_ids, , drop = FALSE]
            } else if (!is.null(n_per_group)) {
              keep_ids <- plot_df %>%
                dplyr::distinct(group_label, sample_id) %>%
                dplyr::group_by(group_label) %>%
                dplyr::slice_head(n = n_per_group) %>%
                dplyr::pull(sample_id)
              plot_df <- plot_df[plot_df$sample_id %in% keep_ids, , drop = FALSE]
            }

            # ------------------------------------------------------------------
            # 5. Determine Colour Scale
            # ------------------------------------------------------------------
            color_scale <- if (is.null(colors)) {
              scale_color_brewer(palette = "Set1", drop = FALSE)
            } else {
              scale_color_manual(values = colors, drop = FALSE)
            }

            # ------------------------------------------------------------------
            # 6. Render with ggplot2
            # ------------------------------------------------------------------
            p <- ggplot(plot_df, aes(x = x, y = y, color = .col)) +
              geom_point(size = actual_size, alpha = alpha, stroke = 0) +
              facet_wrap(~ group_label + sample_id, scales = "free", ncol = ncol) +
              theme_bw(base_size = 12) +
              theme(
                panel.grid       = element_blank(),
                strip.background = element_rect(fill = "white", color = "black"),
                strip.text       = element_text(face = "bold", size = 10),
                legend.position  = "right"
              ) +
              color_scale +
              labs(
                title    = final_title,
                subtitle = final_subtitle,
                x        = "Spatial X",
                y        = "Spatial Y",
                color    = legend_title
              ) +
              # Make legend points slightly larger than plot points for readability.
              guides(color = guide_legend(override.aes = list(size = max(4, actual_size))))

            # ------------------------------------------------------------------
            # 7. Apply Y-axis Reversal if requested
            # ------------------------------------------------------------------
            if (isTRUE(reverse_y)) {
              p <- p + scale_y_reverse()
            }

            return(p)
          }
)
