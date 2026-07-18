# ==============================================================================
# Internal Helper: Calculate Fold Change Statistics
# ==============================================================================
# ==============================================================================
# Internal Helper: Calculate Fold Change Statistics
# ==============================================================================
#' @keywords internal
#' @noRd
#' @importFrom parallel mclapply
#' @importFrom stats median
geneSummary <- function(count_matrix, loc, n_cores, verbose = FALSE) {
  log_count <- log1p(suppressWarnings(as.matrix(count_matrix)))

  # [수정된 부분] 중복된 유전자 이름(row.names) 방어 로직 추가
  rn <- rownames(log_count)
  if (is.null(rn)) {
    rn <- paste0("Gene_", seq_len(nrow(log_count)))
  } else if (any(duplicated(rn))) {
    if (verbose) message("Warning: Duplicated gene names detected. Making them unique...")
    rn <- make.unique(rn)
  }
  rownames(log_count) <- rn
  # -----------------------------------------------------------

  domains <- sort(unique(as.character(loc$domain)))

  if (verbose) message(sprintf("Calculating fold changes across %d domains...", length(domains)))

  # OS-safe core configuration for mclapply
  if (.Platform$OS.type == "windows" && n_cores > 1L) {
    n_cores <- 1L
  }

  fc_results <- parallel::mclapply(domains, function(domain) {
    spot_idx <- which(as.character(loc$domain) == domain)

    # Edge case: Domain doesn't exist or takes up the entire sample
    if (length(spot_idx) == 0 || length(spot_idx) == ncol(log_count)) return(NULL)

    mat_in <- log_count[, spot_idx, drop = FALSE]
    mat_out <- log_count[, -spot_idx, drop = FALSE]

    m_in <- rowMeans(mat_in, na.rm = TRUE)
    m_out <- rowMeans(mat_out, na.rm = TRUE)

    # Vectorized out_low calculation
    if (requireNamespace("matrixStats", quietly = TRUE)) {
      med_out <- matrixStats::rowMedians(mat_out, na.rm = TRUE)
    } else {
      med_out <- apply(mat_out, 1, stats::median, na.rm = TRUE)
    }

    mask_low <- mat_out <= med_out

    mat_out_low <- mat_out
    mat_out_low[!mask_low] <- 0

    n_low <- rowSums(mask_low, na.rm = TRUE)
    n_low[n_low == 0] <- 1 # Prevent division by zero

    m_out_low <- rowSums(mat_out_low, na.rm = TRUE) / n_low

    # 이제 rownames가 고유(unique)하므로 에러가 발생하지 않습니다.
    data.frame(
      logFC = m_in - m_out,
      logFC_low = m_in - m_out_low,
      mean_in = m_in,
      mean_out = m_out,
      row.names = rownames(log_count),
      stringsAsFactors = FALSE
    )

  }, mc.cores = n_cores)

  names(fc_results) <- domains
  return(fc_results)
}

# ==============================================================================
# Exported Feature Selection Functions
# ==============================================================================

#' @title Select Domain-Informative Genes via Multi-Context Meta-Analysis
#'
#' @description
#' Identifies robust marker genes by evaluating spatial expression across multiple
#' biological contexts (e.g., NAT vs. CRC). This meta-analysis mitigates representation
#' bias, ensuring the final feature space reflects the full transcriptomic heterogeneity.
#'
#' @param object A \code{spaCraft} object.
#' @param rep_sample_ids Character vector of representative samples. At least one
#'   sample per condition should be provided. Defaults to the first sample if \code{NULL}.
#' @param logfc_cutoff Numeric. Minimum absolute log fold change (default: 0.5).
#' @param mean_in_cutoff Numeric. Minimum mean log-expression (default: 1.0).
#' @param max_num_gene Integer. Max top genes to extract per domain (default: 20).
#' @param n_cores Integer. CPU cores for parallel processing (default: 1).
#' @param combine_method Character. \code{"union"} (default) or \code{"intersection"}.
#' @param verbose Logical. Print progress.
#'
#' @return Updated \code{spaCraft} object with integrated features.
#' @importFrom utils head
#' @export
#'
#' @examples
#' \dontrun{
#' data("mini_spaCraft_obj", package = "spaCraft")
#' spaCraft_obj <- featureSelection(
#'   object = mini_spaCraft_obj,
#'   logfc_cutoff = 0.5,
#'   mean_in_cutoff = 0.8,
#'   max_num_gene = 3,
#'   n_cores = 8
#' )
#' }
featureSelection <- function(object,
                             rep_sample_ids = NULL,
                             logfc_cutoff = 0.5,
                             mean_in_cutoff = 1.0,
                             max_num_gene = 20,
                             n_cores = 1,
                             combine_method = c("union", "intersection"),
                             verbose = TRUE) {

  combine_method <- match.arg(combine_method)
  pilot_data <- object@pilot_data

  if (is.null(rep_sample_ids)) {
    rep_sample_ids <- names(pilot_data)[1]
    if (is.null(rep_sample_ids)) {
      names(pilot_data) <- paste0("Sample_", seq_along(pilot_data))
      rep_sample_ids <- names(pilot_data)[1]
    }
  }

  if (verbose) {
    message(sprintf(">>> Spatial Feature Selection: Processing %d Reference Context(s)", length(rep_sample_ids)))
    message(sprintf("    Using representative sample(s): %s", paste(rep_sample_ids, collapse = ", ")))
  }

  if (.Platform$OS.type == "windows" && n_cores > 1L) n_cores <- 1L
  meta_genes_list <- list()

  for (samp_id in rep_sample_ids) {
    if (verbose) message(sprintf("  -> Profiling sample: %s", samp_id))
    if (is.null(pilot_data[[samp_id]])) stop(sprintf("Sample ID '%s' not found in pilot_data.", samp_id))

    rep_sample <- pilot_data[[samp_id]]
    FC_list <- geneSummary(rep_sample$counts, rep_sample$coords, n_cores, verbose = FALSE)

    selected_genes_list <- lapply(FC_list, function(DF) {
      if (is.null(DF)) return(character(0))

      idx <- which(DF$mean_in >= mean_in_cutoff & abs(DF$logFC_low) >= logfc_cutoff)

      if (length(idx) == 0) {
        idx <- union(which(DF$mean_in >= mean_in_cutoff), which(abs(DF$logFC_low) >= logfc_cutoff))
      }

      subset_df <- DF[idx, , drop = FALSE]

      if (nrow(subset_df) > 0) {
        subset_df <- subset_df[order(abs(subset_df$logFC_low), decreasing = TRUE), , drop = FALSE]
        subset_df <- utils::head(subset_df, n = max_num_gene)
        return(rownames(subset_df))
      }
      return(character(0))
    })

    if (combine_method == "union") {
      meta_genes_list[[samp_id]] <- unique(unlist(selected_genes_list))
    } else {
      meta_genes_list[[samp_id]] <- if (length(selected_genes_list) > 0) Reduce(intersect, selected_genes_list) else character(0)
    }
  }

  final_genes <- unique(unlist(meta_genes_list))

  if (length(final_genes) == 0) {
    warning("No genes selected. Check your logFC and mean_in cutoffs.", call. = FALSE)
  } else if (verbose) {
    message(sprintf(">>> Integration Complete: Selected %d unique genes across all contexts.", length(final_genes)))
  }

  object@params_expression$top_genes <- final_genes
  return(object)
}

#' @title Select Spatially and Group-wise Stable Genes
#'
#' @description
#' Highly optimized version utilizing pre-filtering and vectorized operations
#' to bypass the computational bottleneck of gene summaries and nested loops.
#'
#' @param object A \code{spaCraft} object.
#' @param control_group_value Character. Value indicating the control group (default "0").
#' @param check_group_diff Logical. Check expression difference between Control and Case (default TRUE).
#' @param case_group_value Character. Value indicating the case group (default "1").
#' @param max_group_logfc_cut Numeric. Maximum allowed absolute logFC between Control and Case (default 0.1).
#' @param mean_global_cut Numeric. Minimum global mean log-expression (default 0.3).
#' @param pct_in_cut Numeric. Minimum detection rate within a domain (default 0.05).
#' @param min_domain_prop Numeric. Minimum proportion of domains that must meet \code{pct_in_cut} (default 0.8).
#' @param max_abs_logfc_cut Numeric. Maximum allowed spatial logFC variance (default 0.2).
#' @param cv_cut Numeric. Maximum allowed Robust CV across domain means (default 0.3).
#' @param max_num_gene Integer. Maximum number of genes to return (default 50).
#' @param use_low Logical. Use robust lower-bound logFC instead of standard logFC (default TRUE).
#' @param n_cores Integer. Number of cores for parallel processing.
#' @param verbose Logical. Print progress.
#'
#' @return Updated \code{spaCraft} object with \code{stable_genes} and \code{stable_gene_stats} populated.
#' @importFrom stats sd setNames
#' @importFrom utils head
#' @export
#'
#' @examples
#' \dontrun{
#' data("mini_spaCraft_obj", package = "spaCraft")
#' spaCraft_obj <- featureSelectionStable(
#'   object = mini_spaCraft_obj,
#'   max_group_logfc_cut = 0.1,
#'   mean_global_cut = 0.1,
#'   min_domain_prop = 0.2,
#'   max_abs_logfc_cut = 0.1,
#'   cv_cut = 0.2,
#'   max_num_gene = 10,
#'   use_low = FALSE,
#'   n_cores = 8
#' )
#' }
featureSelectionStable <- function(object,
                                   control_group_value = "0",
                                   check_group_diff = TRUE,
                                   case_group_value = "1",
                                   max_group_logfc_cut = 0.1,
                                   mean_global_cut = 0.3,
                                   pct_in_cut = 0.05,
                                   min_domain_prop = 0.8,
                                   max_abs_logfc_cut = 0.2,
                                   cv_cut = 0.3,
                                   max_num_gene = 50,
                                   use_low = TRUE,
                                   n_cores = 1,
                                   verbose = TRUE) {

  pilot_data <- object@pilot_data
  num_total_samples <- length(pilot_data)

  if (num_total_samples == 0) stop("pilot_data is empty.")

  # 1. Identify Reference sample(s)
  ctrl_idx <- which(vapply(pilot_data, function(x) as.character(x$group) == as.character(control_group_value), logical(1)))

  if (length(ctrl_idx) == 0) {
    if (num_total_samples == 1) {
      warning(sprintf("Control group '%s' not found. Using the only available sample as reference.", control_group_value), call. = FALSE)
      ctrl_idx <- 1
    } else {
      stop(sprintf("Control group (group='%s') sample not found.", control_group_value))
    }
  }

  num_ctrl <- length(ctrl_idx)
  if (verbose) message(sprintf(">>> Spatial Stability: Processing %d Reference sample(s)", num_ctrl))

  sample_results <- list()

  for (i in ctrl_idx) {
    samp <- pilot_data[[i]]
    if (verbose) message(sprintf("  -> Profiling sample: %s", samp$sample_id))

    count_matrix <- suppressWarnings(as.matrix(samp$counts))
    loc <- samp$coords
    genes <- rownames(count_matrix)
    domains <- sort(unique(as.character(loc$domain)))

    log_count <- log1p(count_matrix)
    mean_global <- rowMeans(log_count, na.rm = TRUE)
    names(mean_global) <- genes

    mu_d <- matrix(NA_real_, nrow = length(genes), ncol = length(domains), dimnames = list(genes, domains))
    pct_in_mat <- matrix(NA_real_, nrow = length(genes), ncol = length(domains), dimnames = list(genes, domains))

    for (d in domains) {
      idx <- which(as.character(loc$domain) == d)
      if (length(idx) < 2) next
      mu_d[, d] <- rowMeans(log_count[, idx, drop = FALSE], na.rm = TRUE)
      pct_in_mat[, d] <- rowMeans(count_matrix[, idx, drop = FALSE] > 0, na.rm = TRUE)
    }

    mu_bar <- rowMeans(mu_d, na.rm = TRUE)
    mu_sd  <- apply(mu_d, 1, stats::sd, na.rm = TRUE)
    mu_cv_robust <- mu_sd / (mu_bar + 0.1)

    pass_global <- (mean_global >= mean_global_cut)
    pass_domain_expr <- apply(pct_in_mat, 1, function(x) {
      sum(is.finite(x) & x >= pct_in_cut) / length(domains) >= min_domain_prop
    })
    pass_cv <- is.finite(mu_cv_robust) & (mu_cv_robust <= cv_cut)

    keep_fast_idx <- which(pass_global & pass_domain_expr & pass_cv)
    keep_fast_genes <- genes[keep_fast_idx]

    if (length(keep_fast_genes) == 0) {
      sample_results[[samp$sample_id]] <- data.frame()
      next
    }

    if (verbose) message(sprintf("     Running logFC calculations on %d pre-filtered genes...", length(keep_fast_genes)))
    count_matrix_sub <- count_matrix[keep_fast_genes, , drop = FALSE]
    FC_list <- geneSummary(count_matrix_sub, loc, n_cores = n_cores, verbose = FALSE)

    pick_fc_col <- if (use_low) "logFC_low" else "logFC"

    fc_mat <- sapply(domains, function(d) {
      df <- FC_list[[d]]
      vals <- stats::setNames(rep(NA_real_, length(keep_fast_genes)), keep_fast_genes)
      if (!is.null(df)) {
        common <- intersect(keep_fast_genes, rownames(df))
        vals[common] <- as.numeric(df[common, pick_fc_col])
      }
      return(vals)
    })

    if (!is.matrix(fc_mat)) fc_mat <- matrix(fc_mat, nrow = 1, dimnames = list(keep_fast_genes, domains))

    max_abs_logfc <- apply(abs(fc_mat), 1, max, na.rm = TRUE)
    pass_stable <- is.finite(max_abs_logfc) & (max_abs_logfc <= max_abs_logfc_cut)
    keep_final_genes <- keep_fast_genes[pass_stable]

    if (length(keep_final_genes) > 0) {
      cand <- data.frame(
        gene = keep_final_genes,
        mean_global = mean_global[keep_final_genes],
        max_abs_logfc = max_abs_logfc[keep_final_genes],
        cv_robust = mu_cv_robust[keep_final_genes],
        stringsAsFactors = FALSE
      )
      rownames(cand) <- cand$gene
      sample_results[[samp$sample_id]] <- cand
    } else {
      sample_results[[samp$sample_id]] <- data.frame()
    }
  }

  # 2. Integrate multi-sample results safely
  if (num_ctrl == 1) {
    avg_stats <- sample_results[[1]]
    common_genes <- rownames(avg_stats)
  } else {
    common_genes <- Reduce(intersect, lapply(sample_results, rownames))
    if (length(common_genes) > 0) {
      avg_stats <- data.frame(gene = common_genes, stringsAsFactors = FALSE, row.names = common_genes)
      avg_stats$mean_global <- rowMeans(do.call(cbind, lapply(sample_results, function(res) res[common_genes, "mean_global"])))
      avg_stats$max_abs_logfc <- rowMeans(do.call(cbind, lapply(sample_results, function(res) res[common_genes, "max_abs_logfc"])))
      avg_stats$cv_robust <- rowMeans(do.call(cbind, lapply(sample_results, function(res) res[common_genes, "cv_robust"])))
    }
  }

  if (length(common_genes) == 0) {
    warning("No spatially stable genes satisfied the given criteria.", call. = FALSE)
    object@params_expression$stable_genes <- character(0)
    return(object)
  }

  # 3. Group LogFC validation against Case group
  if (check_group_diff) {
    case_idx <- which(vapply(pilot_data, function(x) as.character(x$group) == as.character(case_group_value), logical(1)))

    if (length(case_idx) > 0 && !identical(ctrl_idx, case_idx)) {
      if (verbose) message(sprintf(">>> Group Stability: Checking LogFC against %d Case sample(s)", length(case_idx)))

      case_means_list <- list()
      for (i in case_idx) {
        samp <- pilot_data[[i]]
        valid_genes <- intersect(common_genes, rownames(samp$counts))
        X <- log1p(suppressWarnings(as.matrix(samp$counts[valid_genes, , drop = FALSE])))
        case_means_list[[samp$sample_id]] <- stats::setNames(rowMeans(X, na.rm = TRUE), valid_genes)
      }

      mean_case_global <- if (length(case_idx) == 1) {
        case_means_list[[1]]
      } else {
        rowMeans(do.call(cbind, case_means_list), na.rm = TRUE)
      }

      valid_genes_final <- intersect(names(mean_case_global), rownames(avg_stats))
      avg_stats <- avg_stats[valid_genes_final, , drop = FALSE]

      group_logfc <- abs(avg_stats$mean_global - mean_case_global[valid_genes_final])
      avg_stats$group_logfc <- group_logfc

      pass_group <- is.finite(group_logfc) & (group_logfc <= max_group_logfc_cut)
      avg_stats <- avg_stats[pass_group, , drop = FALSE]
      common_genes <- rownames(avg_stats)

      if (length(common_genes) == 0) {
        warning(sprintf("All candidate genes failed the Case vs Control logFC condition (LogFC > %.2f).", max_group_logfc_cut), call. = FALSE)
        object@params_expression$stable_genes <- character(0)
        return(object)
      }
    } else {
      if (verbose) message(">>> Group Stability: Case samples not found (or N=1). Skipping logFC check.")
      avg_stats$group_logfc <- NA_real_
    }
  } else {
    avg_stats$group_logfc <- NA_real_
  }

  # 4. Final sorting and ranking
  has_group_logfc <- all(!is.na(avg_stats$group_logfc))

  if (check_group_diff && has_group_logfc) {
    avg_stats <- avg_stats[order(avg_stats$group_logfc, avg_stats$cv_robust, avg_stats$max_abs_logfc), , drop = FALSE]
  } else {
    avg_stats <- avg_stats[order(avg_stats$cv_robust, avg_stats$max_abs_logfc, -avg_stats$mean_global), , drop = FALSE]
  }

  final_genes <- utils::head(avg_stats$gene, max_num_gene)

  if (verbose) {
    message(sprintf(">>> Integration Complete: Selected top %d stable genes.", length(final_genes)))
  }

  object@params_expression$stable_genes <- final_genes
  object@params_expression$stable_gene_stats <- avg_stats[avg_stats$gene %in% final_genes, ]
  return(object)
}

#' @title Create Custom Gene Sets for Simulation and Testing
#'
#' @description
#' Constructs customized gene sets (SVG, Null, and Spike) based on a specified
#' target and reference domain within a selected pilot sample. Identifies the
#' top differentially expressed genes from a pool of stable candidate genes.
#'
#' @param object A \code{spaCraft} object.
#' @param G_svg_base Character vector. Base marker genes known a priori.
#' @param G_stable Character vector. Stable genes identified via feature selection.
#' @param target_domain Character. The name of the target domain (e.g., "WM").
#' @param reference_domain Character. The name of the reference domain (e.g., "Layer6").
#' @param sample_index Integer or Character. The index or name of the pilot sample to use (default: 1L).
#' @param n_de Integer. Number of top differentially expressed genes to extract (default: 5L).
#' @param verbose Logical. Print progress and summary messages (default: TRUE).
#'
#' @return A list containing:
#'   \item{G_svg}{Character vector. Union of \code{G_svg_base} and \code{G_stable}.}
#'   \item{G_null}{Character vector. Unique stable genes.}
#'   \item{G_spike}{Character vector. Top \code{n_de} differentially expressed genes.}
#'   \item{TR_stats}{Data frame. Mean expressions and differences for all tested genes.}
#' @importFrom utils head
#' @export
#'
makeCustomGeneSets <- function(object,
                               G_svg_base,
                               G_stable,
                               target_domain,
                               reference_domain,
                               sample_index = 1L,
                               n_de = 5L,
                               verbose = TRUE) {

  pilot_data <- object@pilot_data
  if (is.null(pilot_data) || length(pilot_data) == 0) {
    stop("The spaCraft object contains no pilot_data.")
  }

  if (is.numeric(sample_index)) {
    if (sample_index < 1 || sample_index > length(pilot_data)) {
      stop(sprintf("sample_index %d is out of bounds (1 to %d).", sample_index, length(pilot_data)))
    }
  } else if (is.character(sample_index)) {
    if (!(sample_index %in% names(pilot_data))) {
      stop(sprintf("sample_index '%s' not found in pilot_data.", sample_index))
    }
  } else {
    stop("sample_index must be an integer or character string.")
  }

  samp   <- pilot_data[[sample_index]]
  coords <- samp$coords

  # Expression source. Keep it sparse for now so the row subset below is cheap;
  # densify only the handful of genes we actually score. log1p is applied only
  # when no precomputed logcounts are available.
  if (!is.null(samp$logcounts)) {
    X        <- samp$logcounts
    take_log <- FALSE
  } else if (!is.null(samp$counts)) {
    X        <- samp$counts
    take_log <- TRUE
  } else {
    stop("Selected sample does not contain 'counts' or 'logcounts'.")
  }

  idx_T <- which(as.character(coords$domain) == target_domain)
  idx_R <- which(as.character(coords$domain) == reference_domain)

  if (length(idx_T) == 0) {
    stop(sprintf("Target domain '%s' not found in sample %s.", target_domain, as.character(sample_index)))
  }
  if (length(idx_R) == 0) {
    stop(sprintf("Reference domain '%s' not found in sample %s.", reference_domain, as.character(sample_index)))
  }

  G_svg  <- unique(c(G_svg_base, G_stable))
  G_null <- unique(G_stable)
  valid_genes <- intersect(G_null, rownames(X))

  if (length(valid_genes) == 0) {
    stop("None of the genes in G_stable are present in the sample's expression matrix.")
  }

  # Subset to the candidate genes on the (sparse) matrix, then densify that
  # small submatrix so base rowMeans() works regardless of storage class
  # (dgCMatrix is not a base array, which is what caused the original error).
  Xv <- as.matrix(X[valid_genes, , drop = FALSE])
  if (take_log) Xv <- log1p(Xv)

  mean_T  <- rowMeans(Xv[, idx_T, drop = FALSE], na.rm = TRUE)
  mean_R  <- rowMeans(Xv[, idx_R, drop = FALSE], na.rm = TRUE)

  TR_diff <- mean_T - mean_R
  abs_TR  <- abs(TR_diff)

  st_df <- data.frame(
    gene   = valid_genes,
    mean_T = mean_T,
    mean_R = mean_R,
    TR     = TR_diff,
    abs_TR = abs_TR,
    stringsAsFactors = FALSE
  )

  st_sorted <- st_df[order(st_df$abs_TR, decreasing = TRUE), , drop = FALSE]
  G_spike   <- utils::head(st_sorted$gene, n = n_de)

  if (verbose) {
    message(">>> Custom Gene Set Summary:")
    message(sprintf("    G_svg (Base + Stable) : %d genes", length(G_svg)))
    message(sprintf("    G_null (Stable only)  : %d genes", length(G_null)))
    message(sprintf("    G_spike (Top DE)      : %d genes extracted", length(G_spike)))
  }

  list(
    G_svg    = G_svg,
    G_null   = G_null,
    G_spike  = G_spike,
    TR_stats = st_sorted
  )
}
