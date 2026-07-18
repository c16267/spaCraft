#' @title spaCraft Class
#' @description
#' S4 class for multi-sample spatial transcriptomics power analysis. Stores the
#' pilot data, the hierarchical model parameters fitted from it, the simulated
#' synthetic data, and the downstream testing results.
#'
#' @slot pilot_data list. Pilot samples; each element contains \code{counts}
#'   (genes x spots), \code{coords} (x, y, domain), \code{group} (0/1), and
#'   \code{sample_id}.
#' @slot params_geometry list. Hierarchical FGKMM parameters (centroids, shapes,
#'   pooling info).
#' @slot params_expression list. Hierarchical NNGP parameters (mean, covariance,
#'   sigma_bio).
#' @slot params_composition list. Beta-Binomial parameters for domain proportions.
#' @slot synthetic_data list. Simulated/recovered samples produced by the
#'   simulation and clustering pipeline. Each element holds \code{counts},
#'   \code{coords}, and \code{meta}.
#' @slot testing_result list. Endpoint testing results, keyed by endpoint
#'   (e.g., \code{salfc}, \code{composition}).
#'
#' @exportClass spaCraft
setClass("spaCraft",
         slots = c(
           pilot_data         = "list",
           params_geometry    = "list",
           params_expression  = "list",
           params_composition = "list",
           synthetic_data     = "list",
           testing_result     = "list"
         ))

# Internal null-coalescing helper
#' @keywords internal
#' @noRd
`%||%` <- function(a, b) if (is.null(a)) b else a

#' @title Create spaCraft Object
#' @description
#' Creates a \code{spaCraft} object from a pre-processed list of pilot data.
#' If the dataset contains only single replicates (K=1) globally or per group,
#' this function can automatically generate pseudo-replicates via minor random
#' perturbations to ensure downstream multi-sample estimators function correctly.
#'
#' @param pilot_list A list of processed samples. Each element must contain
#'   \code{counts} (genes x spots), \code{coords} (with columns 'x', 'y',
#'   'domain'), and \code{group}.
#' @param generate_pseudo_reps Logical. If TRUE, expands any group with K=1 to
#'   K=3 by adding minor spatial and count perturbations. Default is TRUE.
#' @param noise_level Numeric. Fractional uniform noise added to counts for
#'   pseudo-replicates. Default is 0.05 (+/- 5%).
#'
#' @return A \code{spaCraft} object ready for parameter estimation.
#'
#' @importFrom stats runif
#' @importFrom methods new
#' @export
#'
createspaCraftObject <- function(pilot_list,
                                   generate_pseudo_reps = TRUE,
                                   noise_level = 0.05) {

  if (!is.list(pilot_list) || length(pilot_list) == 0) {
    stop("Input must be a non-empty list of pilot samples.")
  }

  required_fields <- c("counts", "coords", "group")
  for (i in seq_along(pilot_list)) {
    missing_fields <- setdiff(required_fields, names(pilot_list[[i]]))
    if (length(missing_fields) > 0) {
      stop(sprintf("Sample %d is missing required fields: %s",
                   i, paste(missing_fields, collapse = ", ")))
    }
  }

  message(">>> Initializing spaCraft object...")

  # Handle single pilot dataset (no groups) by duplicating into control/case
  if (length(pilot_list) == 1) {
    message("   ! Only 1 total pilot sample detected.")
    message("   -> Duplicating to create '0' (Control) and '1' (Case) groups...")

    s_ctrl <- pilot_list[[1]]
    s_ctrl$group     <- 0
    s_ctrl$sample_id <- paste0(s_ctrl$sample_id %||% "sample1", "_ctrl")

    s_case <- pilot_list[[1]]
    s_case$group     <- 1
    s_case$sample_id <- paste0(s_case$sample_id %||% "sample1", "_case")

    pilot_list <- list(s_ctrl, s_case)
  }

  # Handle pseudo-replication for K=1 groups
  if (generate_pseudo_reps) {
    groups     <- vapply(pilot_list, function(x) as.character(x$group), character(1))
    grp_counts <- table(groups)

    new_pilot_list <- list()

    for (grp in names(grp_counts)) {
      grp_samples <- pilot_list[groups == grp]
      K <- length(grp_samples)

      if (K == 1) {
        warning(sprintf(
          "Group '%s' has K=1. Generating 2 pseudo-replicates to force K=3.\n*** CAUTION: Pseudo-replicates use artificial noise (+/- %d%%) and do not replace biological variance. Power analyses will be highly dependent on the `noise_level` parameter. ***",
          grp, round(noise_level * 100)), call. = FALSE)

        orig <- grp_samples[[1]]
        new_pilot_list[[length(new_pilot_list) + 1]] <- orig

        for (rep_idx in 1:2) {
          pseudo <- orig
          pseudo$sample_id <- paste0(orig$sample_id %||% "sample", "_pseudo_", rep_idx)

          # Perturb coordinates
          n_spots <- nrow(pseudo$coords)
          pseudo$coords$x <- pseudo$coords$x + stats::runif(n_spots, -1e-4, 1e-4)
          pseudo$coords$y <- pseudo$coords$y + stats::runif(n_spots, -1e-4, 1e-4)

          # Perturb counts
          noise_mat <- matrix(
            stats::runif(length(pseudo$counts), 1 - noise_level, 1 + noise_level),
            nrow = nrow(pseudo$counts), ncol = ncol(pseudo$counts))
          pseudo$counts <- round(pseudo$counts * noise_mat)

          new_pilot_list[[length(new_pilot_list) + 1]] <- pseudo
        }
      } else {
        new_pilot_list <- c(new_pilot_list, grp_samples)
      }
    }
    pilot_list <- new_pilot_list
  }

  new_obj <- methods::new("spaCraft",
                          pilot_data         = pilot_list,
                          params_geometry    = list(),
                          params_expression  = list(),
                          params_composition = list(),
                          synthetic_data     = list(),
                          testing_result     = list())

  message(sprintf("Successfully created spaCraft object with %d total samples.",
                  length(pilot_list)))

  return(new_obj)
}

#' @title Get or Set Synthetic Data
#' @description Accessor and replacement methods for the \code{synthetic_data}
#'   slot of a \code{spaCraft} object.
#' @param object A \code{spaCraft} object.
#' @param value A named list of simulated samples to store.
#' @return \code{syntheticData} returns the stored list; the replacement method
#'   returns the updated \code{spaCraft} object.
#' @export
setGeneric("syntheticData", function(object) standardGeneric("syntheticData"))

#' @rdname syntheticData
#' @export
setMethod("syntheticData", "spaCraft", function(object) object@synthetic_data)

#' @rdname syntheticData
#' @export
setGeneric("syntheticData<-", function(object, value) standardGeneric("syntheticData<-"))

#' @rdname syntheticData
#' @importFrom methods validObject
#' @export
setReplaceMethod("syntheticData", "spaCraft", function(object, value) {
  if (!is.list(value)) stop("synthetic_data must be a list.")
  object@synthetic_data <- value
  methods::validObject(object)
  return(object)
})

#' @title Get or Set Testing Results
#' @description Accessor and replacement methods for the \code{testing_result}
#'   slot of a \code{spaCraft} object.
#' @param object A \code{spaCraft} object.
#' @param value A named list of endpoint results to store.
#' @return \code{testingResult} returns the stored list; the replacement method
#'   returns the updated \code{spaCraft} object.
#' @export
setGeneric("testingResult", function(object) standardGeneric("testingResult"))

#' @rdname testingResult
#' @export
setMethod("testingResult", "spaCraft", function(object) object@testing_result)

#' @rdname testingResult
#' @export
setGeneric("testingResult<-", function(object, value) standardGeneric("testingResult<-"))

#' @rdname testingResult
#' @importFrom methods validObject
#' @export
setReplaceMethod("testingResult", "spaCraft", function(object, value) {
  if (!is.list(value)) stop("testing_result must be a list.")
  object@testing_result <- value
  methods::validObject(object)
  return(object)
})

#' @title Add a Single Endpoint Result
#' @description Convenience helper that inserts or overwrites one named endpoint in
#'   the \code{testing_result} slot, leaving existing endpoints intact.
#' @param object A \code{spaCraft} object.
#' @param name Character. Endpoint key (e.g., "salfc", "composition").
#' @param result The result object to store under \code{name}.
#' @return The updated \code{spaCraft} object.
#' @importFrom methods is validObject
#' @export
addTestingResult <- function(object, name, result) {
  stopifnot(methods::is(object, "spaCraft"),
            is.character(name), length(name) == 1L)

  tr <- object@testing_result
  if (!is.list(tr)) tr <- list()

  tr[[name]] <- result
  object@testing_result <- tr
  methods::validObject(object)

  return(object)
}
