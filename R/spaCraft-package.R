#'spaCraft: Analysis-aware power calculation and
#'sample size planning for multi-sample spatial transcriptomics
#'
#' @description
#' A statistical framework for the prospective design of comparative spatial
#' transcriptomics (ST) studies. \pkg{spaCraft} learns a cohort-level generative
#' model from pilot data to estimate statistical power under a full
#' \emph{generate--recover--test} pipeline. By explicitly modeling tissue
#' architecture, spatial dependence, and domain-labeling uncertainty, it treats
#' computational domain recovery as an integral part of the study design.
#'
#' @details
#' The framework formulates sample-size planning as a pilot-based cohort design
#' problem. The generative model consists of three pilot-fitted layers:
#' \itemize{
#'   \item \strong{Domain composition:} A Binomial-logit model for Target/Reference
#'         proportions, summarizing case-control contrasts via odds ratios.
#'   \item \strong{Spatial geometry:} A Fisher-Gaussian kernel mixture model (FGKMM)
#'         capturing complex, non-convex domain shapes and centroid placements.
#'   \item \strong{Spatial gene expression:} An additive model combining domain means,
#'         global shifts, and a Gaussian process (GP) residual field.
#' }
#'
#' Synthetic multi-sample cohorts are generated under user-specified effect sizes.
#' Domain labels are then recovered using \strong{pilot-guided BANKSY (pBANKSY)},
#' and power is evaluated for two primary endpoints:
#' \itemize{
#'   \item \strong{SaLFC}: A spatially adjusted, empirical-Bayes moderated log-fold-change
#'         test for differential expression.
#'   \item \strong{LOR}: A baseline-anchored log-odds-ratio test for compositional shifts.
#' }
#' Finally, shape-constrained additive models (SCAM) are used to smooth Monte-Carlo
#' power estimates and recommend optimal per-group sample sizes.
#'
#' @section Typical Workflow:
#' A typical study design chains the following functions:
#' \enumerate{
#'   \item \code{\link{createspaCraftObject}}: Initialize the S4 design object from pilot data.
#'   \item \code{\link{featureSelection}} & \code{\link{makeCustomGeneSets}}: Define
#'         labeling, empirical-null, and DE-injection gene sets.
#'   \item \strong{Estimation}: Fit the generative layers using \code{\link{estimateCompositionParams}},
#'         \code{\link{estimateGeometryParams}}, and \code{\link{estimateExpressionParams}}.
#'   \item \code{\link{simulatespaCraft}}: Generate synthetic multi-sample cohorts.
#'   \item \code{\link{pBANSKY}}: Recover spatial domains on synthetic samples.
#'   \item \strong{Endpoint Testing}: Evaluate endpoints via \code{\link{SaLFC}} or \code{\link{LOR}}.
#'   \item \strong{Power Evaluation}: Compute power grids and fit SCAM curves using
#'         \code{\link{evaluatePowerSaLFC}} or \code{\link{evaluatePowerLOR}}.
#' }
#'
#' @keywords internal
"_PACKAGE"
#'
#' @aliases spaCraft-package spaCraft
#' @author Jungmin Shin, Dongjun Chung \cr
#'   Maintainer: Jungmin Shin \email{c16267@@gmail.com}
#'
#' @references
#' Mukhopadhyay, M., Li, D. and Dunson, D. B. (2020).
#' \emph{Estimating densities with non-linear support by using Fisher-Gaussian
#' kernels, JRSS-B 82(5): 1249--1271}. \cr
#' Saha, A. and Datta, A. (2018).
#' \emph{BRISC: bootstrap for rapid inference on spatial covariances, Stat 7: e184}. \cr
#' Singhal, V. et al. (2024).
#' \emph{BANKSY unifies cell typing and tissue domain segmentation for scalable
#' spatial omics data analysis, Nature Genetics 56: 431--441}. \cr
#' Smyth, G. K. (2004).
#' \emph{Linear models and empirical Bayes methods for assessing differential
#' expression in microarray experiments, SAGMB 3(1): Article 3}. \cr
#' McCarthy, D. J. and Smyth, G. K. (2009).
#' \emph{Testing significance relative to a fold-change threshold is a TREAT,
#' Bioinformatics 25(6): 765--771}. \cr
#' Pya, N. and Wood, S. N. (2015).
#' \emph{Shape constrained additive models, Statistics and Computing 25(3): 543--559}. \cr
#'
#' @import methods stats
#' @importFrom MASS mvrnorm
#' @importFrom BRISC BRISC_estimation
#' @importFrom FNN get.knn get.knnx knn.index
#' @importFrom limma eBayes squeezeVar
#' @importFrom scam scam
#' @importFrom clue solve_LSAP
#' @importFrom Rcpp sourceCpp
#' @useDynLib spaCraft, .registration = TRUE
NULL
