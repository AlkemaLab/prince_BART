#' princeBART: Principal Stratification with BART
#'
#' @description
#' The princeBART package implements principal stratification with Bayesian
#' Additive Regression Trees (BART) for causal inference when treatment uptake
#' is endogenous and identification relies on an instrument that is as-good-as-
#' random conditional on covariates. It targets the average causal effect among
#' compliers (ATE_C, also known as the Local Average Treatment Effect, LATE),
#' as well as conditional complier effects given covariates (CATE_C(x)).
#'
#' @section Main Functions:
#' \itemize{
#'   \item \code{\link{prince_BART}}: Fit the principal stratification BART model
#'   \item \code{\link{segment_heterogeneity}}: Summarize heterogeneous effects by segments
#'   \item \code{\link{general_BART}}: Generalize effects to external populations (PATE)
#' }
#'
#' @section Model:
#' The model applies to settings with an instrument Z, an endogenous treatment W,
#' and an outcome Y. It requires conditional randomization of Z given covariates X,
#' and allows for both one-sided and two-sided noncompliance.
#' Three principal strata are identified:
#' \itemize{
#'   \item \strong{Compliers}: W(1) = 1, W(0) = 0 — treatment responds to instrument
#'   \item \strong{Never-takers}: W(1) = W(0) = 0 — never treated regardless of Z
#'   \item \strong{Always-takers}: W(1) = W(0) = 1 — always treated regardless of Z
#' }
#'
#' BART is used to flexibly model both the stratum membership probabilities
#' and the potential outcomes within each stratum, allowing for heterogeneous
#' conditional complier effects CATE_C(x) across covariate values.
#'
#' @docType package
#' @name princeBART-package
#' @aliases princeBART
#'
#' @importFrom dbarts dbarts bart2 dbartsControl dbartsData
#' @importFrom posterior as_draws_array
#' @importFrom abind abind
#' @importFrom future.apply future_lapply
#' @importFrom stats pnorm qnorm rbinom lm coef fitted model.frame model.matrix model.response formula
#' @importFrom methods new
#' @importFrom parallel detectCores
#' @importFrom future plan multisession
#' @importFrom Formula Formula
#' @importFrom utils getFromNamespace
#' @importFrom data.table := .N .SD
"_PACKAGE"

# data.table checks this flag for NSE operations (:=, .N, .SD) inside package code.
.datatable.aware <- TRUE
