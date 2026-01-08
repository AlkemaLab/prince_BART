#' Principal Stratification using BART
#'
#' Fits a Bayesian principal stratification model using Bayesian Additive
#' Regression Trees (BART) for causal inference with endogenous treatments.
#' The method is designed for instrumental variable or encouragement designs
#' with noncompliance and targets causal effects for compliers.
#'
#' @param formula A formula following 2SLS convention:
#'   \code{Y ~ X1 + X2 | Z | W}, where Y is the outcome, X1 + X2 are covariates,
#'   Z is the instrument (treatment assignment), and W is treatment uptake.
#' @param data A data.frame containing the variables in the formula.
#' @param X A matrix or data.frame of covariates. Used when \code{formula = NULL}.
#' @param Y A binary outcome vector (0/1). Used when \code{formula = NULL}.
#' @param Z A binary instrument or treatment assignment vector (0/1).
#' @param W A binary treatment uptake/received vector (0/1).
#' @param propensity Optional pre-computed instrument propensity scores
#'   \eqn{e = P(Z \mid X)}. If NULL (default), propensity scores are estimated
#'   internally using BART.
#' @param instrument_overlap Optional bounds for instrument propensity trimming
#'   to enforce overlap. If provided as a length-2 vector (e.g., \code{c(0.1, 0.9)}),
#'   observations with estimated propensity scores outside this range are excluded.
#'   Default is NULL (no trimming). See Crump et al. (2009).
#' @param n_chains Number of parallel MCMC chains (default: 4).
#' @param n_warmup Number of warmup iterations per chain (default: 1000).
#' @param n_samples Number of posterior samples per chain (default: 1000).
#' @param keep_trees Logical; save fitted BART tree structures for downstream
#'   prediction and generalization (default: FALSE).
#' @param k Prior hyperparameter controlling node shrinkage in BART (default: 2).
#' @param n_trees Number of trees in each BART ensemble (default: 200).
#' @param workers Number of parallel workers. If NULL (default), uses all
#'   available cores up to \code{n_chains}. Set to 1 for sequential execution.
#' @param verbose Logical; print progress messages (default: FALSE).
#'
#' @details
#' The model jointly estimates:
#' \itemize{
#'   \item Principal stratum membership probabilities
#'     (\eqn{P(\text{complier} \mid X)}, \eqn{P(\text{never-taker} \mid X)},
#'      \eqn{P(\text{always-taker} \mid X)}),
#'   \item Potential outcome regressions within strata, including
#'     \eqn{E[Y(0) \mid \text{complier}, X]} and
#'     \eqn{E[Y(1) \mid \text{complier}, X]}.
#' }
#'
#' The primary identified causal estimand is the conditional average treatment
#' effect among compliers, \eqn{\mathrm{CATE}_C(x)}. Mixed (sample-based) averages
#' of these conditional effects, including LATE-like estimands, can be obtained
#' using \code{summary()} and \code{coef()}.
#'
#' Identification relies on standard instrumental variable assumptions:
#' conditional independence of the instrument given X, exclusion restriction,
#' and monotonicity (no defiers).
#'
#' Parallel computation is handled via
#' \code{future.apply::future_lapply}. Before calling this function, set your
#' preferred parallel backend:
#'
#' \preformatted{
#' # For Unix/Mac (forked processes)
#' future::plan(future::multicore)
#'
#' # For Windows or any platform
#' future::plan(future::multisession)
#' }
#'
#' @return A list of class \code{"princebart"} containing posterior draws from
#'   all chains, including:
#'   \itemize{
#'     \item \code{imp}: Imputed principal stratum indicators
#'       (iteration x chain x stratum x unit).
#'     \item \code{probs}: Posterior draws of stratum probabilities and
#'       stratum-specific outcome means, including
#'       \code{p_a}, \code{p_n}, \code{m_y0c}, \code{m_y1c}, etc.
#'     \item \code{trees}: Fitted BART trees (if \code{keep_trees = TRUE}).
#'     \item \code{data}: Processed input data, including covariates,
#'       instrument propensity scores, and outcomes.
#'   }
#'
#' @examples
#' \dontrun{
#' library(future)
#' plan(multisession, workers = 4)
#'
#' fit <- prince_BART(
#'   Y ~ X1 + X2 + X3 | Z | W,
#'   data = mydata,
#'   n_chains = 4,
#'   n_warmup = 1000,
#'   n_samples = 1500
#' )
#'
#' summary(fit)
#' }
#'
#' @seealso \code{\link{segment_heterogeneity}},
#'   \code{\link{general_BART}}
#'
#' @export
prince_BART <- function(
  formula = NULL,
  data = NULL,
  X = NULL,
  Y = NULL,
  Z = NULL,

  W = NULL,
  propensity = NULL,
  instrument_overlap = NULL,
  n_warmup = 1000L,
  n_samples = 1000L,
  n_chains = 4L,
  keep_trees = FALSE,
  k = 2,
  n_trees = 200L,
  workers = NULL,
  verbose = FALSE
) {

  # Hardcode n_initial = 0 (MoM offsets not clearly helpful)
  n_initial <- 0L

  # Handle formula interface
  if (!is.null(formula)) {
    parsed <- parse_psbart_formula(formula, data)
    X <- parsed$X
    Y <- parsed$Y
    Z <- parsed$Z
    W <- parsed$W
  }

  # Input validation
  if (is.null(X) || is.null(Y) || is.null(Z) || is.null(W)) {
    stop("Must provide either a formula + data, or X, Y, Z, W directly")
  }

  # Validate and prepare data
  X <- validate_and_prepare_X(X)
  Y <- validate_binary(Y, "Y")
  Z <- validate_binary(Z, "Z")
  W <- validate_binary(W, "W")

  n <- nrow(X)
  if (length(Y) != n || length(Z) != n || length(W) != n) {
    stop("X, Y, Z, W must all have the same number of observations")
  }

  # Store scaling attributes before adding e
  scaled_center <- attr(X, "scaled:center")
  scaled_scale <- attr(X, "scaled:scale")

  # Compute propensity scores
  if (is.null(propensity)) {
    if (verbose) message("Computing propensity scores...")
    e <- dbarts::bart2(X, Z, verbose = FALSE) |> stats::fitted() |> stats::qnorm()
  } else {
    validate_propensity(propensity, n)
    e <- stats::qnorm(propensity)
  }

  # Enforce instrument overlap by trimming extreme propensity scores
  if (!is.null(instrument_overlap)) {
    if (length(instrument_overlap) != 2) {
      stop("instrument_overlap must be a length-2 vector, e.g., c(0.1, 0.9)")
    }
    e_prob <- stats::pnorm(e)
    keep <- e_prob >= instrument_overlap[1] & e_prob <= instrument_overlap[2]
    n_trimmed <- sum(!keep)
    if (verbose) {
      message("Trimming ", n_trimmed, " observations (",
              round(100 * n_trimmed / n, 1), "%) outside propensity range [",
              instrument_overlap[1], ", ", instrument_overlap[2], "]")
    }
    if (sum(keep) < 10) {
      stop("Too few observations remain after instrument overlap trimming")
    }
    X <- X[keep, , drop = FALSE]
    Y <- Y[keep]
    Z <- Z[keep]
    W <- W[keep]
    e <- e[keep]
    n <- sum(keep)
  }

  # Append propensity to X
  X <- cbind(X, e = e)
  
  # Preserve scaling attributes on X (for general_BART)
  attr(X, "scaled:center") <- scaled_center
  attr(X, "scaled:scale") <- scaled_scale

  # Determine number of workers
  if (is.null(workers)) {
    workers <- min(parallel::detectCores(), n_chains)
  }

  # Set up future plan if not already configured
  old_plan <- future::plan()
  if (inherits(old_plan, "sequential") && workers > 1) {
    if (verbose) message("Setting up multisession plan with ", workers, " workers")
    future::plan(future::multisession, workers = workers)
    on.exit(future::plan(old_plan), add = TRUE)
  }

  if (verbose) message("Running ", n_chains, " chains...")

  # Run chains in parallel
  res0 <- future.apply::future_lapply(
    seq_len(n_chains),
    function(chain_id) {
      .fit_psbart(
        X = X,
        Y = Y,
        Z = Z,
        W = W,
        n_warmup = n_warmup,
        n_samples = n_samples,
        save_trees = keep_trees,
        k = k,
        n_trees = n_trees,
        n_initial = n_initial,
        verbose = FALSE
      )
    },
    future.seed = TRUE,
    future.packages = "princeBART"
  )

  # Combine results
  res <- list()

  # Combine trees
  if (keep_trees) {
    list_tree <- lapply(res0, function(x) x$trees)
    res$trees <- do.call(rbind, Map(function(df, id) {
      df$chain <- id
      df
    }, list_tree, seq_along(list_tree)))
  } else {
    res$trees <- NULL
  }

  # Combine imputations
  list_imp <- lapply(res0, function(x) x$imputed)
  res$imp <- abind::abind(list_imp, along = 4)
  res$imp <- aperm(res$imp, c(1, 4, 3, 2))
  dimnames(res$imp) <- list(
    iteration = NULL,
    chain = NULL,
    variable = c("nt", "at"),
    unit = NULL
  )

  # Combine probabilities
  # .fit_psbart outputs in order: p_at, p_nt, m_y0co, m_y1co, m_y0nt, m_y1at
  list_probs <- lapply(res0, function(x) x$probs)
  res$probs <- abind::abind(list_probs, along = 4)
  res$probs <- aperm(res$probs, c(1, 4, 3, 2))
  dimnames(res$probs) <- list(
    iteration = NULL,
    chain = NULL,
    variable = c("p_a", "p_n", "m_y0c", "m_y1c", "m_y0n", "m_y1a"),
    unit = NULL
  )

  # Store data reference with UNSCALED X (for general_BART compatibility)
  # X is currently scaled; unscale before storing
  X_unscaled <- X[, -ncol(X), drop = FALSE]  # Remove e column first
  if (!is.null(scaled_center) && !is.null(scaled_scale)) {
    # Reverse scaling: X_original = X_scaled * scale + center
    for (j in seq_len(ncol(X_unscaled))) {
      v <- colnames(X_unscaled)[j]
      if (v %in% names(scaled_center) && v %in% names(scaled_scale)) {
        X_unscaled[, j] <- X_unscaled[, j] * scaled_scale[v] + scaled_center[v]
      }
    }
  }
  # Add back e (propensity is not scaled)
  X_unscaled <- cbind(X_unscaled, e = X[, "e"])
  
  res$data <- list(X = X_unscaled, Y = Y, Z = Z, W = W, e = X[, "e"])
  
  # Store scaling parameters separately for use in predictions
  res$scaling <- list(
    center = scaled_center,
    scale = scaled_scale
  )
  
  res$call <- match.call()

  class(res) <- "princebart"
  
  if (verbose) message("Done.")
  res
}
