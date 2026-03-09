#' @name prince_bart-class
#' @title S3 Methods for prince_bart Objects
#'
#' @description Print, summary, and coef methods for objects returned
#' by \code{prince_BART}. Behavior is specialized by modality: binary fits
#' report principal strata and treatment effects for compliers, while ordinal
#' fits report contrasts for the monotone compliance group and stratified by
#' baseline uptake level.
#'
#' @param x,object A prince_bart object (either prince_bart_binary or prince_bart_ordinal).
#' @param type Character; type of estimand: \code{"mixed"} (default) for mixed
#'   estimands averaging over principal stratum uncertainty, or \code{"sample"}
#'   for sample estimands using imputed principal strata (experimental).
#'   For ordinal fits, this parameter is accepted but may be ignored depending
#'   on implementation.
#' @param treated_only Logical; if \code{TRUE}, compute treatment effect only
#'   among the treated units (ATT). Default is \code{FALSE} (ATE). For ordinal
#'   fits, this parameter is accepted but may be ignored or handled differently.
#' @param ... Additional arguments (currently ignored).
#'
#' @return
#' \code{print}: Invisibly returns the object.
#' \code{summary}: Prints and invisibly returns a list with posterior summaries
#'   appropriate to the modality.
#' \code{coef}: Posterior summary of the requested estimand.
#'
#' @examples
#' \dontrun{
#' # Binary fit
#' fit_binary <- prince_BART(Y ~ X | Z | W, data = df)
#' summary(fit_binary)
#' coef(fit_binary)
#'
#' # Ordinal fit
#' fit_ordinal <- prince_BART(Y ~ X | Z | W, data = df, uptake_type = "ordinal")
#' summary(fit_ordinal)
#' coef(fit_ordinal)
#' }
NULL


#' @rdname prince_bart-class
#' @export
print.prince_bart <- function(x, ...) {
  cat("Principal Stratification BART Fit\n")
  cat("---------------------------------\n")

  if (!is.null(dim(x$probs))) {
    dims <- dim(x$probs)
    if (length(dims) == 4) {
      cat("Chains:      ", dims[2], "\n")
      cat("Iterations:  ", dims[1], "\n")
      cat("Units:       ", dims[4], "\n")
    } else if (length(dims) == 3) {
      cat("Iterations:  ", dims[1], "\n")
      cat("Units:       ", dims[2], "\n")
    }
  }

  if (!is.null(x$trees)) {
    cat("Trees:        saved\n")
  }

  uptake_label <- if (inherits(x, "prince_bart_binary")) {
    "Binary uptake (3 principal strata: complier, never-taker, always-taker)"
  } else {
    "Ordinal/count uptake"
  }
  cat("Uptake model: ", uptake_label, "\n")

  cat("\nUse summary() or coef() to extract contrasts and treatment effects.\n")
  invisible(x)
}


#' @rdname prince_bart-class
#' @export
summary.prince_bart <- function(object
  , type = c("mixed", "sample")
  , treated_only = FALSE
  , ...
) {
  if (inherits(object, "prince_bart_binary")) {
    summary_prince_bart_binary(object, type = type, treated_only = treated_only)
  } else if (inherits(object, "prince_bart_ordinal")) {
    summary_prince_bart_ordinal(object, type = type, treated_only = treated_only)
  } else {
    stop("Unknown prince_bart class: ", paste(class(object), collapse = ", "))
  }
}


#' @rdname prince_bart-class
#' @export
coef.prince_bart <- function(object
  , type = c("mixed", "sample")
  , treated_only = FALSE
  , ...
) {
  if (inherits(object, "prince_bart_binary")) {
    coef_prince_bart_binary(object, type = type, treated_only = treated_only)
  } else if (inherits(object, "prince_bart_ordinal")) {
    coef_prince_bart_ordinal(object, type = type, treated_only = treated_only)
  } else {
    stop("Unknown prince_bart class: ", paste(class(object), collapse = ", "))
  }
}


# =============================================================================
# Internal helper functions for modality-specific summary/coef
# =============================================================================

#' @keywords internal
summary_prince_bart_binary <- function(object, type = c("mixed", "sample"), treated_only = FALSE) {
  type <- match.arg(type)

  cat("Principal Stratification BART Summary (Binary Uptake)\n")
  cat("=====================================================\n\n")

  # Compute treatment effect using existing binary estimand functions
  effect_label <- if (treated_only) "ATT" else "ATE"
  type_label <- if (type == "mixed") "Mixed" else "Sample"

  if (type == "mixed") {
    if (treated_only) {
      effect <- matt_c(object)
    } else {
      effect <- mate_c(object)
    }
  } else {
    if (treated_only) {
      effect <- satt_c(object)
    } else {
      effect <- sate_c(object)
    }
  }

  # effect[[1]] = strata distribution, effect[[2]] = outcomes & treatment effect
  strata_summary <- effect[[1]]
  effect_summary <- effect[[2]]

  cat("Principal Strata Distribution:\n")
  print(as.data.frame(strata_summary))
  cat("\n")

  cat(sprintf("%s %s for Compliers:\n", type_label, effect_label))
  print(as.data.frame(effect_summary))

  invisible(list(
    strata = strata_summary,
    effect = effect_summary,
    type = type,
    treated_only = treated_only
  ))
}

#' @keywords internal
coef_prince_bart_binary <- function(object, type = c("mixed", "sample"), treated_only = FALSE) {
  type <- match.arg(type)

  if (type == "mixed") {
    if (treated_only) {
      result <- matt_c(object)[[2]]
    } else {
      result <- mate_c(object)[[2]]
    }
  } else {
    if (treated_only) {
      result <- satt_c(object)[[2]]
    } else {
      result <- sate_c(object)[[2]]
    }
  }

  result
}

#' @keywords internal
summary_prince_bart_ordinal <- function(object, type = c("mixed", "sample"), treated_only = FALSE) {
  cat("Principal Stratification BART Summary (Ordinal Uptake)\n")
  cat("=====================================================\n\n")
  cat("Type and treated_only parameters are currently interpreted in overall contrast computation.\n\n")

  # Delegate to ordinal estimand functions (will be implemented in estimands.R)
  results <- estimands_ordinal_mixed(object)

  cat("Complier-Like Contrast [W(0) - W(1) = 1]:\n")
  print(as.data.frame(results$overall))
  cat("\n")

  cat("Contrasts by Baseline Uptake W(0):\n")
  print(as.data.frame(results$by_w0))

  if (any(is.na(results$overall$mean)) || any(is.na(results$by_w0$mean))) {
    cat("\nNote: NA rows indicate no posterior units matched that contrast in some draws/strata.\n")
  }

  invisible(list(
    overall = results$overall,
    by_w0 = results$by_w0
  ))
}

#' @keywords internal
coef_prince_bart_ordinal <- function(object, type = c("mixed", "sample"), treated_only = FALSE) {
  results <- estimands_ordinal_mixed(object)
  # Return combined draws for posterior summary
  c(results$overall, results$by_w0)
}
