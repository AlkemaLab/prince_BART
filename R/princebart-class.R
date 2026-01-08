#' @name princebart-class
#' @title S3 Methods for princebart Objects
#'
#' @description Print, summary, and coef methods for objects returned
#' by \code{prince_BART}.
#'
#' @param x,object A princebart object.
#' @param type Character; type of estimand: \code{"mixed"} (default) for mixed
#'   estimands averaging over principal stratum uncertainty, or \code{"sample"}
#'   for sample estimands using imputed principal strata (experimental).
#' @param treated_only Logical; if \code{TRUE}, compute treatment effect only
#'   among the treated units (ATT). Default is \code{FALSE} (ATE).
#' @param ... Additional arguments (currently ignored).
#'
#' @return
#' \code{print}: Invisibly returns the object.
#' \code{summary}: Prints and invisibly returns a list with strata distribution
#'   and treatment effect estimates.
#' \code{coef}: Posterior summary of the requested estimand.
#'
#' @examples
#' \dontrun{
#' fit <- prince_BART(Y ~ X | Z | W, data = df)
#'
#' # Mixed ATE for compliers (default)
#' summary(fit)
#' coef(fit)
#'
#' # Mixed ATT for compliers
#' summary(fit, treated_only = TRUE)
#' coef(fit, treated_only = TRUE)
#'
#' # Sample estimands (experimental)
#' coef(fit, type = "sample")
#' }
NULL


#' @rdname princebart-class
#' @export
print.princebart <- function(x, ...) {
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

  cat("\nUse summary() or coef() to extract treatment effects.\n")
  invisible(x)
}


#' @rdname princebart-class
#' @export
summary.princebart <- function(object,
                                type = c("mixed", "sample"),
                                treated_only = FALSE,
                                ...) {
  type <- match.arg(type)

  cat("Principal Stratification BART Summary\n")
  cat("=====================================\n\n")

  # Compute treatment effect
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
  print(strata_summary)
  cat("\n")

  cat(sprintf("%s %s for Compliers:\n", type_label, effect_label))
  print(effect_summary)

  invisible(list(
    strata = strata_summary,
    effect = effect_summary,
    type = type,
    treated_only = treated_only
  ))
}


#' @rdname princebart-class
#' @export
coef.princebart <- function(object,
                            type = c("mixed", "sample"),
                            treated_only = FALSE,
                            ...) {
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
