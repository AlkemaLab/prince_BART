#' Impute Potential Outcomes for Compliers
#'
#' Imputes potential outcomes Y(0) and Y(1) for compliers based on
#' posterior samples from a fitted princebart model. Used internally
#' by sample-based estimand functions.
#'
#' @param prince_bart_fit A fitted object from \code{prince_BART}.
#'
#' @return A 4D array with dimensions (iteration, chain, variable, unit)
#'   containing imputed potential outcomes:
#'   \item{y0}{Imputed Y(0) for compliers}
#'   \item{y1}{Imputed Y(1) for compliers}
#'   \item{cy0}{Imputed Y(0) with induced correlation}
#'   \item{cy1}{Imputed Y(1) with induced correlation}
#'
#' @keywords internal
imput_potentialoutcomes_c <- function(prince_bart_fit) {
  data  <- prince_bart_fit$data
  probs <- prince_bart_fit$probs

  Y <- data$Y
  Z <- data$Z
  N <- length(Y)

  my0co <- probs[, , "m_y0c", ]
  my1co <- probs[, , "m_y1c", ]

  # Standard imputations
  y0c <- apply(my0co, 1:2, function(p) {
    (Z == 1) * stats::rbinom(N, 1, p) + (Z == 0) * Y
  })
  y1c <- apply(my1co, 1:2, function(p) {
    (Z == 1) * Y + (Z == 0) * stats::rbinom(N, 1, p)
  })

  # Imputations with induced correlation
  ycc <- apply(probs[, , c("m_y0c", "m_y1c"), ], 1:2, function(m) {
    my0co <- m["m_y0c", ]
    my1co <- m["m_y1c", ]

    s0 <- cbind(my0co, 1 - my0co) < cbind(my1co, 1 - my1co)
    m0 <- cbind(my0co, 1 - my0co) / (cbind(my1co, 1 - my1co) + .0001)
    k0 <- rowSums(m0 * s0)

    s1 <- cbind(my1co, 1 - my1co) < cbind(my0co, 1 - my0co)
    m1 <- cbind(my1co, 1 - my1co) / (cbind(my0co, 1 - my0co) + .0001)
    k1 <- rowSums(m1 * s1)

    mu0cc <- my0co + (Y - my1co) * k0
    mu1cc <- my1co + (Y - my0co) * k1

    y0cc <- (Z == 1) * stats::rbinom(N, 1, mu0cc) + (Z == 0) * Y
    y1cc <- (Z == 1) * Y + (Z == 0) * stats::rbinom(N, 1, mu1cc)

    data.frame(y0cc, y1cc)
  }, simplify = FALSE)

  y0cc <- apply(ycc, 1:2, function(x) x[[1]]$y0cc)
  y1cc <- apply(ycc, 1:2, function(x) x[[1]]$y1cc)

  imp <- abind::abind(list(y0c, y1c, y0cc, y1cc), along = 4)
  imp <- aperm(imp, c(2, 3, 4, 1))

  dimnames(imp) <- list(
    iteration = NULL,
    chain = NULL,
    variable = c("y0", "y1", "cy0", "cy1"),
    unit = NULL
  )

  imp
}
