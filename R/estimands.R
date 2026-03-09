#' Treatment Effect Estimands for Compliers
#'
#' Internal functions to compute various treatment effect estimands for the complier
#' stratum from a fitted princebart model.
#'
#' @param prince_bart_fit A fitted object from \code{prince_BART}.
#' @param induce_residual_corr Logical; for sample estimands, whether to induce
#'   residual correlation between potential outcomes (default: FALSE).
#'
#' @return A \code{posterior} summary object with posterior mean, median,
#'   standard deviation, and quantiles.
#'
#' @name estimands
#' @keywords internal
NULL


#' @describeIn estimands Mixed Average Treatment Effect for Compliers (MATE_C)
#' @keywords internal
mate_c <- function(prince_bart_fit) {
  prob <- prince_bart_fit$probs
  res <- get_mix_tau(prob)
  lapply(res, summary)
}


#' @describeIn estimands Mixed Average Treatment Effect on the Treated
#'   for Compliers (MATT_C)
#' @keywords internal
matt_c <- function(prince_bart_fit) {
  Z <- prince_bart_fit$data$Z
  treated <- Z == 1
  prob <- prince_bart_fit$probs
  res <- get_mix_tau(prob, treated)
  lapply(res, summary)
}


#' @describeIn estimands Sample Average Treatment Effect for Compliers (SATE_C)
#' @keywords internal
sate_c <- function(prince_bart_fit, induce_residual_corr = FALSE) {
  impo <- imput_potentialoutcomes_c(prince_bart_fit)
  impg <- prince_bart_fit$imp
  res <- get_sample_tau(impg, impo, include_corr = induce_residual_corr)
  lapply(res, summary)
}


#' @describeIn estimands Sample Average Treatment Effect on the Treated
#'   for Compliers (SATT_C)
#' @keywords internal
satt_c <- function(prince_bart_fit, induce_residual_corr = FALSE) {
  Z <- prince_bart_fit$data$Z
  treated <- Z == 1
  impo <- imput_potentialoutcomes_c(prince_bart_fit)
  impg <- prince_bart_fit$imp
  res <- get_sample_tau(impg, impo, treated, include_corr = induce_residual_corr)
  lapply(res, summary)
}


#' Compute Mixed Treatment Effects
#'
#' Internal function to compute mixture-based estimands from posterior
#' probability arrays.
#'
#' @param p_arr 4D probability array from princebart fit.
#' @param treated Optional logical vector indicating treated units.
#'
#' @return List with strata probabilities and mean effects as draws arrays.
#'
#' @keywords internal
get_mix_tau <- function(p_arr, treated = NULL) {
  if (is.null(treated)) {
    treated <- rep(TRUE, dim(p_arr)[4])
    outname <- c("Y(0) | compliers", "Y(1) | compliers", "Y(0) | never-takers", 
                 "Y(1) | always-takers", "MATE | compliers")
    sname <- c("compliers", "never-takers", "always-takers")
  } else {
    outname <- c("Y(0) | compliers, Z=1", "Y(1) | compliers, Z=1", 
                 "Y(0) | never-takers, Z=1", "Y(1) | always-takers, Z=1", 
                 "MATT | compliers")
    sname <- c("compliers, Z=1", "never-takers, Z=1", "always-takers, Z=1")
  }

  # Extract strata probabilities for treated units
  # p_arr dimensions: (iter, chain, var, units)
  # Variable indices: 1=p_a, 2=p_n, 3:6=outcomes
  p_a_arr <- p_arr[, , 1, treated, drop = FALSE]  # (iter, chain, 1, n_treated)
  p_n_arr <- p_arr[, , 2, treated, drop = FALSE]  # (iter, chain, 1, n_treated)
  
  # Compute p_c = 1 - p_n - p_a
  p_c_arr <- 1 - p_n_arr - p_a_arr
  
  # Reshape: aperm to move variable dim to end, then select that dim
  # From (iter, chain, 1, units) -> (iter, chain, units, 1)
  p_a <- aperm(p_a_arr, c(1, 2, 4, 3))[, , , 1, drop = FALSE]
  p_n <- aperm(p_n_arr, c(1, 2, 4, 3))[, , , 1, drop = FALSE]
  p_c <- aperm(p_c_arr, c(1, 2, 4, 3))[, , , 1, drop = FALSE]
  
  # Now p_a, p_n, p_c are (iter, chain, units, 1). Use abind to stack along dim 4
  p_g <- abind::abind(p_a, p_n, p_c, along = 4)  # Result: (iter, chain, units, 3)
  
  strata_prob <- apply(p_g, c(1:2, 4), mean)  # Keep iter, chain, strata; average units

  m_y_arr <- p_arr[, , 3:6, treated, drop = FALSE]  # (iter, chain, 4_vars, n_treated)
  # Reshape: aperm to move var dim to end
  # From (iter, chain, 4, units) -> (iter, chain, units, 4)
  m_y <- aperm(m_y_arr, c(1, 2, 4, 3))
  
  str <- c(1, 1, 2, 3)

  mean_pout <- lapply(1:4, function(g) {
    # m_y[, , , g] is (iter, chain, units)
    # p_g[, , , str[g]] is (iter, chain, units)
    numer <- apply(m_y[, , , g] * p_g[, , , str[g]], 1:2, mean)
    numer / strata_prob[, , str[g]]
  })

  mean_effect <- mean_pout
  mean_effect$tau <- mean_pout[[2]] - mean_pout[[1]]
  mean_effect <- abind::abind(mean_effect, along = 3)

  dimnames(mean_effect) <- list(
    iteration = NULL,
    chain = NULL,
    variable = outname
  )
  mean_effect <- posterior::as_draws_array(mean_effect)

  dimnames(strata_prob) <- list(
    iteration = NULL,
    chain = NULL,
    variable = sname
  )
  strata_prob <- posterior::as_draws_array(strata_prob)

  list(strata_prob, mean_effect)
}


#' Compute Sample-Based Treatment Effects
#'
#' Internal function to compute sample-based estimands from imputed
#' compliance classes and outcomes.
#'
#' @param imp_g Imputed compliance class array.
#' @param imp_o Imputed potential outcomes array.
#' @param treated Optional logical vector indicating treated units.
#' @param include_corr Logical; use correlated imputations.
#'
#' @return A draws array with posterior samples.
#'
#' @keywords internal
get_sample_tau <- function(imp_g, imp_o, treated = NULL, include_corr = TRUE) {
  if (is.null(treated)) {
    treated <- rep(TRUE, dim(imp_g)[4])
    outname <- c("Y(0) | compliers", "Y(1) | compliers", "Y(0) | never-takers", 
                 "Y(1) | always-takers", "SATE | compliers")
    sname <- c("compliers", "never-takers", "always-takers")
  } else {
    outname <- c("Y(0) | compliers, Z=1", "Y(1) | compliers, Z=1", 
                 "Y(0) | never-takers, Z=1", "Y(1) | always-takers, Z=1", 
                 "SATT | compliers")
    sname <- c("compliers, Z=1", "never-takers, Z=1", "always-takers, Z=1")
  }

  # Complier indicator
  nt <- imp_g[, , "nt", ]
  at <- imp_g[, , "at", ]
  co <- 1 - nt - at
  p_g <- abind::abind(co, nt, at, along = 4)
  strata_prob <- apply(p_g, c(1:2, 4), mean)

  # Select outcome imputations
  if (include_corr) {
    Y0 <- imp_o[, , "cy0", ]
    Y1 <- imp_o[, , "cy1", ]
  } else {
    Y0 <- imp_o[, , "y0", ]
    Y1 <- imp_o[, , "y1", ]
  }

  Y0n <- Y0
  Y1a <- Y1
  Y0n[co == 1] <- NA
  Y1a[co == 1] <- NA

  # Mask non-compliers
  Y0[co == 0] <- NA
  Y1[co == 0] <- NA

  mY0c <- apply(Y0[, , treated], 1:2, mean, na.rm = TRUE)
  mY1c <- apply(Y1[, , treated], 1:2, mean, na.rm = TRUE)
  mY0n <- apply(Y0n[, , treated], 1:2, mean, na.rm = TRUE)
  mY1a <- apply(Y1a[, , treated], 1:2, mean, na.rm = TRUE)

  tau <- Y1 - Y0
  tau_c <- apply(tau[, , treated], 1:2, mean, na.rm = TRUE)

  mean_effect <- abind::abind(list(mY0c, mY1c, mY0n, mY1a, tau_c), along = 3)

  dimnames(mean_effect) <- list(
    iteration = NULL,
    chain = NULL,
    variable = outname
  )
  mean_effect <- posterior::as_draws_array(mean_effect)

  dimnames(strata_prob) <- list(
    iteration = NULL,
    chain = NULL,
    variable = sname
  )
  strata_prob <- posterior::as_draws_array(strata_prob)

  list(strata_prob, mean_effect)
}


# =============================================================================
# Ordinal Uptake Estimands
# =============================================================================

#' Ordinal Complier-Like Contrast Estimands
#'
#' Compute posterior summaries for contrasts in ordinal uptake settings.
#' Focuses on the monotone compliance group (W(0) - W(1) = 1) and
#' stratifies by baseline uptake W(0) up to level 5.
#'
#' @param prince_bart_fit A fitted ordinal prince_bart object.
#'
#' @return A list with the following elements:
#'   \item{overall}{Data frame of posterior summary for the overall W(0)-W(1)=1 contrast.}
#'   \item{by_w0}{Data frame of posterior summaries stratified by W(0)=j for j=0..5.}
#'
#' @keywords internal
estimands_ordinal_mixed <- function(prince_bart_fit) {
  imp <- prince_bart_fit$imp
  probs <- prince_bart_fit$probs

  w0 <- imp[, , "w0", , drop = FALSE][, , 1, , drop = FALSE]
  w1 <- imp[, , "w1", , drop = FALSE][, , 1, , drop = FALSE]
  y0 <- probs[, , "m_y0", , drop = FALSE][, , 1, , drop = FALSE]
  y1 <- probs[, , "m_y1", , drop = FALSE][, , 1, , drop = FALSE]

  # Overall contrast: complier-like units where W(0) - W(1) = 1
  co <- (w0 - w1) == 1

  # Treatment effect-like contrast: mean Y(1) - mean Y(0) | complier-like
  d_p <- y0 - y1
  d_p[!co] <- NA

  # Average over units per iteration-chain pair (keep iter x chain shape)
  md_overall <- apply(d_p, c(1, 2), function(x) {
    if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)
  })
  md_overall <- array(md_overall, dim = c(dim(md_overall)[1], dim(md_overall)[2], 1))
  md_overall <- posterior::as_draws_array(md_overall)
  dimnames(md_overall)[[3]] <- "W(0) - W(1) = 1 contrast"

  overall_summary <- summary(md_overall)

  # Stratified by W(0) = j for j = 0, 1, ..., 5
  w0_vals <- w0
  max_w0 <- min(5, max(w0_vals, na.rm = TRUE))
  
  md_by_w0_list <- lapply(0:max_w0, function(j) {
    co_j <- ((w0 - w1) == 1) & (w0 == j)
    d_p_j <- y0 - y1
    d_p_j[!co_j] <- NA
    md_j <- apply(d_p_j, c(1, 2), function(x) {
      if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)
    })
    md_j
  })

  # Combine into draws array with variable names
  md_by_w0 <- abind::abind(md_by_w0_list, along = 3)
  dimnames(md_by_w0) <- list(
    iteration = NULL,
    chain = NULL,
    variable = paste0("W(0)=", 0:max_w0, " | W(0)-W(1)=1")
  )
  md_by_w0 <- posterior::as_draws_array(md_by_w0)
  by_w0_summary <- summary(md_by_w0)

  list(
    overall = overall_summary,
    by_w0 = by_w0_summary
  )
}
