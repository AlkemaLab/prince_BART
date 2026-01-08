library(parallel)
library(CVXR)

#standard bayesian bootstrap
#Rubin (1981)'s algorithm
bb <- \(n) {
  U <- runif(n - 1)
  G <- diff(c(0, sort(U), 1))
  G
}

# toy example
# c(1, 1, 2) * bb(3)

inc <- \(x, ci) {
  ci <- as.vector(ci)
  lw <- ci[1]
  up <- ci[2]
  x > lw & x < up
}

sm <- \(r) {
  out <- c(apply(r, 1, mean), sqrt(mean(r[2, ] ^ 2)))
  names(out) <- c("cvg", "bias", "sd", "rmse")
  round(out, 3)
}

gen_b <- \(y
  , subpop = dhs$sp
  , psu = dhs$psu
  , wts = dhs$wt
  , seed = 0203
) {
  y <- as.matrix(y)
  y[!subpop, ] <- NA
  #all psus, some maybe empty
  n_psu <- length(unique(psu))
  #size of psu, regardless of missingness
  # note that the weigths do not vary within
  p_c  <- tapply(wts, psu, mean)
  n_c  <- tapply(!is.na(psu), psu, sum)
  D    <- ncol(y)
  set.seed(seed)
  bb_list  <- lapply(seq_len(D), \(r) bb(n_psu))
  yp <- mclapply(seq_len(D), \(k) {
    # this average by psu, generates NA if psu is empty
    psu_y <- tapply(y[, k], psu, mean, na.rm = TRUE)
    # weighted mena considers the relative importnace of each psu
    # it does not matter how much the weights addd up to
    weighted.mean(psu_y, bb_list[[k]] * p_c  * n_c, na.rm = TRUE)
  }, mc.cores = 20) %>% unlist
  est <- mean(yp)
  ci <- as.numeric(quantile(yp, c(.025, .975)))
  psd <- sd(yp)
  c(m = est, s = psd, l = ci[1], u = ci[2])
}



gen_simple_sim <- \(y
  , subpop = dhs$sp
  , psu = dhs$psu
  , wts = dhs$wt
  , seed = 0203
  , prev = .1
  , kappa = .66
) {
  y <- as.matrix(y)
  y[!subpop, ] <- NA
  #all psus, some maybe empty
  n_psu <- length(unique(psu))
  #size of psu, regardless of missingness
  # note that the weigths do not vary within
  p_c  <- tapply(wts, psu, mean)
  n_c  <- tapply(!is.na(psu), psu, sum)
  D    <- ncol(y)
  set.seed(seed)
  bb_list  <- lapply(seq_len(D), \(r) bb(n_psu))
  u <- replicate(ncol(y), rbinom(nrow(y), 1, prev))
  yp <- mclapply(seq_len(D), \(k) {
    psu_y <- tapply(y[, k] - u[, k] * kappa, psu, mean, na.rm = TRUE)
    weighted.mean(psu_y, bb_list[[k]] * p_c  * n_c, na.rm = TRUE)
  }, mc.cores = 20) %>% unlist
  est <- mean(yp)
  ci <- as.numeric(quantile(yp, c(.025, .975)))
  psd <- sd(yp)
  c(m = est, s = psd, l = ci[1], u = ci[2])
}

find_shift_wts <- function(inf, γ, y) {
  wts <- rep(1, length(y))
  new_wts <- rep(0, length(y))
  # uses only not missing data
  # but will return vector of same length
  # completing with zeros
  s <- !is.na(y)
  wts <- as.numeric(wts) [s]
  sum_wts <- sum(wts)
  wts <- wts / sum(wts)
  y <- y[s]
  n <- sum(s)
  r               <- Variable(n)  #shift
  if (inf) {
    objtive       <- Minimize(sum(y * wts * r))
  } else {
    objtive       <- Maximize(sum(y * wts * r))
  }
  constraints   <- list(
    sum(wts * r) == 1
    , r <= γ
    , r >= 1 / γ
  )
  problem       <- Problem(objtive, constraints = constraints)
  result        <- solve(problem, solver = "CLARABEL")
  r             <- round(result$getValue(r), 4)
  new_wts[s] <- wts * r * sum_wts
  print(result$status)
  return(new_wts)
}

###toy example
#ytest <- c(1, 2, 3, NA)
#wmin <- find_shift_wts(inf = TRUE, γ = 1.1, y = ytest)
#wmax <- find_shift_wts(inf = FALSE, γ = 1.1, y = ytest)
#weighted.mean(ytest, w, na.rm = TRUE)
#weighted.mean(ytest, wmin, na.rm = TRUE)
#weighted.mean(ytest, wmax, na.rm = TRUE)
#wmin
#wmax

gen_shift_sim <- \(y
  , subpop = dhs$sp
  , psu = dhs$psu
  , wts = dhs$wt
  , seed = 0203
  , γ = 1.1
  , inf = TRUE
) {
  y <- as.matrix(y)
  y[!subpop, ] <- NA
  #all psus, some maybe empty
  n_psu <- length(unique(psu))
  #size of psu, regardless of missingness
  # note that the weigths do not vary within
  p_c  <- tapply(wts, psu, mean)
  n_c  <- tapply(!is.na(psu), psu, sum)
  D    <- ncol(y)
  set.seed(seed)
  bb_list  <- lapply(seq_len(D), \(r) bb(n_psu))
  yp <- mclapply(seq_len(D), \(k) {
    shift_wts <- find_shift_wts(inf = inf, γ = γ, y = y[, k])
    psu_y <- tapply(y[, k] * shift_wts, psu, mean, na.rm = TRUE)
    weighted.mean(psu_y, bb_list[[k]] * p_c  * n_c, na.rm = TRUE)
  }, mc.cores = 20) %>% unlist
  est <- mean(yp)
  ci <- as.numeric(quantile(yp, c(.025, .975)))
  psd <- sd(yp)
  c(m = est, s = psd, l = ci[1], u = ci[2])
}