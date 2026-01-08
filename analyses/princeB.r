library(tidyverse)
library(dbarts)
library(posterior)

#==============================================================================
#'data must contain
#' Z (assigment),
#' W (treatment uptake),
#' Y(outcome), and
#' some covariates
#==============================================================================
#data <- readRDS("/Users/Lucas/Documents/GitHub/Untitled/fp_impact/debug_NA_co_2025-11-17 08:28:36.102286.rds")
#with(data, table(Z, W))

#data <- tdt
## fit_psbart is now internal (.fit_psbart)
  data = NULL
  , n_warmup = 20
  , n_samples = 20
  , save_trees = FALSE
  , k = 2
  , ntrees = 200
  , verbose = FALSE
  , n_initial = 20
  , ...
) {

  # scale covariates
  X <- subset(data, select = -c(Y, Z, W)) |>
    scale() |>
    as.matrix()
  Y <- data$Y
  Z <- data$Z
  W <- data$W


  e <- bart2(X, Z, verbose = FALSE) |> fitted() |> qnorm()
  X <- cbind(X, e)

  n_total   <- n_warmup + n_samples

  # imputed class
  s_nt <- s_at <- matrix(NA_real_, n_samples, nrow(data))

  # trees for new predcitions
  trees <- list()

  #class probabilities & outcome expected values
  p_at <- p_nt <- m_y0co <- m_y1co <-
    m_y0nt <- m_y1at <- matrix(
      NA_real_, n_samples, nrow(data)
    )

  clip <- function(x, lower = 0.001) {
    upper <- 1 - lower
    x[x < lower] <- lower
    x[x > upper] <- upper
    x
  }
  if (n_initial > 0) {
    # Method-of-moment based estiamtes for reference
    intat <- (coef(lm(W ~ ., data.frame(W, X)[Z == 0, ]))[1]) |> clip()
    intnt <- (coef(lm(W == 0 ~ ., data.frame(W, X)[Z == 1, ]))[1]) |> clip()
    intco <- (1 - intat - intnt) |> clip()

    get_int <- function(a, b) {
      s <- W == a & Z == b
      if (sum(s) > 4) {
        coef(lm(Y ~ ., data.frame(Y, X)[s, ]))[1] |> clip()
      } else {
        coef(lm(Y ~ ., data.frame(Y, X)[W == a, ]))[1] |> clip()
      }
    }
    inty1at <- get_int(1, 0) #at
    inty0nt <- get_int(0, 1) #nt
    inty1   <- get_int(1, 1) #co or at
    inty0   <- get_int(0, 0) #co or nt

    intco_noat <- (intco / (intco + intnt)) |> clip()
    intco_nont <- (intco / (intco + intat)) |> clip()
    inty1co <- ((inty1 - inty1at * (1 - intco_nont)) / intco_nont) |> clip()
    inty0co <- ((inty0 - inty0nt * (1 - intco_noat)) / intco_noat) |> clip()
    c(intat, intnt, intco, inty1at, inty0nt, inty1co, inty0co)
  } else {
    intat <- intnt <- intco <- inty1at <- inty0nt <- inty1co <- inty0co <- 0.5
  }



  # Starting values for G nt, at, co
  # everyone with Z == 1 W == 1 and Z == 0 W == 0 is initially complier
  # to endure there is a sizable number of compliers
  nt <- (Z == 1) * (W == 0)
  at <- (Z == 0) * (W == 1)
  co <- 1 - nt - at

  # Based on simulations a bit of thinning is needed.
  control <- dbartsControl(
    updateState = TRUE
    , keepTrees = save_trees
    , verbose = FALSE
    , n.burn = 0L
    , n.samples = 1L
    , n.thin = 20L
    , n.chains = 1L
    , n.trees = ntrees
    # together with mc.reset.stream() this should ensure reproducibility
    , n.threads = 1
  )

  sampler_co <- dbarts(
    X
    , co
    , test = X
    , offset = qnorm(intco)
    , control = control
    , node.prior = normal(k)
  )
  sampler_atnoco <- dbarts(
    X
    , at
    , test = X
    , subset =  co == 0
    , offset = qnorm(intat / (intat + intnt))
    , control = control
    , node.prior = normal(k)
  )

  if (sum(nt) > 0) nt1 <- nt else nt1 <- W == 0
  sampler_y0nt <- dbarts(
    X
    , Y
    , test = X
    , subset =  (nt1 == 1)
    , offset = qnorm(inty0nt)
    , control = control
    , node.prior = normal(k)
  )

  if (sum(at) > 0) at1 <- at else at1 <- W == 1
  sampler_y1at <- dbarts(
    X
    , Y
    , test = X
    , subset =  (at1 == 1)
    , offset = qnorm(inty1at)
    , control = control
    , node.prior = normal(k)
  )
  sampler_y1co <- dbarts(
    X
    , Y
    , test = X
    , subset =  (co == 1) & (Z == 1)
    , offset = qnorm(inty1co)
    , control = control
    , node.prior = normal(k)
  )
  sampler_y0co <- dbarts(
    X
    , Y
    , test = X
    , subset =  (co == 1) & (Z == 0)
    , offset = qnorm(inty0co)
    , control = control
    , node.prior = normal(k)
  )

  ### run from prior
  sampler_co$sampleTreesFromPrior()
  sampler_atnoco$sampleTreesFromPrior()
  sampler_y0nt$sampleTreesFromPrior()
  sampler_y1at$sampleTreesFromPrior()
  sampler_y0co$sampleTreesFromPrior()
  sampler_y1co$sampleTreesFromPrior()


  for (i in seq_len(n_total)) {

    # Draw a sample from the posterior
    sample_co     <- sampler_co$run()
    sample_atnoco <- sampler_atnoco$run()
    sample_y0nt   <- sampler_y0nt$run()
    sample_y1at   <- sampler_y1at$run()
    sample_y0co   <- sampler_y0co$run()
    sample_y1co   <- sampler_y1co$run()

    # probabilities based on current run

    pco     <- pnorm(sample_co$test)
    patnoco <- pnorm(sample_atnoco$test)

    pat     <- patnoco * (1 - pco)
    pnt     <- (1 - patnoco) * (1 - pco)

    my0nt   <- pnorm(sample_y0nt$test)
    my1at   <- pnorm(sample_y1at$test)
    my0co   <- pnorm(sample_y0co$test)
    my1co   <- pnorm(sample_y1co$test)

    #' update class probabilities using Bayes rule
    #' based on observed Y
    #' and current estimates of p's and m's
    pcoy0  <-
      Y * (pco * my0co /
        (pco * my0co + pnt * my0nt)
      ) +
      (1 - Y) * (pco * (1 - my0co) /
        (pco * (1 - my0co) + pnt * (1 - my0nt))
      )

    pcoy1  <-
      Y * (pco * my1co /
        (pco * my1co + pat * my1at)
      ) +
      (1 - Y) * (pco * (1 - my1co) /
        (pco * (1 - my1co) + pat * (1 - my1at))
      )

    #' Draw a sample from the posterior class probabilities
    #' NOTE: these are imputed discrete values

    nty <- rbinom(nrow(data), 1, (1 - pcoy0))
    aty <- rbinom(nrow(data), 1, (1 - pcoy1))
    # impute new valeus when missing
    nt <- (Z == 1) * (W == 0) + (Z == 0) * (W == 0) * nty
    at <- (Z == 0) * (W == 1) + (Z == 1) * (W == 1) * aty
    co <- 1 - nt - at
    #c(sum(co), sum(nt), sum(at)) |> print()


    #' use updated values to create new data for sampler
    #' (becouse the low proportion of compliers
    #' sometimes there may be none
    #' we skip the update in those cases
    #if (is.na(sum(co))) {
    #  saveRDS(data, file = paste0("debug_NA_co_", Sys.time(), ".rds"))
    #  stop("NA in co")
    #}
    if (sum(co) > 0) {
      if (i < n_initial) {
        dt_co     <- dbartsData(X, co, test = X)
        dt_atnoco <- dbartsData(X, at, test = X, subset = co == 0)
      } else {
        dt_co     <- dbartsData(X, co, test = X, offset = 0)
        dt_atnoco <- dbartsData(X, at, test = X, subset = co == 0, offset = 0)
      }
      sampler_co$setData(dt_co)
      sampler_atnoco$setData(dt_atnoco)
    } else {
      print(paste("skip co update ", i))
    }

    if (sum(co * Z) > 0) {
      if (i < n_initial) {
        dt_y1co   <- dbartsData(X, Y, test = X, subset = co == 1 & (Z == 1))
      } else {
        dt_y1co   <- dbartsData(X, Y, test = X, subset = co == 1 & (Z == 1), offset = 0)
      }
      sampler_y1co$setData(dt_y1co)
    } else {
      print(paste("skip co update ", i))
    }

    if (sum(co * (Z == 0)) > 0) {
      if (i < n_initial) {
        dt_y0co   <- dbartsData(X, Y, test = X, subset = co == 1 & (Z == 0))
      } else {
        dt_y0co   <- dbartsData(X, Y, test = X, subset = co == 1 & (Z == 0), offset = 0)
      }
      sampler_y0co$setData(dt_y0co)
    } else {
      print(paste("skip co update ", i))
    }

    if (sum(nt) > 0) {
      if (i < n_initial) {
        dt_y0nt   <- dbartsData(X, Y, test = X, subset = nt == 1)
      } else {
        dt_y0nt   <- dbartsData(X, Y, test = X, subset = nt == 1, offset = 0)
      }
      sampler_y0nt$setData(dt_y0nt)
    } else {
      print(paste("skip nt update ", i))
    }

    if (sum(at) > 0) {
      if (i < 40) {
        dt_y1at   <- dbartsData(X, Y, test = X, subset = at == 1)
      } else {
        dt_y1at   <- dbartsData(X, Y, test = X, subset = at == 1, offset = 0)
      }
      sampler_y1at$setData(dt_y1at)
    } else {
      print(paste("skip at update ", i))
    }


    #' Instead of saving imputations
    #' we can we can save probabilities

    # Store samples if no longer warming up.
    if (i > n_warmup) {
      offset <- i - n_warmup
      s_nt[offset, ] <- nt
      s_at[offset, ] <- at

      p_at[offset, ]   <- pat
      p_nt[offset, ]   <- pnt

      m_y0co[offset, ] <- my0co
      m_y1co[offset, ] <- my1co
      m_y0nt[offset, ] <- my0nt
      m_y1at[offset, ] <- my1at

      if (save_trees) {
        list_trees <- list(
          cbind(m = "co", sampler_co$getTrees())
          , cbind(m = "atnoco", sampler_atnoco$getTrees())
          , cbind(m = "y1co", sampler_y1co$getTrees())
          , cbind(m = "y0co", sampler_y0co$getTrees())
          , cbind(m = "y1at", sampler_y1at$getTrees())
          , cbind(m = "y0nt", sampler_y0nt$getTrees())
        )
        trees[[offset]] <- do.call(rbind, list_trees)
        trees[[offset]]$iteration <- offset
      }
    }

    if (verbose == TRUE && i %% 20 == 0) print(i)
  }
  if (save_trees) {
    trees  <- do.call(rbind, trees)
    trees$sample  <- NULL
  } else {
    trees <- NULL
  }

  imputed <- array(c(s_nt, s_at)
    , dim = c(n_samples, nrow(data), 2)
  )
  probs   <- array(c(p_at, p_nt, m_y0co, m_y1co, m_y0nt, m_y1at)
    , dim = c(n_samples, nrow(data), 6)
  )


  results <- list(
    imputed = imputed
    , trees = trees
    , probs = probs
  )
  print("done")
  results
}



## mc_psbart is now prince_BART
  dt = NULL
  , n_warmup = 5L
  , n_samples = 5L
  , n_chains = 4
  , keep_trees = FALSE
  , k = 2
  , n_trees = 200
  , initial = 0
) {
  res0 <- mclapply(1:n_chains
    , .fit_psbart
    , data = dt
    , n_warmup = n_warmup
    , n_samples = n_samples
    , mc.cores = min(parallel::detectCores(), n_chains)
    , save_trees = keep_trees
    , k = k
    , ntrees = n_trees
    , n_initial = initial
  )
  res <- list()
  if (keep_trees) {
    list_tree <- lapply(res0, \(x) x$trees)
    res$trees <- bind_rows(list_tree, .id = "chain")
    res$trees$chain <- as.integer(res$trees$chain)
  } else {
    res$trees <- NULL
  }

  list_imp <- lapply(res0, \(x) x$imputed)
  res$imp <- abind::abind(list_imp, along = 4)
  res$imp <- aperm(res$imp, c(1, 4, 3, 2))
  names(attributes(res$imp)$dimnames) <-
    c("iteration", "chain", "variable", "unit")
  dimnames(res$imp)$variable <- c("nt", "at")

  list_probs <- lapply(res0, \(x) x$probs)
  res$probs <- abind::abind(list_probs, along = 4)
  res$probs <- aperm(res$probs, c(1, 4, 3, 2))
  names(attributes(res$probs)$dimnames) <-
    c("iteration", "chain", "variable", "unit")
  dimnames(res$probs)$variable <- c(
    "p_n", "p_a", "m_y0c", "m_y1c", "m_y0n", "m_y1a"
  )
  #save data
  res$data <- dt

  res
}



# impute potential outcomes amomng compliers
# for sample estimands
imput_potentialoutcomes_c <- \(prince_bart_fit) {
  data  <- prince_bart_fit$data
  probs <- prince_bart_fit$probs

  Y <- data$Y
  Z <- data$Z
  N <- nrow(data)


  my0co <- probs[, , "m_y0c", ]
  my1co <- probs[, , "m_y1c", ]

  ### 'draw sample for outcome posterior predcitvie distribution
  y0c <- apply(my0co, 1:2, \(p) (Z == 1) * rbinom(N, 1, p) + (Z == 0) * Y)
  y1c <- apply(my1co, 1:2, \(p) (Z == 1) * Y + (Z == 0) * rbinom(N, 1, p))


  #now with induced correlation
  ycc <- apply(probs[, , c("m_y0c", "m_y1c"), ], 1:2, \(m) {
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
    s0 <- m0 <- k0 <- s1 <- m1 <- k1 <- NULL

    y0cc <- (Z == 1) * rbinom(N, 1, mu0cc) + (Z == 0) * Y
    y1cc <- (Z == 1) * Y + (Z == 0) * rbinom(N, 1, mu1cc)

    data.frame(y0cc, y1cc)
  }, simplify = FALSE)

  y0cc  <- apply(ycc, 1:2, \(x) x[[1]]$y0cc)
  y1cc  <- apply(ycc, 1:2, \(x) x[[1]]$y1cc)

  imp <- abind::abind(list(y0c, y1c, y0cc, y1cc), along = 4)
  imp <- aperm(imp, c(2, 3, 4, 1))
  names(attributes(imp)$dimnames) <-
    c("iteration", "chain", "variable", "unit")
  dimnames(imp)$variable <- c("y0", "y1", "cy0", "cy1")
  imp
}

# compute quantities of interest
get_mix_tau <- \(p_arr, treated = NULL) {
  if (is.null(treated)) {
    treated <- rep(TRUE, dim(p_arr)[4])
    outname <- c("y0c", "y1c", "y0n", "y1a", "mate_c")
    sname <- c("c", "n", "a")
  } else {
    outname <- c("y0c1", "y1c1", "y0n1", "y1a1", "matt_c")
    sname <- c("c1", "n1", "a1")
  }

  p_n <- p_arr[, , "p_n", treated]
  p_a <- p_arr[, , "p_a", treated]
  p_c <- 1 - p_n - p_a

  p_g <- abind::abind(p_c, p_n, p_a, along = 4)
  strata_prob   <- apply(p_g, c(1:2, 4), mean)

  m_y <- p_arr[, , 3:6, treated]
  str <- c(1, 1, 2, 3)
  mean_pout <- lapply(1:4, \(g) {
    numer <- apply(m_y[, , g, ] * p_g[, , , str[g]], 1:2, mean)
    numer / strata_prob[, , str[g]]
  })

  mean_effect     <- mean_pout
  mean_effect$tau <- mean_pout[[2]] - mean_pout[[1]] #!!!!
  mean_effect <- abind::abind(mean_effect, along = 3)

  names(attributes(mean_effect)$dimnames) <-
    c("iteration", "chain", "variable")
  dimnames(mean_effect)$variable <- outname
  mean_effect <- as_draws_array(mean_effect)

  names(attributes(strata_prob)$dimnames) <-
    c("iteration", "chain", "variable")
  dimnames(strata_prob)$variable <- sname
  strata_prob <- as_draws_array(strata_prob)

  list(strata_prob, mean_effect)
}


get_sample_tau <- \(imp_g, imp_o
  , treated = NULL
  , include_corr = FALSE
) {
  if (is.null(treated)) {
    treated <- rep(TRUE, nrow(imp_g))
    resvarname <- c("co", "y0c", "y1c", "sate_c")
  } else {
    resvarname <- c("co1", "y0c1", "y1c1", "satt_c")
  }


  #arrray "iteration" by "chain"
  co <-  (1 - imp_g[, , "nt", ] - imp_g[, , "at", ])
  pc <- apply(co[, , treated], 1:2, mean)

  Y0 <- Y1 <- NULL

  if (include_corr) {
    print("including residual correlation between potential outcomes")
    Y0 <- imp_o[, , "cy0", ]
    Y1 <- imp_o[, , "cy1", ]
  } else {
    Y0 <- imp_o[, , "y0", ]
    Y1 <- imp_o[, , "y1", ]
  }

  Y0[co == 0] <- NA
  Y1[co == 0] <- NA
  tau <- Y1 - Y0

  Y0c <- apply(Y0[, , treated], 1:2, mean, na.rm = TRUE)
  Y1c <- apply(Y1[, , treated], 1:2, mean, na.rm = TRUE)
  tau_c <- apply(tau[, , treated], 1:2, mean, na.rm = TRUE)


  kqts <- abind::abind(list(
    pc, Y0c, Y1c, tau_c
  ), along = 3)

  names(attributes(kqts)$dimnames) <-
    c("iteration", "chain", "variable")

  dimnames(kqts)$variable <- resvarname

  post <- as_draws_array(kqts)
  post
}


#================================================================
###* predict function
#================================================================
#predict.trees <- \(trees
#  , scaled_x
#  , n_cores = 12
#) {
#  trees    <- as.data.frame(trees)
#  scaled_x <- as.matrix(scaled_x)
#  n_trees  <- max(trees$tree)

#  predictions   <-  mclapply(unique(trees$sample), \(sidx) {
#    these_trees <- trees[trees$sample == sidx, ]
#    pre_mat     <- sapply(1:n_trees, \(i) {
#      single_tree <- these_trees[these_trees$tree == i, ]
#      getPredictionsForTree(single_tree, scaled_x)
#    })
#    apply(pre_mat, 1, sum)
#  }, mc.cores = n_cores)
#  pred_arr <- simplify2array(predictions)
#  p <- pnorm(pred_arr)
#  p
#}
predict_one_draw <- function(trees_df, x, draw_idx) {
  these <- trees_df[trees_df$sample == draw_idx, ]    # trees for this draw
  n_trees <- max(these$tree)
  pre_mat <- sapply(seq_len(n_trees), \(i) {
    getPredictionsForTree(these[these$tree == i, ], x)
  })
  pnorm(apply(pre_mat, 1, sum))                       # vector of predictions
}

predict.trees <- function(trees, scaled_x, n_cores = 12) {
  trees <- as.data.frame(trees)
  samples <- unique(trees$sample)

  # allow a single matrix or a list of matrices, one per sample
  x_list <- if (is.list(scaled_x)) scaled_x else rep(list(scaled_x), length(samples))
  stopifnot(length(x_list) == length(samples))

  preds <- mapply(samples, x_list, SIMPLIFY = FALSE, USE.NAMES = FALSE,
    FUN = function(sidx, xmat) {
      xmat <- as.matrix(xmat)
      n_trees <- max(trees$tree)
      these <- trees[trees$sample == sidx, ]
      pre_mat <- sapply(seq_len(n_trees), \(i) {
        getPredictionsForTree(these[these$tree == i, ], xmat)
      })
      apply(pre_mat, 1, sum)
    })

  pnorm(simplify2array(preds))
}
# imputed draws: n x length(imputed_cols) x n_draws
# trees_df: has one row per split with a 'sample' column


getPredictionsForTree <- function(tree, x) {
  predictions <- rep(NA_real_, nrow(x))
  getPredictionsForTreeRecursive <- function(tree, indices) { 
    if (tree$var[1] == -1) {
      # Assigns in the calling environment by using <<-
      predictions[indices] <<- tree$value[1]
      return(1)
    }
  goesLeft <- x[indices, tree$var[1]] <= tree$value[1]
  headOfLeftBranch <- tree[-1,]
  n_nodes.left <- getPredictionsForTreeRecursive(
  headOfLeftBranch, indices[goesLeft])
  headOfRightBranch <- tree[seq.int(2 + n_nodes.left, nrow(tree)),]
  n_nodes.right <- getPredictionsForTreeRecursive(
  headOfRightBranch, indices[!goesLeft])
  return(1 + n_nodes.left + n_nodes.right)
  }
  getPredictionsForTreeRecursive(tree, seq_len(nrow(x)))
  return(predictions)
}

#getPredictionsForTree(treeOfInterest, bartFit$fit$data@x[1:5,])

#=============================================================
#summary functions
#=============================================================
mate_c <- \(prince_bart_fit) {
  prob <- prince_bart_fit$probs
  res <- get_mix_tau(prob)
  lapply(res, summary)
}

matt_c <- \(prince_bart_fit) {
  treated <- prince_bart_fit$data$Z == 1
  prob <- prince_bart_fit$probs
  res <- get_mix_tau(prob, treated)
  lapply(res, summary)
}


satt_c <- \(prince_bart_fit, induce_residual_corr = FALSE) {
  treated <- prince_bart_fit$data$Z == 1
  impo <- imput_potentialoutcomes_c(prince_bart_fit)
  impg <- prince_bart_fit$imp
  post_satt_c <- get_sample_tau(impg, impo, treated
    , include_corr = induce_residual_corr
  )
  summary(post_satt_c)
}

sate_c <- \(prince_bart_fit, induce_residual_corr = FALSE) {
  impo <- imput_potentialoutcomes_c(prince_bart_fit)
  impg <- prince_bart_fit$imp
  post_satt_c <- get_sample_tau(impg, impo
    , include_corr = induce_residual_corr
  )
  summary(post_satt_c)
}

validateArgumentsInEnvironment <- dbarts:::validateArgumentsInEnvironment
redirectCall <- dbarts:::redirectCall
addCallArgument <- dbarts:::addCallArgument
quoteInNamespace <- dbarts:::quoteInNamespace
parsePriors <- dbarts:::parsePriors
setDefaultsFromFormals <- dbarts:::setDefaultsFromFormals

dbarts <- function (formula, data, test, subset, weights, offset, offset.test = offset, 
    verbose = FALSE, n.samples = 800L, tree.prior = cgm, node.prior = normal, 
    resid.prior = chisq, proposal.probs = c(birth_death = 0.5, 
        swap = 0.1, change = 0.4, birth = 0.5), control = dbarts::dbartsControl(), 
    sigma = NA_real_) 
{
    matchedCall <- match.call()
    evalEnv <- parent.frame(1L)
    validateCall <- redirectCall(matchedCall, quoteInNamespace(validateArgumentsInEnvironment))
    validateCall <- addCallArgument(validateCall, 1L, sys.frame(sys.nframe()))
    validateCall <- addCallArgument(validateCall, 2L, dbarts::dbarts)
    eval(validateCall, evalEnv, getNamespace("dbarts"))

    if (length(control@call) == 1L && control@call == call("NA")) 
        control@call <- matchedCall
    control@verbose <- verbose
    dataCall <- redirectCall(matchedCall, quoteInNamespace(dbartsData))
    data <- eval(dataCall, evalEnv)
    data@n.cuts <- rep_len(attr(control, "n.cuts"), ncol(data@x))


    data@sigma <- sigma
    attr(control, "n.cuts") <- NULL
    uniqueResponses <- unique(data@y)

    if (length(uniqueResponses) == 2 && all(sort(uniqueResponses) == c(0, 1))) control@binary <- TRUE
    if (length(uniqueResponses) == 1 && all(uniqueResponses == 0 | uniqueResponses == 1)) control@binary <- TRUE
    if (is.na(data@sigma) && !control@binary) 
        data@sigma <- summary(lm(data@y ~ data@x, weights = data@weights, offset = data@offset))$sigma
    if (!control@binary && !is.null(data@offset) && all(data@offset == 0)) {
        data@offset <- NULL
    }
    if (!control@binary && !is.null(data@offset.test) && all(data@offset.test == 
        0)) {
        data@offset.test <- NULL
    }
    parsePriorsCall <- redirectCall(matchedCall, quoteInNamespace(parsePriors))
    parsePriorsCall <- setDefaultsFromFormals(parsePriorsCall, 
        formals(dbarts), "tree.prior", "node.prior", "resid.prior")
    parsePriorsCall$control <- control
    parsePriorsCall$data <- data
    parsePriorsCall$parentEnv <- evalEnv
    if (control@binary) {
        if (any(names(parsePriorsCall) == "resid.prior")) {
            parsePriorsCall[[which(names(parsePriorsCall) == 
                "resid.prior")]] <- quote(fixed(1))
        }
        else {
            parsePriorsCall[[length(parsePriorsCall) + 1L]] <- quote(fixed(1))
            names(parsePriorsCall)[length(parsePriorsCall)] <- "resid.prior"
        }
    }
    priors <- eval(parsePriorsCall)
    model <- new("dbartsModel", priors$tree.prior, priors$node.prior, 
        priors$node.hyperprior, priors$resid.prior, proposal.probs = proposal.probs, 
        node.scale = if (control@binary) 
            3
        else 0.5)
    result <- new("dbartsSampler", control, model, data)
    result
}