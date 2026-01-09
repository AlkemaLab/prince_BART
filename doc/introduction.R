## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  eval = TRUE,
  fig.width = 8,
  fig.height = 5,
  out.width = "100%"
)

## ----simulate-data------------------------------------------------------------
library(princeBART)
set.seed(42)
n <- 800

# Covariates
X <- data.frame(
  age = rnorm(n, 40, 10),
  education = rnorm(n, 12, 3),
  income = rnorm(n, 50000, 15000)
)

# Principal strata (latent)
# Probability of being a complier depends on covariates
p_complier <- plogis(-1.2 + 0.03 * X$age + 0.12 * X$education - 0.00001 * X$income)
p_always   <- plogis(-2 + 0.01 * X$age)
p_never    <- pmax(1 - p_complier - p_always, 0)
denom <- p_complier + p_never + p_always
p_complier <- p_complier / denom
p_never <- p_never / denom
p_always <- p_always / denom


strata <- sapply(1:n, function(i) {
  sample(c("complier", "never", "always"), 1,
         prob = c(p_complier[i], p_never[i], p_always[i]))
})

# Instrument (randomized, possibly conditional on X)
Z <- rbinom(n, 1, 0.5)

# Treatment uptake depends on strata and instrument
W <- Z
W[strata == "always"] <- 1L
W[strata == "never"]  <- 0L

# Heterogeneous complier treatment effect
tau <- 0.2 +
  0.15 * (X$education > 12) -
  0.10 * (X$age > 50) +
  0.08 * (X$income > 60000)

# Potential outcomes (only Y(W) is observed)
lin0 <- -0.6 + 0.01 * X$age + 0.04 * X$education + 0.00001 * X$income
Y0 <- plogis(lin0)
Y1 <- plogis(lin0 + tau)

Y <- rbinom(n, 1, ifelse(W == 1, Y1, Y0))

# Combine into data frame
simdata <- data.frame(X, Z = Z, W = W, Y = Y)

# True ATE among compliers
true_ate_c <- mean(tau[strata == "complier"])

## ----formula-usage, eval = FALSE----------------------------------------------
# 
# library(future)
# 
# # Set up parallel backend (works on all platforms)
# plan(multisession, workers = 4)
# 
# # Fit the model
# fit <- prince_BART(
#   Y ~ age + education + income | Z | W,
#   data = simdata,
#   n_chains = 2,
#   n_warmup = 100,
#   n_samples = 200,
#   instrument_overlap = c(0.1, 0.9),
#   keep_trees = TRUE  # Save trees for generalization
# )
# 

## ----estimands----------------------------------------------------------------
fit <- readRDS(system.file("extdata", "fit_intro.rds", package = "princeBART"))

# Print summary
print(fit)

# Summary shows strata distribution and treatment effect
summary(fit)

# Mixed ATE for Compliers (default)
coef(fit)

## ----heterogeneity, fig.width=8, fig.height=5, out.width="100%"---------------
res <- segment_heterogeneity(
  fit,
  min_compliers_bucket = 50,
  plot = TRUE,
  contrast = TRUE
)

res$contrast$summary
res$plot$diff


## ----heterogeneity_effects, fig.width=8, fig.height=5, out.width="100%"-------
res$effects
res$plot$effect

