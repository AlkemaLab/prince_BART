## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  eval = TRUE,
  fig.width = 8,
  fig.height = 5,
  out.width = "100%"
)

## ----simulate-external--------------------------------------------------------
library(princeBART)
set.seed(123)
n_external <- 500

# External population overlaps with source but shifts slightly
external_data <- data.frame(
  age = rnorm(n_external, 43, 11),        # Slightly older on average
  education = rnorm(n_external, 11.5, 3.5) # Slightly lower education
  # Note: income is NOT available in external data
)
# Survey design: 100 PSUs with varying sizes
external_data$psu <- sample(1:100, n_external, replace = TRUE)
external_data$weight <- runif(n_external, 0.5, 2)  # Survey weights

# Define a target subpopulation (e.g., adults under 60)
external_data$eligible <- external_data$age < 60

## ----general-bart, eval = FALSE-----------------------------------------------
# fit <- readRDS(system.file("extdata", "fit_intro.rds", package = "princeBART"))
# 
# # Generalize treatment effects to external population
# # Missing variables are auto-detected and imputed
# pate_result <- general_BART(
#   princebart_fit = fit,
#   newdata = external_data,
#   subpop = external_data$eligible,       # Target subpopulation
#   psu = external_data$psu,               # PSU for survey inference
#   weights = external_data$weight,        # Survey weights
#   n_cores = 4,
#   verbose = TRUE
# )
# 
# # View results
# print(pate_result)

## ----pate-summary-------------------------------------------------------------
pate_result <- readRDS(system.file("extdata", "pate_result.rds", package = "princeBART"))


# Detailed summary
summary(pate_result)

# Access posterior draws for custom analyses
hist(pate_result$draws,
     main = "Posterior Distribution of PATE",
     xlab = "Treatment Effect")
abline(v = pate_result$pate, col = "red", lwd = 2)
abline(v = pate_result$ci, col = "red", lty = 2)

## ----overlap-analysis, eval = FALSE-------------------------------------------
# # Compute overlap scores from the fitted general_pate object
# overlap <- general_BART_overlap(
#   object = pate_result,
#   verbose = TRUE
# )
# 

## ----overlap-analysis-plot----------------------------------------------------
overlap <- readRDS(system.file("extdata", "overlap_result.rds", package = "princeBART"))

# View overlap metrics
head(overlap$overlap)

# Plot overlap diagnostics
plot_overlap(overlap)

## ----overlap-trimming, eval = FALSE-------------------------------------------
# # Trim observations with |e_s_tilde| > 2 (i.e., outliers)
# overlap_trimmed <- general_BART_overlap(
#   object = pate_result,
#   threshold = 2,           # Standardized score threshold
#   overlap_value = "zero",  # Set tau = 0 for trimmed units
#   verbose = TRUE
# )

## ----overlap-trimming-results-------------------------------------------------
overlap_trimmed <- readRDS(system.file("extdata", "overlap_trimmed_result.rds", package = "princeBART"))
# Compare trimmed vs original PATE
cat("Original PATE:", round(pate_result$pate, 4), "\n")
cat("Trimmed PATE: ", round(overlap_trimmed$pate_trimmed, 4), "\n")
cat("Units trimmed:", overlap_trimmed$n_trimmed, "\n")

## ----sensitivity, eval = FALSE------------------------------------------------
# # Sensitivity analysis with gamma = 2
# # This computes bounds on PATE under worst-case weight perturbations
# sens <- general_BART_transportability(
#   object = pate_result,
#   gamma = 2,        # Allow weights to vary by factor of 2
#   n_sample = 100,   # Number of posterior draws to use
#   verbose = TRUE
# )
# 

## ----sensitivity_results------------------------------------------------------
sens <- readRDS(system.file("extdata", "sensitivity_result.rds", package = "princeBART"))
# View sensitivity bounds
cat("Gamma =", sens$gamma, "\n")
cat("Lower bound:", round(sens$lower$estimate, 4), 
    "[", round(sens$lower$ci[1], 4), ",", round(sens$lower$ci[2], 4), "]\n")
cat("Upper bound:", round(sens$upper$estimate, 4),
    "[", round(sens$upper$ci[1], 4), ",", round(sens$upper$ci[2], 4), "]\n")

## ----sensitivity-curve, eval = FALSE------------------------------------------
# # Sensitivity curve for multiple gamma values
# gammas <- c(1.1, 1.25, 1.5, 2, 3)
# sens_results <- lapply(gammas, function(g) {
#   general_BART_transportability(pate_result, gamma = g)
# })
# 

## ----sensitivity-curve-results------------------------------------------------
sens_results <- readRDS(system.file("extdata", "sensitivity_curve_results.rds", package = "princeBART"))
gammas <- c(1.1, 1.25, 1.5, 2, 3)
# Extract bounds
lower_bounds <- sapply(sens_results, function(x) x$lower$estimate)
upper_bounds <- sapply(sens_results, function(x) x$upper$estimate)

# Plot sensitivity curve
plot(gammas, upper_bounds, type = "l", col = "red", 
     ylim = range(c(lower_bounds, upper_bounds)),
     xlab = "Gamma", ylab = "PATE Bounds", main = "Sensitivity Analysis")
lines(gammas, lower_bounds, col = "blue")
abline(h = 0, lty = 2)
abline(h = pate_result$pate, lty = 3)
legend("topright", c("Upper", "Lower", "Null", "Point est."), 
       col = c("red", "blue", "black", "black"), lty = c(1, 1, 2, 3))

