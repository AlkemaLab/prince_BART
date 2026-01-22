# Generalize Treatment Effects to External Populations (Transportability)

Estimates Population Average Treatment Effects (PATE) in an external
population using a fitted `princebart` model and external data, possibly
from complex sample surveys. Assumes that covariates X capture all
sources of effect heterogeneity, allowing the conditional complier
effect CATE_C(x) to generalize to the target population.

## Usage

``` r
general_BART(
  princebart_fit,
  newdata,
  subpop = NULL,
  psu = NULL,
  weights = NULL,
  fast_propensity = TRUE,
  n_cores = 1L,
  seed = NULL,
  verbose = FALSE
)
```

## Arguments

- princebart_fit:

  A fitted `princebart` object with `keep_trees = TRUE`.

- newdata:

  A data.frame containing the external population (e.g., survey data).
  Covariates present in the source data but missing from `newdata` will
  be automatically detected and multiply imputed using auxiliary BART
  models.

- subpop:

  Logical vector of length `nrow(newdata)` indicating which units belong
  to the target subpopulation. Default is all units.

- psu:

  Vector of primary sampling unit identifiers for `newdata`. Required
  for complex survey inference.

- weights:

  Vector of survey weights for `newdata`. Default is equal weights.

- fast_propensity:

  Logical; if TRUE (default), compute instrument propensity e = P(Z\|X)
  once using only covariates common to both source and external data,
  before imputation. This is much faster. If FALSE, compute e for each
  MI-completed dataset (slower but more statistically more faithful when
  the propensity depends strongly on imputed variables).

- n_cores:

  Number of cores for parallel computation. Default is 1.

- seed:

  Random seed for reproducibility.

- verbose:

  Logical; print progress messages. Default is FALSE

## Value

An object of class `general_pate` containing:

- `pate`: Point estimate (posterior mean) of the PATE

- `ci`: 95\\

- `sd`: Posterior standard deviation

- `draws`: Vector of posterior draws of the PATE

- `y0`: Array of predicted Y(0) values (units x iterations x chains)

- `y1`: Array of predicted Y(1) values (units x iterations x chains)

- `subpop`: Subpopulation indicator used

- `psu`: PSU identifiers used

- `weights`: Survey weights used

## Details

This function implements a multi-step procedure:

1.  Multiple imputation of missing covariates in external data using
    auxiliary BART models fit on source data.

2.  Feature expansion: compute instrument propensity e = P(Z\|X) in
    external data using BART fit on source data.

3.  Predict potential outcomes Y(0) and Y(1) using saved princeBART
    trees.

4.  Estimate PATE using Bayesian bootstrap for complex survey data.

The resulting PATE is a population-level estimand defined over the
specified target population.

The key identifying assumption is that \\\mathrm{CATE}\_C(x)\\ is
transportable, meaning that conditional on X, treatment effects for
compliers in the source study equal conditional effects in the target
population.

For sensitivity analyses (overlap trimming, confounding bounds), use
[`general_BART_overlap`](general_BART_overlap.md) and
[`general_BART_transportability`](general_BART_transportability.md) on
the returned object.

## See also

[`general_BART_overlap`](general_BART_overlap.md),
[`general_BART_transportability`](general_BART_transportability.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Fit princeBART on source study
fit <- prince_BART(Y ~ X1 + X2 + X3 | Z | W, data = source_data,
                 keep_trees = TRUE, n_samples = 1000)

# Generalize to external survey population
pate <- general_BART(
  princebart_fit = fit,
  newdata = survey_data,
  subpop = survey_data$eligible == 1,
  psu = survey_data$cluster_id,
  weights = survey_data$survey_weight,
  n_cores = 4
)

# Run sensitivity analyses on the result
overlap <- general_BART_overlap(pate, threshold = 0.05)
sens <- general_BART_transportability(pate, gamma = 2)

print(pate)
} # }
```
