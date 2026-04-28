# Principal Stratification using BART

Fits a Bayesian principal stratification model using Bayesian Additive
Regression Trees (BART) for causal inference with endogenous treatments.
The method is designed for instrumental variable or encouragement
designs with noncompliance and targets causal effects for compliers.

## Usage

``` r
prince_BART(
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
  uptake_type = c("auto", "binary", "ordinal"),
  rho = 0,
  verbose = FALSE
)
```

## Arguments

- formula:

  A formula following 2SLS convention: `Y ~ X1 + X2 | Z | W`, where Y is
  the outcome, X1 + X2 are covariates, Z is the instrument (treatment
  assignment), and W is treatment uptake.

- data:

  A data.frame containing the variables in the formula.

- X:

  A matrix or data.frame of covariates. Used when `formula = NULL`.

- Y:

  A binary outcome vector (0/1). Used when `formula = NULL`.

- Z:

  A binary instrument or treatment assignment vector (0/1).

- W:

  Treatment uptake/received vector.

- propensity:

  Optional pre-computed instrument propensity scores \\e = P(Z \mid
  X)\\. If NULL (default), propensity scores are estimated internally
  using BART.

- instrument_overlap:

  Optional bounds for instrument propensity trimming to enforce overlap.
  If provided as a length-2 vector (e.g., `c(0.1, 0.9)`), observations
  with estimated propensity scores outside this range are excluded.
  Default is NULL (no trimming).

- n_warmup:

  Number of warmup iterations per chain (default: 1000).

- n_samples:

  Number of posterior samples per chain (default: 1000).

- n_chains:

  Number of parallel MCMC chains (default: 4).

- keep_trees:

  Logical; save fitted BART tree structures for downstream prediction
  and generalization (default: FALSE).

- k:

  Prior hyperparameter controlling node shrinkage in BART (default: 2).

- n_trees:

  Number of trees in each BART ensemble (default: 200).

- workers:

  Number of parallel workers. If NULL (default), uses all available
  cores up to `n_chains`. Set to 1 for sequential execution.

- uptake_type:

  Character scalar controlling uptake model: `"auto"` (default) chooses
  binary if `W in {0,1}` and ordinal otherwise; `"binary"` enforces
  binary uptake model; `"ordinal"` enforces ordinal/count uptake model.

- rho:

  Numeric in (-1, 1); residual correlation between the latent propensity
  scores for \\W(1)\\ and \\W(0)\\ in the bivariate ordinal sampler.
  Only used when `uptake_type = "ordinal"`. Default is `0`
  (conditionally independent potential treatments given X). Varying this
  parameter is a natural sensitivity analysis since the assumption is
  not testable from the observed data.

- verbose:

  Logical; print progress messages (default: FALSE).

## Value

A list of class `"prince_bart"` containing posterior draws from all
chains, including:

- `imp`: Imputed latent quantities (iteration x chain x variable x
  unit).

- `probs`: Posterior draws of probabilities/outcome means.

- `check`: Ordinal diagnostic array (ordinal mode only).

- `trees`: Fitted BART trees (if `keep_trees = TRUE`).

- `data`: Stored input data, including `X_model` (processed covariates
  used by the fitted model), `X_raw` (raw covariates when available),
  `X` (compatibility alias to `X_model`), instrument propensity scores,
  and outcomes.

## Details

Shared preprocessing is applied in all modes: formula parsing,
validation, covariate scaling, propensity estimation for the binary
instrument Z, optional overlap trimming, and propensity augmentation of
X. After preprocessing, model fitting dispatches to either binary or
ordinal PS-BART chain runners based on `uptake_type`.

For fits created from a formula or data.frame input, the returned object
stores both a model-space covariate matrix and, when available, a raw
covariate data.frame. The model-space representation is used internally
for fitting and prediction; the raw representation is retained for
downstream tasks that benefit from original factor/ordered-factor
classes, such as [`segment_heterogeneity()`](segment_heterogeneity.md).

The model jointly estimates:

- Principal stratum membership probabilities (\\P(\text{complier} \mid
  X)\\, \\P(\text{never-taker} \mid X)\\, \\P(\text{always-taker} \mid
  X)\\),

- Potential outcome regressions within strata, including \\E\[Y(0) \mid
  \text{complier}, X\]\\ and \\E\[Y(1) \mid \text{complier}, X\]\\.

The primary identified causal estimand is the conditional average
treatment effect among compliers, \\\mathrm{CATE}\_C(x)\\. Mixed
(sample-based) averages of these conditional effects, including
LATE-like estimands, can be obtained using
[`summary()`](https://rdrr.io/r/base/summary.html) and
[`coef()`](https://rdrr.io/r/stats/coef.html).

Identification relies on standard instrumental variable assumptions:
conditional independence of the instrument given X, exclusion
restriction, and monotonicity (no defiers).

Parallel computation is handled via
[`future.apply::future_lapply`](https://future.apply.futureverse.org/reference/future_lapply.html).
Before calling this function, set your preferred parallel backend:

    # For Unix/Mac (forked processes)
    future::plan(future::multicore)

    # For Windows or any platform
    future::plan(future::multisession)

## See also

[`segment_heterogeneity`](segment_heterogeneity.md),
[`general_BART`](general_BART.md)

## Examples

``` r
if (FALSE) { # \dontrun{
library(future)
plan(multisession, workers = 4)

fit <- prince_BART(
  Y ~ X1 + X2 + X3 | Z | W,
  data = mydata,
  uptake_type = "auto",
  n_chains = 4,
  n_warmup = 1000,
  n_samples = 1500
)

summary(fit)
} # }
```
