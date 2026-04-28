# Fit Principal Stratification BART for Ordinal Uptake (Single Chain)

Internal function that fits one MCMC chain of the ordinal principal
stratification BART model. Users should typically call
[`prince_BART`](prince_BART.md).

## Usage

``` r
.fit_psbart_ordinal(
  data = NULL,
  n_warmup = 20L,
  n_samples = 20L,
  k = 2,
  n_trees = 200L,
  lambda_z1 = 1,
  lambda_z0 = 1,
  init_w_poisson_lambda = 3,
  n_thin = 1L,
  n_threads = 4L,
  monotonicity = TRUE,
  rho = 0,
  save_trees = TRUE,
  save_trees_interval = 10L,
  verbose = FALSE,
  ...
)
```

## Arguments

- data:

  A data.frame/data.table containing columns `Y`, `Z`, `W`, and
  covariates.

- n_warmup:

  Number of warmup iterations.

- n_samples:

  Number of posterior samples to retain.

- k:

  BART node prior hyperparameter.

- n_trees:

  Number of trees per BART ensemble.

- lambda_z1:

  Box-Cox parameter for latent \\z_1\\.

- lambda_z0:

  Box-Cox parameter for latent \\z_0\\.

- init_w_poisson_lambda:

  Poisson rate used to jitter initial `w0`/`w1` assignments.

- n_thin:

  Thinning interval for retained posterior samples.

- n_threads:

  Number of threads used by dbarts.

- monotonicity:

  Logical; if TRUE enforce \\w_1 \le w_0\\ when constructing candidate
  strata.

- rho:

  Correlation used in bivariate latent sampling.

- save_trees:

  Logical; if TRUE include sampled trees.

- save_trees_interval:

  Save trees every this many retained samples.

- verbose:

  Logical; print progress updates.

- ...:

  Unused.

## Value

A list with elements:

- imputed:

  Array `(n_samples, n, 2)` for sampled `w0`, `w1`.

- probs:

  Array `(n_samples, n, 2)` for `m_y0`, `m_y1`.

- check:

  Array `(n_samples, n, 2)` with posterior predictive draws for `w0`,
  `w1` diagnostics.

- trees:

  Tree structures when `save_trees = TRUE`.
