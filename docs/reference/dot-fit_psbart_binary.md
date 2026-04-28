# Fit Principal Stratification BART (Single Chain)

Internal function that fits a single chain of the principal
stratification BART model. Users should typically use
[`prince_BART`](prince_BART.md) instead, which handles data preparation
and runs multiple chains.

## Usage

``` r
.fit_psbart_binary(
  X,
  Y,
  Z,
  W,
  n_warmup = 1000L,
  n_samples = 1000L,
  save_trees = FALSE,
  k = 2,
  n_trees = 200L,
  verbose = FALSE,
  n_initial = 0L
)
```

## Arguments

- X:

  A scaled numeric matrix of covariates with propensity `e` already
  appended as the last column.

- Y:

  A binary outcome vector (0/1).

- Z:

  A binary treatment assignment/instrument vector (0/1).

- W:

  A binary treatment uptake/received vector (0/1).

- n_warmup:

  Number of warmup/burn-in iterations (default: 1000).

- n_samples:

  Number of posterior samples to collect (default: 1000).

- save_trees:

  Logical; save tree structures for prediction (default: FALSE).

- k:

  Prior hyperparameter for node mean prior (default: 2).

- n_trees:

  Number of trees in the BART ensemble (default: 200).

- verbose:

  Logical; print progress (default: FALSE).

- n_initial:

  Number of initial iterations using MoM offsets (internal, default: 0).

## Value

A list containing:

- imputed:

  Array of imputed compliance class memberships with dimensions
  (n_samples, n, 2) and named variables: "nt" (never-takers), "at"
  (always-takers)

- probs:

  Array of posterior probabilities and outcome means with dimensions
  (n_samples, n, 6) and named variables: "p_a" (P(always-taker)), "p_n"
  (P(never-taker)), "m_y0c" (\\E\[Y(0) \mid complier\]\\), "m_y1c"
  (\\E\[Y(1) \mid complier\]\\), "m_y0n" (\\E\[Y(0) \mid
  never-taker\]\\), "m_y1a" (\\E\[Y(1) \mid always-taker\]\\)

- trees:

  Tree structures if `save_trees = TRUE`
