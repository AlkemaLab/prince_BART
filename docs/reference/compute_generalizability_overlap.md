# Compute Generalizability Overlap

Estimates the selection score s = P(complier\|X, in_source) \*
P(in_source\|X) for assessing generalizability from source study to
external population. Missing covariates are single-imputed using BART
mean predictions.

## Usage

``` r
compute_generalizability_overlap(
  princebart_fit,
  newdata,
  weights = NULL,
  verbose = FALSE
)
```

## Arguments

- princebart_fit:

  A fitted `princebart` object with saved trees.

- newdata:

  External population data. Covariates missing from source model will be
  auto-detected and single-imputed using BART.

- weights:

  Survey weights for external data.

- verbose:

  Print progress.

## Value

A list with:

- `pi_c`: P(complier\|X, in_source) for each external unit

- `pi_t`: P(in_source\|X) for each external unit

- `pi_s`: Selection score s = pi_c \* pi_t

- `e_s_tilde`: Standardized selection score (logit, then z-score)

- `e_s_tilde_source`: Standardized scores for source sample
