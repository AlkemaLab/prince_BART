# Predict from Saved Trees

Generate predictions for new data using saved tree structures from a
prince_bart fit.

## Usage

``` r
predict_trees(trees, newdata, scaling = NULL, n_cores = 1)
```

## Arguments

- trees:

  A data.frame of tree structures from a prince_bart fit.

- newdata:

  A matrix of covariates for prediction (can be unscaled if scaling is
  provided), or a list of matrices (one per posterior sample).

- scaling:

  Optional list with `center` and `scale` named vectors for
  standardizing newdata before prediction. Typically from
  `prince_bart_fit$scaling`. If NULL, newdata is assumed to already be
  scaled.

- n_cores:

  Number of cores for parallel prediction (default: 1).

## Value

A matrix of predicted probabilities with rows for observations and
columns for posterior samples.
