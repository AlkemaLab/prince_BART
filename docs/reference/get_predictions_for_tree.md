# Get Predictions for a Single Tree

Recursive function to traverse a single BART tree and generate
predictions.

## Usage

``` r
get_predictions_for_tree(tree, x)
```

## Arguments

- tree:

  A data.frame representing a single tree structure.

- x:

  A matrix of covariates.

## Value

A numeric vector of predictions for each row in x.
