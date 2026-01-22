# Internal Binary BART Wrapper

A custom wrapper for dbarts that properly handles binary outcomes,
including edge cases where a stratum may have only 0s or only 1s.

## Usage

``` r
dbarts_binary(
  x,
  y,
  test = NULL,
  subset = NULL,
  offset = NULL,
  control = dbarts::dbartsControl(),
  k = 2
)
```

## Arguments

- x:

  Predictor matrix.

- y:

  Response vector (binary 0/1).

- test:

  Test predictor matrix.

- subset:

  Logical vector for subsetting training data.

- offset:

  Offset values for training data.

- control:

  A dbartsControl object.

- k:

  Prior hyperparameter for node prior (default 2).

## Value

A dbartsSampler object configured for binary probit BART.

## Details

This addresses a limitation in dbarts::dbarts where binary detection
only triggers when there are exactly 2 unique values (0 and 1). In
principal stratification, subsets may sometimes contain only one unique
outcome value, which would cause dbarts to treat it as continuous and
fail when estimating sigma.

Based on the fix from princeB.r lines 679-739.
