# Impute Potential Outcomes for Compliers

Imputes potential outcomes Y(0) and Y(1) for compliers based on
posterior samples from a fitted prince_bart model. Used internally by
sample-based estimand functions.

## Usage

``` r
imput_potentialoutcomes_c(prince_bart_fit)
```

## Arguments

- prince_bart_fit:

  A fitted object from `prince_BART`.

## Value

A 4D array with dimensions (iteration, chain, variable, unit) containing
imputed potential outcomes:

- y0:

  Imputed Y(0) for compliers

- y1:

  Imputed Y(1) for compliers

- cy0:

  Imputed Y(0) with induced correlation

- cy1:

  Imputed Y(1) with induced correlation
