# Compute Sample-Based Treatment Effects

Internal function to compute sample-based estimands from imputed
compliance classes and outcomes.

## Usage

``` r
get_sample_tau(imp_g, imp_o, treated = NULL, include_corr = TRUE)
```

## Arguments

- imp_g:

  Imputed compliance class array.

- imp_o:

  Imputed potential outcomes array.

- treated:

  Optional logical vector indicating treated units.

- include_corr:

  Logical; use correlated imputations.

## Value

A draws array with posterior samples.
