# Treatment Effect Estimands for Compliers

Internal functions to compute various treatment effect estimands for the
complier stratum from a fitted prince_bart model.

## Usage

``` r
mate_c(prince_bart_fit)

matt_c(prince_bart_fit)

sate_c(prince_bart_fit, induce_residual_corr = FALSE)

satt_c(prince_bart_fit, induce_residual_corr = FALSE)
```

## Arguments

- prince_bart_fit:

  A fitted object from `prince_BART`.

- induce_residual_corr:

  Logical; for sample estimands, whether to induce residual correlation
  between potential outcomes (default: FALSE).

## Value

A `posterior` summary object with posterior mean, median, standard
deviation, and quantiles.

## Functions

- `mate_c()`: Mixed Average Treatment Effect for Compliers (MATE_C)

- `matt_c()`: Mixed Average Treatment Effect on the Treated for
  Compliers (MATT_C)

- `sate_c()`: Sample Average Treatment Effect for Compliers (SATE_C)

- `satt_c()`: Sample Average Treatment Effect on the Treated for
  Compliers (SATT_C)
