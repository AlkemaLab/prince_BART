# Sensitivity Analysis for Unmeasured Confounding

Compute sensitivity bounds for the PATE estimate under potential
unmeasured confounding using weight-shift optimization.

## Usage

``` r
general_BART_transportability(object, gamma, n_sample = 100, verbose = FALSE)
```

## Arguments

- object:

  A `general_pate` object from [`general_BART()`](general_BART.md)

- gamma:

  Numeric sensitivity parameter \> 1. Represents the maximum ratio by
  which observation weights can be shifted. Larger values allow more
  severe confounding.

- n_sample:

  Number of posterior draws to use for sensitivity analysis. Default is
  100.

- verbose:

  Logical; print progress messages. Default is FALSE.

## Value

A list with components:

- gamma:

  The sensitivity parameter used

- lower:

  Lower bound estimate with CI

- upper:

  Upper bound estimate with CI

## Details

The sensitivity analysis follows the weight-shift framework where
observation weights can be multiplied by a factor between 1/gamma and
gamma. The lower and upper bounds represent the extremes of the PATE
estimate under this weight perturbation.

Requires the CVXR package for optimization.

## See also

[`general_BART()`](general_BART.md),
[`general_BART_overlap()`](general_BART_overlap.md)

## Examples

``` r
if (FALSE) { # \dontrun{
fit <- general_BART(princebart_fit, newdata, ...)
sens <- general_BART_transportability(fit, gamma = 1.5)
} # }
```
