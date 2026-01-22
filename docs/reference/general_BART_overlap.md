# Generalizability Overlap Analysis for PATE Estimates

Compute generalizability overlap scores and optionally trim observations
with low overlap to produce a trimmed PATE estimate.

## Usage

``` r
general_BART_overlap(
  object,
  threshold = NULL,
  overlap_value = c("zero", "NA"),
  verbose = FALSE
)
```

## Arguments

- object:

  A `general_pate` object from [`general_BART()`](general_BART.md)

- threshold:

  Numeric threshold for trimming based on standardized selection scores.
  Observations with `|e_s_tilde| > threshold` are excluded. Default is
  NULL (no trimming, just compute overlap).

- overlap_value:

  Value to replace `tau` with for trimmed observations. Either `"zero"`
  (set to 0) or `"NA"` (exclude completely). Default is `"zero"`.

- verbose:

  Logical; print progress messages. Default is FALSE.

## Value

A list with components:

- overlap:

  Data frame with overlap metrics: pi_c, pi_t, pi_s, e_s_tilde

- e_s_tilde_source:

  Standardized selection scores for source data

- n_trimmed:

  Number of observations trimmed (if threshold used)

- pate_trimmed:

  Trimmed PATE estimate (if threshold used)

- ci_trimmed:

  95 percent CI for trimmed PATE (if threshold used)

## See also

[`general_BART()`](general_BART.md),
[`general_BART_transportability()`](general_BART_transportability.md)

## Examples

``` r
if (FALSE) { # \dontrun{
fit <- general_BART(princebart_fit, newdata, ...)
overlap <- general_BART_overlap(fit, threshold = 2)
} # }
```
