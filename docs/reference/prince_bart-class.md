# S3 Methods for prince_bart Objects

Print, summary, and coef methods for objects returned by `prince_BART`.
Behavior is specialized by modality: binary fits report principal strata
and treatment effects for compliers, while ordinal fits report contrasts
for the monotone compliance group and stratified by baseline uptake
level.

## Usage

``` r
# S3 method for class 'prince_bart'
print(x, ...)

# S3 method for class 'prince_bart'
summary(
  object,
  type = c("mixed", "sample"),
  treated_only = FALSE,
  adaptive_levels = TRUE,
  cumulative_mass = 0.8,
  level_threshold = NULL,
  ...
)

# S3 method for class 'prince_bart'
coef(
  object,
  type = c("mixed", "sample"),
  treated_only = FALSE,
  adaptive_levels = TRUE,
  cumulative_mass = 0.8,
  level_threshold = NULL,
  ...
)
```

## Arguments

- x, object:

  A prince_bart object (either prince_bart_binary or
  prince_bart_ordinal).

- ...:

  Additional arguments (currently ignored).

- type:

  Character; type of estimand: `"mixed"` (default) for mixed estimands
  averaging over principal stratum uncertainty, or `"sample"` for sample
  estimands using imputed principal strata (experimental). For ordinal
  fits, this parameter is accepted but may be ignored depending on
  implementation.

- treated_only:

  Logical; if `TRUE`, compute treatment effect only among the treated
  units (ATT). Default is `FALSE` (ATE). For ordinal fits, this
  parameter is accepted but may be ignored or handled differently.

- adaptive_levels:

  Logical; for ordinal fits, if `TRUE` (default), show level-specific
  effects up to a data-adaptive threshold based on cumulative
  affected-unit mass.

- cumulative_mass:

  Numeric in (0, 1\]; for ordinal fits with adaptive grouping, levels
  are shown individually until this cumulative mass is reached. Default
  is `0.80`.

- level_threshold:

  Optional integer threshold `K` for ordinal fits. If supplied, levels
  `1..K` are shown individually and higher levels are pooled as
  `"Level > K"`. This overrides adaptive grouping.

## Value

`print`: Invisibly returns the object. `summary`: Prints and invisibly
returns a list with posterior summaries (including uncertainty and
diagnostics) appropriate to the modality. `coef`: Named numeric vector
of posterior mean estimates for the key estimands.

## Examples

``` r
if (FALSE) { # \dontrun{
# Binary fit
fit_binary <- prince_BART(Y ~ X | Z | W, data = df)
summary(fit_binary)
coef(fit_binary)

# Ordinal fit
fit_ordinal <- prince_BART(Y ~ X | Z | W, data = df, uptake_type = "ordinal")
summary(fit_ordinal)
coef(fit_ordinal)
} # }
```
