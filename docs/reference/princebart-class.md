# S3 Methods for princebart Objects

Print, summary, and coef methods for objects returned by `prince_BART`.

## Usage

``` r
# S3 method for class 'princebart'
print(x, ...)

# S3 method for class 'princebart'
summary(object, type = c("mixed", "sample"), treated_only = FALSE, ...)

# S3 method for class 'princebart'
coef(object, type = c("mixed", "sample"), treated_only = FALSE, ...)
```

## Arguments

- x, object:

  A princebart object.

- ...:

  Additional arguments (currently ignored).

- type:

  Character; type of estimand: `"mixed"` (default) for mixed estimands
  averaging over principal stratum uncertainty, or `"sample"` for sample
  estimands using imputed principal strata (experimental).

- treated_only:

  Logical; if `TRUE`, compute treatment effect only among the treated
  units (ATT). Default is `FALSE` (ATE).

## Value

`print`: Invisibly returns the object. `summary`: Prints and invisibly
returns a list with strata distribution and treatment effect estimates.
`coef`: Posterior summary of the requested estimand.

## Examples

``` r
if (FALSE) { # \dontrun{
fit <- prince_BART(Y ~ X | Z | W, data = df)

# Mixed ATE for compliers (default)
summary(fit)
coef(fit)

# Mixed ATT for compliers
summary(fit, treated_only = TRUE)
coef(fit, treated_only = TRUE)

# Sample estimands (experimental)
coef(fit, type = "sample")
} # }
```
