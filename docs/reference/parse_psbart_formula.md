# Parse Principal Stratification Formula

Parses a formula in 2SLS style: Y ~ X1 + X2 \| Z \| W

## Usage

``` r
parse_psbart_formula(formula, data)
```

## Arguments

- formula:

  A formula with three parts separated by \|

- data:

  A data.frame containing the variables

## Value

A list with components X, Y, Z, W
