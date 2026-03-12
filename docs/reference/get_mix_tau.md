# Compute Mixed Treatment Effects

Internal function to compute mixture-based estimands from posterior
probability arrays.

## Usage

``` r
get_mix_tau(p_arr, treated = NULL)
```

## Arguments

- p_arr:

  4D probability array from prince_bart fit.

- treated:

  Optional logical vector indicating treated units.

## Value

List with strata probabilities and mean effects as draws arrays.
