# Ordinal Complier-Like Contrast Estimands

Compute posterior summaries for contrasts in ordinal uptake settings.
Focuses on the monotone compliance group (W(0) - W(1) = 1) and
stratifies by baseline uptake W(0) up to level 5.

## Usage

``` r
estimands_ordinal_mixed(
  prince_bart_fit,
  adaptive_levels = TRUE,
  cumulative_mass = 0.8,
  level_threshold = NULL
)
```

## Arguments

- prince_bart_fit:

  A fitted ordinal prince_bart object.

- adaptive_levels:

  Logical; if TRUE (default), choose reported levels adaptively by
  cumulative affected-unit mass.

- cumulative_mass:

  Numeric in (0, 1\]; target cumulative mass used when
  \codeadaptive_levels = TRUE. Default is \code0.80.

- level_threshold:

  Optional integer threshold K. If supplied, levels \code1..K are shown
  individually and higher levels are pooled as \code"Level \> K"; this
  overrides adaptive selection.

## Value

A list with the following elements:

- overall:

  Data frame of posterior summary for \eqnE\[Y(w=1)-Y(w=0) \mid
  W(0)-W(1)=1\].

- by_w0:

  Data frame of posterior summaries by reported level groups.

- grouping:

  List with grouping metadata used for level reporting.
