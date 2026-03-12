# Effect Heterogeneity by Segments

Partition the sample into covariate-defined segments and summarize how
estimated latent-group effects vary across segments, following the
heterogeneity analysis in the Prince BART paper.

## Usage

``` r
segment_heterogeneity(
  princebart,
  data = NULL,
  vars = NULL,
  rpart_control = rpart::rpart.control(cp = 0, maxdepth = 3L),
  min_compliers_bucket = 100,
  ci_level = 0.9,
  plot = FALSE,
  contrast = FALSE
)
```

## Arguments

- princebart:

  A fitted object from `prince_BART`.

- data:

  Optional data.frame of covariates used for segmentation. If NULL, uses
  stored raw covariates when available, otherwise `princebart$data$X`.

- vars:

  Optional character vector of covariate names to include in the
  segmentation model. If NULL, uses all covariates in `data`.

- rpart_control:

  An `rpart.control` object for tree fitting.

- min_compliers_bucket:

  Minimum effective latent-group count for pruning (sum of posterior
  mean latent-group membership probabilities within a segment). Default
  is 100.

- ci_level:

  Credible interval level (e.g., 0.9 for 90% credible intervals).
  Default is 0.9.

- plot:

  Logical; if TRUE, return ggplot2 objects for segment effects (and the
  max–min contrast plot if `contrast = TRUE`).

- contrast:

  Logical; if TRUE, compute posterior draws for the difference between
  the segments with the highest and lowest estimated mean effects.

## Value

A list with components:

- tree:

  Fitted (and optionally pruned) rpart tree.

- segment_data:

  Input data with `cate` (posterior mean unit-level conditional effect),
  `w` (posterior mean latent-group membership probability), and
  `segment`.

- effects:

  Segment-level effect summaries: posterior mean, sd, credible interval
  bounds, and `p_gt0`, plus `n` (segment size) and `n_group` (estimated
  weighted latent-group size).

- draws:

  List of posterior draw matrices by segment.

- contrast:

  If `contrast = TRUE` , posterior comparison of the highest- vs
  lowest-effect segments, containing `$summary` and posterior draws in
  `$draws`.

- plot:

  If `plot = TRUE` , a list of ggplot objects including `$effect`
  (segment effects) and, when `contrast = TRUE` , `$diff` (difference
  distribution).

## Details

This function fits a shallow `rpart` regression tree to posterior mean
conditional latent-group effects and assigns each unit to a
terminal-node segment. Tree fitting is weighted by each unit's posterior
mean latent-group membership probability, so that splits emphasize
regions of the covariate space with more relevant units.

When raw covariates were stored in the fitted object, they are used by
default so that factor and ordered-factor variables can be represented
more interpretably in segment labels. Legacy fit objects may only store
processed covariates; in that case, users can supply `data =` explicitly
to recover more interpretable labels.

For each segment, the function aggregates posterior draws of conditional
effects using latent-group membership weights to obtain a
segment-specific average effect for the relevant latent group. For
binary fits, this reproduces segment-specific mixed effects among
compliers. For ordinal fits, this targets segment-specific mixed effects
among affected units defined by \\W(0)-W(1)=1\\. These segment-level
summaries are mixed (sample-based) estimands: they average conditional
effects over the empirical covariate distribution within each segment
rather than defining new population-level causal parameters.

The main segment summaries are returned in `res$effects` and visualized
in `res$plot$effect` (if `plot = TRUE`).

Optionally, segments with small effective latent-group counts can be
pruned/merged for stability. If `contrast = TRUE`, heterogeneity is
summarized by the posterior distribution of the difference between the
segments with the largest and smallest estimated mean effects. Numerical
results are in `res$contrast$summary` and the histogram is in
`res$plot$diff` (if `plot = TRUE`).
