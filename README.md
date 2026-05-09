# princeBART

**princeBART** is an R package implementing **Principal Stratification with Bayesian Additive Regression Trees (BART)** for causal inference with endogenous treatments.

The package is designed for instrumental variable and encouragement designs with noncompliance, and provides Bayesian estimation of complier-specific causal effects, effect heterogeneity, and tools for assessing and implementing generalization to external populations under explicit transportability assumptions, including sensitivity analyses.

The methodology is described in:

Godoy Garraza, L., Speizer, I, Alkema, L. (forthcoming). Combining BART and Principal Stratification to estimate the effect of intermediate variables on primary outcomes with application to estimating the effect of family planning on employment in Nigeria and Senegal. Statistics and Public Policy. See <https://arxiv.org/abs/2412.16320>.

Godoy Garraza, L., Speizer, I, Alkema, L. (2024). How to estimate causal effects associated with family planning? An introduction to Prince BART, a new approach to effect estimation based on principal stratification and Bayesian non-parametric models. See <https://verixiv.org/articles/1-5>.

Godoy Garraza, L., Speizer, I, Alkema, L. (2025). Does fertility affect woman's labor force participation in low- and middle-income settings? Findings from a Bayesian nonparametric analysis. See <https://arxiv.org/abs/2508.10787>.

R package repository: <https://github.com/AlkemaLab/prince_BART>

R Package website: <https://alkemalab.github.io/prince_BART>

------------------------------------------------------------------------

## Installation

The package is currently available from GitHub:

``` r
# install.packages("devtools")
devtools::install_github("AlkemaLab/prince_BART")
```

------------------------------------------------------------------------

## Overview

princeBART estimates causal effects when:

-   Treatment uptake is endogenous

-   Identification relies on an instrument that may only be valid after controlling for observed covariates

-   Treatment effects may vary flexibly with covariates

Key features include:

-   Bayesian principal stratification using BART

-   Estimation of average effects among compliers (LATE-like estimands)

-   Conditional complier effects and mixed (sample-based) averages

-   Segment-based heterogeneity summaries using shallow trees

-   Generalization of treatment effects to external populations

-   Support for complex survey designs, overlap diagnostics, and sensitivity analyses

------------------------------------------------------------------------

## Basic Usage

``` r
library(princeBART)

fit <- prince_BART(
  Y ~ age + education + income | Z | W,
  data = mydata,
  n_chains = 4,
  n_warmup = 1000,
  n_samples = 1000
)

summary(fit)
coef(fit)
```

------------------------------------------------------------------------

## Vignettes

-   [**Introduction to princeBART**](https://alkemalab.github.io/prince_BART/articles/introduction.html)\
    Estimation of complier-specific causal effects and summaries of effect heterogeneity.

-   [**Generalizing princeBART results to external populations**](https://alkemalab.github.io/prince_BART/articles/generalizing.html)\
    Transporting estimated effects to a target population using external data, with overlap diagnostics and sensitivity analysis.

-   [**Simulation study**](https://alkemalab.github.io/prince_BART/articles/simulation.html)\
    Comparison of princeBART with standard estimators such as two-stage least squares (2SLS) in simulated settings.

-   [**Prince BART with Ordinal Uptake**](https://alkemalab.github.io/prince_BART/articles/ordinal.html)\
    Extension to treatments with multiple ordered levels, illustrated with a fertility application estimating effects of number of children on employment.

------------------------------------------------------------------------

## Related Materials

-   Data and documentation of MLE data, used in application: <https://doi.org/10.15139/S3/BRLE7L>

------------------------------------------------------------------------

## License

MIT License.
