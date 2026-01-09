# princeBART

**princeBART** is an R package implementing **Principal Stratification with Bayesian Additive Regression Trees (BART)** for causal inference with endogenous treatments.

The package is designed for instrumental variable and encouragement designs with noncompliance, and provides Bayesian estimation of complier-specific causal effects, effect heterogeneity, and generalization to external populations.

The methodology is described in:

Godoy Garraza, L., Speizer, I, Alkema, L. 
*Combining BART and Principal Stratification to estimate the effect of intermediate variables on primary outcomes with application to estimating the effect of family planning on employment in Nigeria and Senegal.*  
(working paper)

Repository: [https://github.com/AlkemaLab/prince_BART](https://github.com/AlkemaLab/prince_BART)

---

## Installation

The package is currently available from GitHub:

```r
# install.packages("devtools")
devtools::install_github("AlkemaLab/prince_BART")
```

---

## Overview

princeBART estimates causal effects when:
- Treatment uptake is endogenous
- Identification relies on an instrument that may only be valid after controlling for observed covariates
- Treatment effects may vary flexibly with covariates

**Key features include:**
- Bayesian principal stratification using BART
- Estimation of average effects among compliers (LATE-like estimands)
- Conditional complier effects and mixed (sample-based) averages
- Segment-based heterogeneity summaries using shallow trees
- Generalization of treatment effects to external populations
- Support for complex survey designs, overlap diagnostics, and sensitivity analyses

---

## Basic Usage

```r
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

---

## Vignettes

The package includes worked examples illustrating core functionality.  
Browse vignettes with:

```r
browseVignettes("princeBART")
```

**Available vignettes include:**
- [Introduction to princeBART](https://htmlpreview.github.io/?https://github.com/AlkemaLab/prince_BART/blob/main/doc/introduction.html)
- [Estimation of complier effects and effect heterogeneity](https://htmlpreview.github.io/?https://github.com/AlkemaLab/prince_BART/blob/main/doc/simulation.html)
- [Generalizing princeBART results to external populations, including diagnostics](https://htmlpreview.github.io/?https://github.com/AlkemaLab/prince_BART/blob/main/doc/generalizing.html)


---

## Related Materials

- Paper preprint: [https://arxiv.org/abs/2408.03777](https://arxiv.org/abs/2408.03777)
- Data and documentation: [https://doi.org/10.15139/S3/BRLE7L](https://doi.org/10.15139/S3/BRLE7L)

---

## License

MIT License.
