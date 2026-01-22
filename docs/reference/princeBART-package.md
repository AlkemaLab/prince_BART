# princeBART: Principal Stratification with BART

The princeBART package implements principal stratification with Bayesian
Additive Regression Trees (BART) for causal inference when treatment
uptake is endogenous and identification relies on an instrument that is
as-good-as- random conditional on covariates. It targets the average
causal effect among compliers (ATE_C, also known as the Local Average
Treatment Effect, LATE), as well as conditional complier effects given
covariates (CATE_C(x)).

## Main Functions

- [`prince_BART`](prince_BART.md): Fit the principal stratification BART
  model

- [`segment_heterogeneity`](segment_heterogeneity.md): Summarize
  heterogeneous effects by segments

- [`general_BART`](general_BART.md): Generalize effects to external
  populations (PATE)

## Model

The model applies to settings with an instrument Z, an endogenous
treatment W, and an outcome Y. It requires conditional randomization of
Z given covariates X, and allows for both one-sided and two-sided
noncompliance. Three principal strata are identified:

- **Compliers**: W(1) = 1, W(0) = 0 — treatment responds to instrument

- **Never-takers**: W(1) = W(0) = 0 — never treated regardless of Z

- **Always-takers**: W(1) = W(0) = 1 — always treated regardless of Z

BART is used to flexibly model both the stratum membership probabilities
and the potential outcomes within each stratum, allowing for
heterogeneous conditional complier effects CATE_C(x) across covariate
values.

## Author

**Maintainer**: Lucas Godoy Garraza <lgodoygarraz@umass.edu>
([ORCID](https://orcid.org/0000-0001-8426-9485))
