# brtps: Bayesian Regularized Thin Plate Splines

<!-- badges: start -->
[![R-CMD-check](https://github.com/giovannitasso/brtps/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/giovannitasso/brtps/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

`brtps` fits Bayesian spatial regression models in which the spatial random
field is approximated with **regTPS** (Cavieres & Krumscheid, 2026): a
truncated Karhunen–Loève expansion built from the eigendecomposition of the
thin plate spline penalty. Models are written in
[TMB](https://github.com/kaskr/adcomp) and sampled with the sparse no-U-turn
sampler of [SparseNUTS](https://noaa-afsc.github.io/SparseNUTS/).

> **Status:** early development. The spatial model with a Gaussian response
> works; other families (Poisson, binomial, Gamma, lognormal) and the
> spatio-temporal model (`regtps_st()`) are planned.

## Installation

The package contains C++ code, so on Windows you need
[Rtools](https://cran.r-project.org/bin/windows/Rtools/) installed first.

```r
# install.packages("remotes")
remotes::install_github("giovannitasso/brtps")

# Sampler used by default (optional for the Laplace approximation):
install.packages("StanEstimators",
                 repos = c("https://andrjohns.r-universe.dev", "https://cloud.r-project.org"))
remotes::install_github("noaa-afsc/SparseNUTS")
install.packages(c("future", "callr"))   # to run chains in parallel
```

## Example

```r
library(brtps)

# Quick check with the Laplace approximation (seconds, no MCMC)
fit0 <- brtps(y ~ x1 + x2, data = brtps_sim,
              model = regtps_s(s1, s2, k_s = 100, m_s_trunc = 0.99),
              engine = "laplace")
fit0

# Full Bayesian fit with sparse NUTS
options(mc.cores = 4)
fit <- brtps(y ~ x1 + x2, data = brtps_sim,
             model = regtps_s(s1, s2, k_s = 100, m_s_trunc = 0.99),
             chains = 4, num_samples = 1000)
summary(fit)

# Map of the spatial field
grid <- expand.grid(s1 = seq(0, 1, length.out = 50), s2 = seq(0, 1, length.out = 50))
field <- predict(fit, newdata = grid, type = "spatial")

# Diagnostics with posterior / bayesplot, and PSIS-LOO
draws <- posterior::as_draws_array(fit, variable = c("x1", "x2", "sigma", "alpha"))
bayesplot::mcmc_trace(draws)
loo::loo(fit)
```

### Main arguments

| Argument | Meaning |
|---|---|
| `regtps_s(s1, s2, ...)` | Coordinates of the spatial field (column names in `data`). |
| `k_s` | Number of thin plate spline basis functions (`K` in the paper). |
| `m_s_trunc` | Truncation: a proportion of prior variance to keep (e.g. `0.99`) or a number of components `M`. |
| `priors` | `brtps_priors(sigma = prior_pc(0.5, 0.05), alpha = prior_lognormal(0, 1), beta = ...)`. |
| `engine` | `"snuts"` (default), `"tmbstan"` or `"laplace"`. |

## How the package is organised

```
R/regtps_s.R          spatial component: TPS basis, eigendecomposition, truncation
R/brtps.R             brtps(): formula -> TMB object -> engine -> fit
R/families.R          observation models (one entry per family)
R/priors.R            prior constructors and brtps_priors()
R/engines.R           SparseNUTS, tmbstan and Laplace
R/methods.R           summary, predict, coef, loo, posterior draws
src/brtps.cpp         one TMB library; dispatches on the model name
src/brtps/*.hpp       regtps_s model, families and priors in C++
tests/testthat/       tests (including the truncation levels of the paper)
```

Adding a family means one entry in `R/families.R` and one `case` in
`src/brtps/families.hpp`. Adding a model (such as `regtps_st`) means a
constructor in `R/`, a header in `src/brtps/` and one line in the dispatcher
of `src/brtps.cpp`.

## Reference

Cavieres, J. and Krumscheid, S. (2026). regTPS: a novel approach to
approximate a Gaussian random field for Bayesian spatial modeling.
*Statistics and Computing*, 36, 261.
<https://doi.org/10.1007/s11222-026-10990-x>
