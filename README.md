# brtps: Bayesian Regularized Thin Plate Splines

R package for Bayesian spatial modelling with **regTPS** (Cavieres &
Krumscheid, 2026): the spatial random field is approximated by a truncated
Karhunen–Loève expansion built from the thin plate spline penalty. Models are
written in [TMB](https://github.com/kaskr/adcomp) and sampled with
[SparseNUTS](https://noaa-afsc.github.io/SparseNUTS/).

> **Status:** package skeleton. The first model (spatial, Gaussian response)
> is being added.

## Reference

Cavieres, J. and Krumscheid, S. (2026). regTPS: a novel approach to
approximate a Gaussian random field for Bayesian spatial modeling.
*Statistics and Computing*, 36, 261.
<https://doi.org/10.1007/s11222-026-10990-x>
