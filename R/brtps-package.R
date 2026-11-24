#' brtps: Bayesian Regularized Thin Plate Splines
#'
#' Bayesian spatial models in which the spatial random field is approximated
#' with the regularized thin plate spline (regTPS) approach of Cavieres and
#' Krumscheid (2026). The main function is [brtps()]; the spatial component is
#' specified with [regtps_s()] and priors with [brtps_priors()].
#'
#' @section How the package is organised:
#' * `R/regtps_s.R` and `src/brtps/regtps_s.hpp`: the spatial model (the R
#'   side builds the basis, the C++ side evaluates the posterior).
#' * `R/families.R` and `src/brtps/families.hpp`: likelihoods.
#' * `R/priors.R` and `src/brtps/priors.hpp`: prior densities.
#' * `R/engines.R`: samplers ('SparseNUTS', 'tmbstan', Laplace).
#' * `R/methods.R`: summaries, predictions and draws.
#'
#' @references Cavieres, J. and Krumscheid, S. (2026). regTPS: a novel
#'   approach to approximate a Gaussian random field for Bayesian spatial
#'   modeling. *Statistics and Computing*, 36, 261.
#'   \doi{10.1007/s11222-026-10990-x}
#'
#' @keywords internal
"_PACKAGE"

## usethis namespace: start
#' @useDynLib brtps, .registration = TRUE
#' @importFrom stats coef complete.cases delete.response dnorm fitted gaussian
#'   median model.frame model.matrix model.offset model.response na.pass nobs
#'   predict quantile residuals rnorm sd terms .getXlevels
## usethis namespace: end
NULL
