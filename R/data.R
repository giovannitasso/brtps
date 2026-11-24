#' Simulated spatial data
#'
#' A small simulated data set for examples and tests: 150 locations on the
#' unit square with a Gaussian random field of Matern covariance
#' (\eqn{\nu = 1}, range \eqn{\rho = 0.3}, standard deviation 1), as in
#' Simulation 2 of Cavieres and Krumscheid (2026), plus two covariates. The
#' response is \eqn{y = 1 + 0.5 x_1 - 0.8 x_2 + u(s) + \varepsilon} with
#' \eqn{\varepsilon \sim N(0, 0.3^2)}.
#'
#' @format A data frame with 150 rows and 6 columns:
#' \describe{
#'   \item{s1, s2}{Coordinates.}
#'   \item{x1}{Standard normal covariate.}
#'   \item{x2}{Uniform covariate on (-1, 1).}
#'   \item{y}{Response.}
#'   \item{u}{True value of the spatial field.}
#' }
#' @source Simulated; see `data-raw/brtps_sim.R` in the source repository.
"brtps_sim"
