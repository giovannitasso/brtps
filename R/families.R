# Registry of observation models.
#
# Each entry describes one family/link combination: the integer code used by
# src/brtps/families.hpp, the names of its dispersion parameters (estimated
# on the log scale), their default priors, and R functions used after fitting
# (inverse link, pointwise log-likelihood and simulation). To add a family,
# add an entry here and the matching case in families.hpp.
brtps_family_registry <- function() {
  list(
    gaussian = list(
      family = "gaussian",
      link = "identity",
      code = 0L,
      phi_names = "sigma",
      default_priors = function() list(sigma = prior_pc(u = 0.5, prob = 0.05)),
      linkinv = function(eta) eta,
      loglik = function(y, mu, phi) dnorm(y, mu, phi[["sigma"]], log = TRUE),
      simulate = function(mu, phi) rnorm(length(mu), mu, phi[["sigma"]]),
      check_y = function(y) {
        if (!is.numeric(y)) stop("The gaussian family needs a numeric response.", call. = FALSE)
      },
      init_phi = function(y) log(max(sd(y), 1e-3) / 2)
    )
  )
}

resolve_family <- function(family) {
  if (is.character(family)) {
    family <- get(family, mode = "function", envir = asNamespace("stats"))
  }
  if (is.function(family)) family <- family()
  if (!inherits(family, "family")) {
    stop("'family' must be a family object such as gaussian().", call. = FALSE)
  }
  registry <- brtps_family_registry()
  spec <- registry[[family$family]]
  if (is.null(spec)) {
    stop("The '", family$family, "' family is not available yet. Available: ",
         paste(names(registry), collapse = ", "), ".", call. = FALSE)
  }
  if (!identical(family$link, spec$link)) {
    stop("The ", family$family, " family is only available with the '",
         spec$link, "' link.", call. = FALSE)
  }
  spec
}
