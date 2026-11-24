#' Priors for brtps models
#'
#' `brtps_priors()` collects the priors of a [brtps()] model. Anything you do
#' not set keeps its default, so you only list what you want to change.
#'
#' @param ... Named priors. The names available depend on the model and
#'   family:
#'   * `beta`: every regression coefficient (default [prior_flat()]).
#'   * `intercept`: the intercept only (default: same as `beta`).
#'   * `alpha`: the regTPS regularization parameter \eqn{\alpha} (default
#'     `prior_lognormal(0, 1)`, as in the paper).
#'   * `sigma`: the residual standard deviation of the gaussian family
#'     (default `prior_pc(u = 0.5, prob = 0.05)`, the PC prior of the paper).
#'
#' @details The defaults follow Cavieres and Krumscheid (2026). The prior on
#'   `sigma` assumes a response on a scale where a noise standard deviation
#'   above 0.5 is unlikely; for other scales change `u` (for example to the
#'   standard deviation of the response).
#'
#' @return A list of class `brtps_priors`.
#' @seealso [prior_pc()] and the other prior constructors.
#' @export
#' @examples
#' brtps_priors(sigma = prior_pc(u = 1, prob = 0.05),
#'              beta = prior_normal(0, 10))
brtps_priors <- function(...) {
  priors <- list(...)
  if (length(priors) && (is.null(names(priors)) || any(names(priors) == ""))) {
    stop("All priors in brtps_priors() must be named, e.g. sigma = prior_pc().",
         call. = FALSE)
  }
  bad <- !vapply(priors, inherits, logical(1), what = "brtps_prior")
  if (any(bad)) {
    stop("Not a prior: ", paste(names(priors)[bad], collapse = ", "),
         ". Use the prior_*() functions.", call. = FALSE)
  }
  structure(priors, class = "brtps_priors")
}

#' @export
print.brtps_priors <- function(x, ...) {
  if (!length(x)) {
    cat("brtps priors: all defaults\n")
  } else {
    for (nm in names(x)) cat(format(nm, width = 10), "~ ", format(x[[nm]]), "\n", sep = "")
  }
  invisible(x)
}

#' Prior distributions
#'
#' Constructors for the priors used in [brtps_priors()]. Priors for positive
#' parameters (`sigma`, `alpha`) are placed on the parameter itself; the
#' package adds the Jacobian needed because the samplers work on the log
#' scale.
#'
#' * `prior_flat()`: improper flat prior (on the log scale for positive
#'   parameters).
#' * `prior_normal()`, `prior_student_t()`: for regression coefficients.
#' * `prior_pc()`: penalized complexity prior for a standard deviation, set
#'   through \eqn{P(\sigma > u) = prob}. It is an exponential prior with rate
#'   \eqn{-\log(prob)/u} (Simpson et al., 2017).
#' * `prior_exponential()`, `prior_half_normal()`, `prior_half_cauchy()`,
#'   `prior_half_t()`, `prior_lognormal()`, `prior_gamma()`: for positive
#'   parameters.
#'
#' @param mean,sd Mean and standard deviation.
#' @param df,location,scale Degrees of freedom, location and scale.
#' @param u,prob Upper value and tail probability of the PC prior.
#' @param rate Rate parameter.
#' @param meanlog,sdlog Mean and standard deviation on the log scale.
#' @param shape Shape parameter.
#' @return An object of class `brtps_prior`.
#' @name priors
#' @examples
#' prior_pc(u = 0.5, prob = 0.05)
#' prior_lognormal(0, 1)
NULL

new_prior <- function(dist, code, par, support, label) {
  structure(list(dist = dist, code = as.integer(code), par = c(par, 0, 0, 0)[1:3],
                 support = support, label = label),
            class = "brtps_prior")
}

check_positive <- function(...) {
  args <- list(...)
  for (nm in names(args)) {
    x <- args[[nm]]
    if (!is.numeric(x) || length(x) != 1 || !is.finite(x) || x <= 0) {
      stop("'", nm, "' must be a single positive number.", call. = FALSE)
    }
  }
}

check_number <- function(...) {
  args <- list(...)
  for (nm in names(args)) {
    x <- args[[nm]]
    if (!is.numeric(x) || length(x) != 1 || !is.finite(x)) {
      stop("'", nm, "' must be a single finite number.", call. = FALSE)
    }
  }
}

#' @rdname priors
#' @export
prior_flat <- function() {
  new_prior("flat", 0L, numeric(0), "any", "flat (improper)")
}

#' @rdname priors
#' @export
prior_normal <- function(mean = 0, sd = 10) {
  check_number(mean = mean); check_positive(sd = sd)
  new_prior("normal", 1L, c(mean, sd), "real", sprintf("normal(%g, %g)", mean, sd))
}

#' @rdname priors
#' @export
prior_student_t <- function(df = 3, location = 0, scale = 2.5) {
  check_positive(df = df, scale = scale); check_number(location = location)
  new_prior("student_t", 2L, c(df, location, scale), "real",
            sprintf("student_t(%g, %g, %g)", df, location, scale))
}

#' @rdname priors
#' @export
prior_pc <- function(u = 0.5, prob = 0.05) {
  check_positive(u = u)
  if (!is.numeric(prob) || length(prob) != 1 || !(prob > 0 && prob < 1)) {
    stop("'prob' must be a probability strictly between 0 and 1.", call. = FALSE)
  }
  rate <- -log(prob) / u
  new_prior("pc", 1L, rate, "positive",
            sprintf("PC prior, P(> %g) = %g [exponential(rate = %.4g)]", u, prob, rate))
}

#' @rdname priors
#' @export
prior_exponential <- function(rate = 1) {
  check_positive(rate = rate)
  new_prior("exponential", 1L, rate, "positive", sprintf("exponential(%g)", rate))
}

#' @rdname priors
#' @export
prior_half_normal <- function(scale = 1) {
  check_positive(scale = scale)
  new_prior("half_normal", 2L, scale, "positive", sprintf("half-normal(0, %g)", scale))
}

#' @rdname priors
#' @export
prior_half_cauchy <- function(scale = 1) {
  check_positive(scale = scale)
  new_prior("half_cauchy", 3L, scale, "positive", sprintf("half-Cauchy(0, %g)", scale))
}

#' @rdname priors
#' @export
prior_half_t <- function(df = 3, scale = 1) {
  check_positive(df = df, scale = scale)
  new_prior("half_t", 4L, c(df, scale), "positive", sprintf("half-t(%g, 0, %g)", df, scale))
}

#' @rdname priors
#' @export
prior_lognormal <- function(meanlog = 0, sdlog = 1) {
  check_number(meanlog = meanlog); check_positive(sdlog = sdlog)
  new_prior("lognormal", 5L, c(meanlog, sdlog), "positive",
            sprintf("lognormal(%g, %g)", meanlog, sdlog))
}

#' @rdname priors
#' @export
prior_gamma <- function(shape = 2, rate = 0.1) {
  check_positive(shape = shape, rate = rate)
  new_prior("gamma", 6L, c(shape, rate), "positive", sprintf("gamma(%g, %g)", shape, rate))
}

#' @export
format.brtps_prior <- function(x, ...) x$label

#' @export
print.brtps_prior <- function(x, ...) {
  cat("<brtps prior>", format(x), "\n")
  invisible(x)
}

# Combines user priors with the defaults and turns them into TMB data.
resolve_priors <- function(priors, family, beta_names, model_priors) {
  if (is.null(priors)) priors <- brtps_priors()
  if (!inherits(priors, "brtps_priors")) {
    stop("'priors' must be created with brtps_priors().", call. = FALSE)
  }
  defaults <- c(list(beta = prior_flat()), model_priors, family$default_priors())
  valid <- c("beta", "intercept", setdiff(names(defaults), "beta"))
  unknown <- setdiff(names(priors), valid)
  if (length(unknown)) {
    stop("Unknown prior(s): ", paste(unknown, collapse = ", "),
         ". Valid names for this model: ", paste(valid, collapse = ", "), ".",
         call. = FALSE)
  }
  used <- defaults
  used[names(priors)] <- priors
  if (is.null(used$intercept)) used$intercept <- used$beta

  check_support <- function(p, name, support) {
    if (!p$support %in% c("any", support)) {
      stop("The prior for '", name, "' must be for a ",
           if (support == "real") "real-valued" else "positive", " parameter.",
           call. = FALSE)
    }
  }
  check_support(used$beta, "beta", "real")
  check_support(used$intercept, "intercept", "real")
  for (nm in setdiff(names(defaults), "beta")) check_support(used[[nm]], nm, "positive")

  beta_priors <- lapply(beta_names, function(nm) {
    if (nm == "(Intercept)") used$intercept else used$beta
  })
  phi_names <- family$phi_names
  list(
    used = used[c(if (length(beta_names)) c("beta", if ("(Intercept)" %in% beta_names) "intercept"),
                  names(model_priors), phi_names)],
    tmb = list(
      prior_beta_code = vapply(beta_priors, `[[`, integer(1), "code"),
      prior_beta_par = matrix(as.numeric(unlist(lapply(beta_priors, `[[`, "par"))),
                              ncol = 3, byrow = TRUE),
      prior_alpha_code = used$alpha$code,
      prior_alpha_par = used$alpha$par,
      prior_phi_code = vapply(used[phi_names], `[[`, integer(1), "code"),
      prior_phi_par = matrix(as.numeric(unlist(lapply(used[phi_names], `[[`, "par"))),
                             ncol = 3, byrow = TRUE)
    )
  )
}
