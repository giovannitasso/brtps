#' Summarise a brtps fit
#'
#' Posterior summaries of the regression coefficients and hyperparameters,
#' with convergence diagnostics (\eqn{\hat R}, bulk and tail effective sample
#' sizes) computed by \pkg{posterior}.
#'
#' @param object,x A `brtps` fit.
#' @param probs Quantiles to report.
#' @param digits Number of significant digits to print.
#' @param ... Not used.
#' @return An object of class `summary.brtps` with tables `fixed` and
#'   `hyper`.
#' @export
summary.brtps <- function(object, probs = c(0.025, 0.5, 0.975), ...) {
  beta_names <- colnames(object$data$X)
  hyper_names <- c(object$family$phi_names, "alpha")
  qfun <- function(x) stats::quantile(x, probs = probs, names = FALSE)
  qnames <- paste0("q", formatC(100 * probs, format = "fg", width = 1))
  summarise <- function(vars) {
    if (!length(vars)) return(NULL)
    d <- posterior::subset_draws(object$draws, variable = vars)
    tab <- as.data.frame(posterior::summarise_draws(d, "mean", "sd", qfun))
    names(tab) <- c("variable", "mean", "sd", qnames)
    if (object$engine != "laplace") {
      conv <- as.data.frame(posterior::summarise_draws(d, "rhat", "ess_bulk", "ess_tail"))
      tab <- cbind(tab, conv[, -1, drop = FALSE])
    }
    rownames(tab) <- tab$variable
    tab$variable <- NULL
    tab
  }
  structure(
    list(
      fixed = summarise(beta_names),
      hyper = summarise(hyper_names),
      basis = object$model_info$basis,
      family = object$family,
      engine = object$engine,
      engine_label = object$engine_label,
      diagnostics = object$diagnostics,
      nobs = length(object$data$y),
      ndraws = posterior::ndraws(object$draws),
      nchains = posterior::nchains(object$draws),
      priors = object$priors,
      formula = object$formula,
      elapsed = object$elapsed
    ),
    class = "summary.brtps"
  )
}

#' @rdname summary.brtps
#' @export
print.summary.brtps <- function(x, digits = 3, ...) {
  b <- x$basis
  cat("Bayesian regularized thin plate spline model (brtps)\n")
  cat(" Formula: ", deparse1(x$formula), "\n", sep = "")
  cat(" Family:  ", x$family$family, " (", x$family$link, " link)\n", sep = "")
  cat(" Spatial: regTPS, k = ", b$k, ", M = ", b$M, " (",
      round(100 * b$cumvar[b$M], 1), "% of the prior variance), null space ",
      b$d_null, if (b$center) " (centred)" else "", "\n", sep = "")
  cat(" Data:    ", x$nobs, " observations\n", sep = "")
  cat(" Engine:  ", x$engine_label, "; ", x$ndraws, " draws",
      if (x$engine != "laplace") paste0(" from ", x$nchains, " chain(s)"), "; ",
      format_time(x$elapsed), "\n", sep = "")
  if (!is.null(x$fixed)) {
    cat("\nRegression coefficients:\n")
    print(signif(as.matrix(x$fixed), digits))
  }
  cat("\nHyperparameters:\n")
  print(signif(as.matrix(x$hyper), digits))
  cat("\nPriors:\n")
  for (nm in names(x$priors)) cat("  ", format(nm, width = 10), "~ ", format(x$priors[[nm]]), "\n", sep = "")
  if (x$engine != "laplace") {
    all <- rbind(x$fixed, x$hyper)
    warn <- character(0)
    if (any(all$rhat > 1.01, na.rm = TRUE)) warn <- c(warn, "some Rhat > 1.01")
    if (any(all$ess_bulk < 400, na.rm = TRUE)) warn <- c(warn, "some bulk ESS < 400")
    if (isTRUE(x$diagnostics$divergences > 0)) {
      warn <- c(warn, paste(x$diagnostics$divergences, "divergent transition(s)"))
    }
    if (length(warn)) {
      cat("\nWarning: ", paste(warn, collapse = "; "),
          ". Consider more draws or a higher adapt_delta.\n", sep = "")
    }
  } else {
    cat("\nNote: draws come from a Gaussian approximation at the mode, not MCMC.\n")
  }
  invisible(x)
}

#' @rdname summary.brtps
#' @export
print.brtps <- function(x, digits = 3, ...) {
  print(summary(x), digits = digits, ...)
  invisible(x)
}

format_time <- function(sec) {
  if (sec < 60) return(sprintf("%.1f seconds", sec))
  if (sec < 3600) return(sprintf("%.1f minutes", sec / 60))
  sprintf("%.1f hours", sec / 3600)
}

#' Posterior means of the regression coefficients
#'
#' @param object A `brtps` fit.
#' @param ... Not used.
#' @return A named numeric vector.
#' @export
coef.brtps <- function(object, ...) {
  vars <- colnames(object$data$X)
  if (!length(vars)) return(numeric(0))
  colMeans(draws_matrix_plain(posterior::subset_draws(object$draws, variable = vars)))
}

#' @export
nobs.brtps <- function(object, ...) length(object$data$y)

#' Predictions from a brtps fit
#'
#' Posterior predictions at the observed locations or at new ones.
#'
#' @param object A `brtps` fit.
#' @param newdata Optional data frame with the covariates and the coordinate
#'   columns used in [regtps_s()]. `NULL` uses the data of the fit.
#' @param type What to predict:
#'   * `"response"`: the mean of the response, \eqn{g^{-1}(\eta)}.
#'   * `"link"`: the linear predictor \eqn{\eta}.
#'   * `"spatial"`: the spatial field \eqn{u(s)} alone (no covariates needed
#'     in `newdata`).
#'   * `"predictive"`: draws of new observations (posterior predictive).
#' @param summary If `TRUE`, return posterior summaries; if `FALSE`, a
#'   matrix of draws (one row per draw, one column per location).
#' @param probs Quantiles for the summaries.
#' @param ndraws Optional number of posterior draws to use (randomly chosen).
#' @param seed Optional seed for `ndraws` and for `type = "predictive"`.
#' @param ... Not used.
#' @return A data frame with columns `mean`, `sd` and the quantiles, or a
#'   matrix of draws.
#' @export
#' @examples
#' fit <- brtps(y ~ x1 + x2, data = brtps_sim,
#'              model = regtps_s(s1, s2, k_s = 60), engine = "laplace")
#' grid <- expand.grid(s1 = seq(0, 1, length.out = 20), s2 = seq(0, 1, length.out = 20))
#' field <- predict(fit, newdata = grid, type = "spatial")
#' head(field)
predict.brtps <- function(object, newdata = NULL,
                          type = c("response", "link", "spatial", "predictive"),
                          summary = TRUE, probs = c(0.025, 0.5, 0.975),
                          ndraws = NULL, seed = NULL, ...) {
  type <- match.arg(type)
  if (!is.null(seed)) set.seed(seed)
  d <- draws_matrix_plain(object$draws)
  if (!is.null(ndraws)) {
    check_count(ndraws, "ndraws")
    d <- d[sample.int(nrow(d), min(ndraws, nrow(d))), , drop = FALSE]
  }
  info <- object$model_info
  Phi <- field_design(object$model, info, newdata)
  v <- info$basis$v
  zt <- d[, sprintf("z_tilde[%d]", seq_along(v)), drop = FALSE]
  z <- zt / sqrt(1 + outer(as.numeric(d[, "alpha"]), v))
  eta <- z %*% t(Phi)                                  # draws x locations

  if (type != "spatial") {
    beta_names <- colnames(object$data$X)
    if (is.null(newdata)) {
      X <- object$data$X
      offset <- object$data$offset
    } else {
      tt <- delete.response(object$terms)
      mf <- model.frame(tt, newdata, xlev = object$xlevels, na.action = na.pass)
      X <- model.matrix(tt, mf, contrasts.arg = object$contrasts)
      offset <- model.offset(mf)
      if (is.null(offset)) offset <- rep(0, nrow(X))
    }
    if (length(beta_names)) {
      eta <- eta + d[, beta_names, drop = FALSE] %*% t(X[, beta_names, drop = FALSE])
    }
    eta <- sweep(eta, 2, offset, "+")
  }

  out <- switch(type,
    link = eta,
    spatial = eta,
    response = object$family$linkinv(eta),
    predictive = {
      mu <- object$family$linkinv(eta)
      phi <- d[, object$family$phi_names, drop = FALSE]
      t(vapply(seq_len(nrow(mu)), function(i) {
        object$family$simulate(mu[i, ], as.list(phi[i, ]))
      }, numeric(ncol(mu))))
    }
  )
  out <- matrix(out, nrow = nrow(d))
  if (!summary) return(out)
  summarise_columns(out, probs)
}

summarise_columns <- function(x, probs) {
  q <- t(apply(x, 2, stats::quantile, probs = probs, names = FALSE))
  if (length(probs) == 1) q <- t(q)
  res <- data.frame(mean = colMeans(x), sd = apply(x, 2, sd), q)
  names(res)[-(1:2)] <- paste0("q", formatC(100 * probs, format = "fg", width = 1))
  rownames(res) <- NULL
  res
}

#' @rdname predict.brtps
#' @export
fitted.brtps <- function(object, ...) {
  predict(object, type = "response")$mean
}

#' @rdname predict.brtps
#' @export
residuals.brtps <- function(object, ...) {
  object$data$y - fitted(object)
}

#' Pointwise log-likelihood
#'
#' Log-likelihood of each observation for each posterior draw, as needed by
#' \pkg{loo} and for comparing models.
#'
#' @param object A `brtps` fit.
#' @param array If `TRUE`, return an iterations x chains x observations
#'   array; otherwise a draws x observations matrix.
#' @return A matrix or array of log-likelihood values.
#' @export
brtps_loglik <- function(object, array = FALSE) {
  eta <- predict(object, type = "link", summary = FALSE)
  mu <- object$family$linkinv(eta)
  d <- draws_matrix_plain(object$draws)
  phi <- d[, object$family$phi_names, drop = FALSE]
  y <- object$data$y
  ll <- t(vapply(seq_len(nrow(mu)), function(i) {
    object$family$loglik(y, mu[i, ], as.list(phi[i, ]))
  }, numeric(length(y))))
  if (!array) return(ll)
  niter <- posterior::niterations(object$draws)
  nchain <- posterior::nchains(object$draws)
  # draws are stacked chain after chain, iterations varying fastest
  array(ll, dim = c(niter, nchain, length(y)))
}

# Draws as a plain numeric matrix (draws x variables), chains stacked.
draws_matrix_plain <- function(draws) {
  d <- posterior::as_draws_matrix(draws)
  matrix(as.numeric(d), nrow = nrow(d), dimnames = list(NULL, colnames(d)))
}

#' Leave-one-out cross-validation
#'
#' PSIS-LOO for a `brtps` fit, computed with [loo::loo()].
#'
#' @param x A `brtps` fit.
#' @param ... Passed to [loo::loo()].
#' @param cores Number of cores.
#' @return A `loo` object.
#' @exportS3Method loo::loo
loo.brtps <- function(x, ..., cores = getOption("mc.cores", 1L)) {
  if (!requireNamespace("loo", quietly = TRUE)) {
    stop("Please install the loo package.", call. = FALSE)
  }
  ll <- brtps_loglik(x, array = TRUE)
  r_eff <- if (x$engine == "laplace") rep(1, dim(ll)[3]) else
    loo::relative_eff(exp(ll), cores = cores)
  loo::loo(ll, r_eff = r_eff, cores = cores, ...)
}

#' Posterior draws of a brtps fit
#'
#' Methods for the \pkg{posterior} generics, so a fit can be passed to
#' \pkg{posterior} and \pkg{bayesplot} functions. The draws contain the
#' regression coefficients (named as in the design matrix), the dispersion
#' parameter (`sigma`), `alpha`, the non-centred coefficients
#' `z_tilde[k]` and, for MCMC engines, `lp__`.
#'
#' @param x A `brtps` fit.
#' @param variable Optional character vector of variables to keep.
#' @param ... Not used.
#' @return A `draws_array` or `draws_df` object.
#' @name as_draws.brtps
#' @examples
#' fit <- brtps(y ~ x1, data = brtps_sim, model = regtps_s(s1, s2, k_s = 40),
#'              engine = "laplace")
#' posterior::as_draws_df(fit, variable = c("x1", "sigma", "alpha"))
NULL

#' @rdname as_draws.brtps
#' @exportS3Method posterior::as_draws_array
as_draws_array.brtps <- function(x, variable = NULL, ...) {
  d <- x$draws
  if (!is.null(variable)) d <- posterior::subset_draws(d, variable = variable)
  d
}

#' @rdname as_draws.brtps
#' @exportS3Method posterior::as_draws
as_draws.brtps <- function(x, variable = NULL, ...) {
  as_draws_array.brtps(x, variable = variable)
}

#' @rdname as_draws.brtps
#' @exportS3Method posterior::as_draws_df
as_draws_df.brtps <- function(x, variable = NULL, ...) {
  posterior::as_draws_df(as_draws_array.brtps(x, variable = variable))
}
