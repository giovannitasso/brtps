#' Fit a Bayesian regularized thin plate spline model
#'
#' Fits a spatial regression model whose spatial random field is
#' approximated with the regTPS construction of Cavieres and Krumscheid
#' (2026). For the gaussian family the model is
#' \deqn{y_i = x_i^\top \beta + u(s_i) + \varepsilon_i, \qquad
#'   \varepsilon_i \sim N(0, \sigma^2),}
#' where \eqn{u(s) = \sum_{k=1}^M z_k \hat\phi_k(s)} and
#' \eqn{z_k = \tilde z_k / \sqrt{1 + \alpha v_k}} with
#' \eqn{\tilde z_k \sim N(0, 1)} (non-centred parameterization, eq. A20).
#'
#' @param formula Formula for the fixed effects, for example `y ~ x1 + x2`.
#'   `offset()` terms are allowed. Use `y ~ 1` for an intercept-only mean, or
#'   `y ~ 0` to reproduce the models of the paper, where the spline null space
#'   carries the mean.
#' @param data Data frame with the response, covariates and coordinates.
#'   Rows with missing values in any of them are dropped.
#' @param model Spatial component, created with [regtps_s()].
#' @param family Observation model. Only [gaussian()] (identity link) is
#'   available in this version.
#' @param priors Priors, created with [brtps_priors()]. Defaults follow the
#'   paper.
#' @param engine Inference method:
#'   * `"snuts"` (default): sparse NUTS from \pkg{SparseNUTS}.
#'   * `"tmbstan"`: NUTS from \pkg{tmbstan} (returns a `stanfit` inside).
#'   * `"laplace"`: Gaussian approximation at the posterior mode, for quick
#'     checks; no MCMC.
#' @param chains Number of MCMC chains.
#' @param num_samples Draws per chain after warmup. For `"laplace"`, the
#'   total number of draws from the approximation.
#' @param num_warmup Warmup iterations per chain. `NULL` uses the engine
#'   default (\pkg{SparseNUTS} chooses it; \pkg{tmbstan} uses `num_samples`).
#' @param cores Number of cores used to run chains in parallel.
#' @param seed Optional random seed.
#' @param control Optional list of NUTS settings such as `adapt_delta` and
#'   `max_treedepth`.
#' @param silent Silence the output of \pkg{TMB}.
#' @param ... Further arguments passed to `SparseNUTS::sample_snuts()` or
#'   `tmbstan::tmbstan()`.
#'
#' @return An object of class `brtps`. Use [summary()][summary.brtps],
#'   [predict()][predict.brtps], [posterior::as_draws_array()] (to work with
#'   \pkg{posterior} and \pkg{bayesplot}) and [loo()][loo.brtps].
#'
#' @references Cavieres, J. and Krumscheid, S. (2026). regTPS: a novel
#'   approach to approximate a Gaussian random field for Bayesian spatial
#'   modeling. *Statistics and Computing*, 36, 261.
#'   \doi{10.1007/s11222-026-10990-x}
#' @export
#' @examples
#' # Quick fit with the Laplace approximation (no MCMC)
#' fit <- brtps(y ~ x1 + x2, data = brtps_sim,
#'              model = regtps_s(s1, s2, k_s = 60), engine = "laplace")
#' fit
#'
#' \donttest{
#' # Full Bayesian fit with sparse NUTS (needs the SparseNUTS package)
#' if (requireNamespace("SparseNUTS", quietly = TRUE)) {
#'   fit_mcmc <- brtps(y ~ x1 + x2, data = brtps_sim,
#'                     model = regtps_s(s1, s2, k_s = 60),
#'                     chains = 2, num_samples = 500)
#'   summary(fit_mcmc)
#' }
#' }
brtps <- function(formula, data, model, family = gaussian(), priors = brtps_priors(),
                  engine = c("snuts", "tmbstan", "laplace"), chains = 4,
                  num_samples = 1000, num_warmup = NULL,
                  cores = getOption("mc.cores", 1L), seed = NULL, control = NULL,
                  silent = TRUE, ...) {
  call <- match.call()
  engine <- match.arg(engine)
  if (missing(data) || !is.data.frame(data)) {
    stop("'data' must be a data frame.", call. = FALSE)
  }
  if (missing(model) || !inherits(model, "brtps_model")) {
    stop("'model' must be a spatial component created with regtps_s().", call. = FALSE)
  }
  check_count(chains, "chains")
  check_count(num_samples, "num_samples")
  if (!is.null(num_warmup)) check_count(num_warmup, "num_warmup", zero = TRUE)
  fam <- resolve_family(family)

  # Response and fixed effects ------------------------------------------------
  formula <- stats::as.formula(formula)
  tt <- terms(formula, data = data)
  if (attr(tt, "response") == 0) {
    stop("The formula needs a response, for example y ~ x1.", call. = FALSE)
  }
  has_intercept <- attr(tt, "intercept") == 1

  built <- build_model(model, data, has_intercept)
  mf_all <- model.frame(tt, data, na.action = na.pass)
  rows <- complete.cases(mf_all) & complete.cases(built$coords)
  if (!any(rows)) stop("No row of 'data' is complete.", call. = FALSE)
  if (any(!rows)) message("Dropping ", sum(!rows), " row(s) with missing values.")
  mf <- model.frame(tt, data[rows, , drop = FALSE])
  y <- model.response(mf)
  fam$check_y(y)
  y <- as.numeric(y)
  X <- model.matrix(attr(mf, "terms"), mf)
  if (ncol(X) > 0 && qr(X)$rank < ncol(X)) {
    stop("The fixed-effects design matrix is rank deficient; check the ",
         "formula for collinear terms.", call. = FALSE)
  }
  offset <- model.offset(mf)
  offset <- if (is.null(offset)) rep(0, length(y)) else as.numeric(offset)

  # Spatial basis and priors -----------------------------------------------------
  built <- finish_model(model, built, rows)
  pri <- resolve_priors(priors, fam, colnames(X), model_default_priors(model))
  check_sigma_scale(priors, fam, y)

  # TMB object -----------------------------------------------------------------
  tmb_data <- c(
    list(model = built$tmb_model, y = y, X = X, offset = offset,
         size = rep(1, length(y))),
    built$tmb_data,
    list(family = fam$code),
    pri$tmb
  )
  init <- initial_values(y, X, offset, fam, built)
  obj <- TMB::MakeADFun(data = tmb_data, parameters = init, random = "z_tilde",
                        DLL = "brtps", silent = silent)

  # Inference --------------------------------------------------------------------
  start <- proc.time()[["elapsed"]]
  res <- switch(engine,
    snuts = engine_snuts(obj, chains = chains, num_samples = num_samples,
                         num_warmup = num_warmup, cores = cores, seed = seed,
                         control = control, ...),
    tmbstan = engine_tmbstan(obj, chains = chains, num_samples = num_samples,
                             num_warmup = num_warmup, cores = cores, seed = seed,
                             control = control, ...),
    laplace = engine_laplace(obj, num_samples = num_samples, seed = seed)
  )
  elapsed <- proc.time()[["elapsed"]] - start

  draws <- tmb_to_draws(res$array, res$lp, obj, colnames(X), fam$phi_names)

  structure(
    list(
      call = call,
      formula = formula,
      family = fam,
      model = model,
      model_info = built,
      priors = pri$used,
      terms = attr(mf, "terms"),
      xlevels = .getXlevels(attr(mf, "terms"), mf),
      contrasts = attr(X, "contrasts"),
      data = list(y = y, X = X, offset = offset, rows = which(rows)),
      tmb = list(obj = obj, data = tmb_data, init = init),
      engine = engine,
      engine_label = res$engine_label,
      fit = res$fit,
      diagnostics = res$diagnostics,
      draws = draws,
      elapsed = elapsed,
      version = utils::packageVersion("brtps")
    ),
    class = "brtps"
  )
}

initial_values <- function(y, X, offset, family, built) {
  M <- length(built$tmb_data$v)
  r <- y - offset
  beta <- if (ncol(X) > 0) stats::lm.fit(X, r)$coefficients else numeric(0)
  res <- if (ncol(X) > 0) r - drop(X %*% beta) else r
  list(
    beta = unname(beta),
    log_alpha = built$basis$log_alpha_init,
    log_phi = family$init_phi(res),
    z_tilde = rep(0, M)
  )
}

# Turns engine output (draws on the TMB parameter scale, as an
# iteration x chain x parameter array in the order of obj$env$par) into a
# posterior draws_array on the user scale.
tmb_to_draws <- function(arr, lp, obj, beta_names, phi_names) {
  par_names <- names(obj$env$par)
  if (dim(arr)[3] != length(par_names)) {
    stop("Internal error: the sampler returned ", dim(arr)[3],
         " parameters but the model has ", length(par_names), ".", call. = FALSE)
  }
  idx_beta <- which(par_names == "beta")
  idx_phi <- which(par_names == "log_phi")
  idx_alpha <- which(par_names == "log_alpha")
  idx_z <- which(par_names == "z_tilde")
  pieces <- list(
    arr[, , idx_beta, drop = FALSE],
    exp(arr[, , idx_phi, drop = FALSE]),
    exp(arr[, , idx_alpha, drop = FALSE]),
    arr[, , idx_z, drop = FALSE]
  )
  vars <- c(beta_names, phi_names, "alpha", sprintf("z_tilde[%d]", seq_along(idx_z)))
  if (!is.null(lp)) {
    pieces <- c(pieces, list(array(lp, dim = c(dim(arr)[1:2], 1))))
    vars <- c(vars, "lp__")
  }
  out <- array(NA_real_, dim = c(dim(arr)[1:2], length(vars)))
  pos <- 0
  for (p in pieces) {
    k <- dim(p)[3]
    if (k > 0) out[, , pos + seq_len(k)] <- p
    pos <- pos + k
  }
  dimnames(out) <- list(NULL, NULL, vars)
  posterior::as_draws_array(out)
}

check_count <- function(x, name, zero = FALSE) {
  ok <- is.numeric(x) && length(x) == 1 && !is.na(x) && x == round(x) &&
    (x > 0 || (zero && x == 0))
  if (!ok) stop("'", name, "' must be a whole number", if (!zero) " of at least 1",
                ".", call. = FALSE)
}

# The default prior on sigma from the paper assumes a unit-scale response.
check_sigma_scale <- function(priors, family, y) {
  if (!"sigma" %in% family$phi_names || "sigma" %in% names(priors)) return(invisible())
  if (sd(y) > 2.5) {
    message("Note: the default prior for sigma, prior_pc(u = 0.5, prob = 0.05), ",
            "assumes a response on a unit scale but sd(y) = ", signif(sd(y), 3),
            ". Consider brtps_priors(sigma = prior_pc(u = ", signif(sd(y), 2),
            ", prob = 0.05)) or rescaling y.")
  }
  invisible()
}
