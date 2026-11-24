# Inference engines.
#
# Every engine receives the TMB object (built with z_tilde as random effects)
# and returns a list with
#   array:        iterations x chains x parameters, TMB scale, order of obj$env$par
#   lp:           iterations x chains matrix of log-posterior values, or NULL
#   fit:          the raw output of the engine
#   diagnostics:  a list with engine-specific diagnostics
#   engine_label: a short description for printing
# Adding an engine means writing one more function with this contract.

engine_snuts <- function(obj, chains, num_samples, num_warmup, cores, seed,
                         control, ...) {
  if (!requireNamespace("SparseNUTS", quietly = TRUE)) {
    stop("engine = \"snuts\" needs the SparseNUTS package. Install it with\n",
         "  install.packages(\"StanEstimators\", repos = c(\"https://andrjohns.r-universe.dev\", ",
         "\"https://cloud.r-project.org\"))\n",
         "  remotes::install_github(\"noaa-afsc/SparseNUTS\")\n",
         "or use engine = \"laplace\" for a quick approximation.", call. = FALSE)
  }
  if (cores > 1 && !(requireNamespace("future", quietly = TRUE) &&
                     requireNamespace("callr", quietly = TRUE))) {
    message("Running the chains one after another: parallel chains with ",
            "SparseNUTS need the 'future' and 'callr' packages.")
    cores <- 1
  }
  fit <- SparseNUTS::sample_snuts(obj, num_samples = num_samples,
                                  num_warmup = num_warmup, chains = chains,
                                  cores = cores, seed = seed, control = control,
                                  laplace = FALSE, print = FALSE, ...)
  post <- fit$samples[-seq_len(fit$warmup), , , drop = FALSE]
  np <- dim(post)[3]
  sp <- SparseNUTS::extract_sampler_params(fit)
  list(
    array = post[, , -np, drop = FALSE],
    lp = matrix(post[, , np], nrow = dim(post)[1]),
    fit = fit,
    diagnostics = list(
      divergences = sum(sp[["divergent__"]]),
      max_treedepth_hits = sum(sp[["treedepth__"]] >= fit$max_treedepth),
      metric = fit$metric
    ),
    engine_label = sprintf("SparseNUTS (%s metric)", fit$metric)
  )
}

engine_tmbstan <- function(obj, chains, num_samples, num_warmup, cores, seed,
                           control, ...) {
  if (!requireNamespace("tmbstan", quietly = TRUE)) {
    stop("engine = \"tmbstan\" needs the tmbstan package: install.packages(\"tmbstan\").",
         call. = FALSE)
  }
  # Start every chain at the posterior mode, as in the paper.
  stats::nlminb(obj$par, obj$fn, obj$gr)
  if (is.null(num_warmup)) num_warmup <- num_samples
  args <- list(obj = obj, chains = chains, iter = num_warmup + num_samples,
               warmup = num_warmup, cores = cores, init = "last.par.best",
               laplace = FALSE, ...)
  if (!is.null(control)) args$control <- control
  if (!is.null(seed)) args$seed <- seed
  fit <- do.call(tmbstan::tmbstan, args)
  post <- as.array(fit)
  np <- dim(post)[3]
  sp <- rstan::get_sampler_params(fit, inc_warmup = FALSE)
  max_depth <- if (is.null(control$max_treedepth)) 10 else control$max_treedepth
  list(
    array = post[, , -np, drop = FALSE],
    lp = matrix(post[, , np], nrow = dim(post)[1]),
    fit = fit,
    diagnostics = list(
      divergences = sum(vapply(sp, function(s) sum(s[, "divergent__"]), numeric(1))),
      max_treedepth_hits = sum(vapply(sp, function(s) sum(s[, "treedepth__"] >= max_depth),
                                      numeric(1)))
    ),
    engine_label = "tmbstan (NUTS)"
  )
}

engine_laplace <- function(obj, num_samples, seed) {
  opt <- stats::nlminb(obj$par, obj$fn, obj$gr,
                       control = list(eval.max = 10000, iter.max = 10000))
  if (opt$convergence != 0) {
    warning("The optimizer did not report convergence: ", opt$message, call. = FALSE)
  }
  sdr <- TMB::sdreport(obj, getJointPrecision = TRUE)
  mode <- obj$env$last.par.best
  Q <- as.matrix(sdr$jointPrecision)
  R <- tryCatch(chol(Q), error = function(e) NULL)
  if (is.null(R)) {
    stop("The Hessian at the mode is not positive definite, so the Laplace ",
         "approximation cannot be sampled. Try an MCMC engine.", call. = FALSE)
  }
  if (!is.null(seed)) set.seed(seed)
  Z <- matrix(rnorm(length(mode) * num_samples), nrow = length(mode))
  draws <- mode + backsolve(R, Z)               # columns are N(mode, Q^-1) draws
  list(
    array = array(t(draws), dim = c(num_samples, 1L, length(mode))),
    lp = NULL,
    fit = list(opt = opt, sdreport = sdr, mode = mode),
    diagnostics = list(convergence = opt$convergence,
                       max_gradient = max(abs(obj$gr(opt$par)))),
    engine_label = "Laplace approximation (Gaussian at the mode)"
  )
}
