# The TMB template must equal the negative log posterior written out in R.

make_obj <- function(priors) {
  mod <- regtps_s(s1, s2, k_s = 40)
  built <- brtps:::build_model(mod, brtps_sim, TRUE)
  built <- brtps:::finish_model(mod, built, rep(TRUE, nrow(brtps_sim)))
  fam <- brtps:::resolve_family(gaussian())
  X <- model.matrix(~ x1 + x2, brtps_sim)
  pri <- brtps:::resolve_priors(priors, fam, colnames(X),
                                brtps:::model_default_priors(mod))
  n <- nrow(brtps_sim)
  dat <- c(list(model = "regtps_s", y = brtps_sim$y, X = X, offset = rep(0, n),
                size = rep(1, n)),
           built$tmb_data, list(family = 0L), pri$tmb)
  init <- brtps:::initial_values(brtps_sim$y, X, rep(0, n), fam, built)
  list(obj = TMB::MakeADFun(dat, init, DLL = "brtps", silent = TRUE), dat = dat)
}

nll_r <- function(par, dat, log_prior_sigma, log_prior_logalpha, log_prior_beta) {
  p <- ncol(dat$X)
  M <- length(dat$v)
  beta <- par[1:p]
  log_alpha <- par[p + 1]
  log_sigma <- par[p + 2]
  z_tilde <- par[p + 2 + seq_len(M)]
  z <- z_tilde / sqrt(1 + exp(log_alpha) * dat$v)
  eta <- drop(dat$X %*% beta + dat$Phi %*% z)
  sigma <- exp(log_sigma)
  -(sum(dnorm(dat$y, eta, sigma, log = TRUE)) + sum(dnorm(z_tilde, log = TRUE)) +
      log_prior_sigma(sigma) + log_sigma +            # Jacobian of log(sigma)
      log_prior_logalpha(log_alpha) + sum(log_prior_beta(beta)))
}

expect_same_density <- function(priors, lp_sigma, lp_logalpha, lp_beta) {
  o <- make_obj(priors)
  set.seed(42)
  for (r in 1:5) {
    par <- o$obj$par + rnorm(length(o$obj$par), sd = 0.5)
    expect_equal(unname(o$obj$fn(par)), unname(nll_r(par, o$dat, lp_sigma, lp_logalpha, lp_beta)),
                 tolerance = 1e-10)
  }
}

test_that("default priors (paper) match", {
  rate <- -log(0.05) / 0.5
  expect_same_density(brtps_priors(),
                      function(s) dexp(s, rate, log = TRUE),
                      function(la) dnorm(la, 0, 1, log = TRUE),
                      function(b) 0 * b)
})

test_that("half-normal, gamma and normal priors match", {
  expect_same_density(
    brtps_priors(sigma = prior_half_normal(2), alpha = prior_gamma(2, 0.5),
                 beta = prior_normal(0, 3)),
    function(s) log(2) + dnorm(s, 0, 2, log = TRUE),
    function(la) dgamma(exp(la), 2, rate = 0.5, log = TRUE) + la,
    function(b) dnorm(b, 0, 3, log = TRUE))
})

test_that("half-Cauchy, half-t, exponential and Student-t priors match", {
  expect_same_density(
    brtps_priors(sigma = prior_half_cauchy(1), beta = prior_student_t(3, 0, 2.5)),
    function(s) log(2) + dcauchy(s, 0, 1, log = TRUE),
    function(la) dnorm(la, 0, 1, log = TRUE),
    function(b) dt(b / 2.5, 3, log = TRUE) - log(2.5))
  expect_same_density(
    brtps_priors(sigma = prior_half_t(4, 1.5), alpha = prior_exponential(0.2),
                 intercept = prior_normal(1, 5)),
    function(s) log(2) + dt(s / 1.5, 4, log = TRUE) - log(1.5),
    function(la) dexp(exp(la), 0.2, log = TRUE) + la,
    function(b) c(dnorm(b[1], 1, 5, log = TRUE), 0, 0))
})

test_that("the gradient from TMB matches finite differences", {
  o <- make_obj(brtps_priors())
  par <- o$obj$par + 0.1
  num <- vapply(seq_along(par), function(j) {
    h <- 1e-6
    e <- replace(numeric(length(par)), j, h)
    (o$obj$fn(par + e) - o$obj$fn(par - e)) / (2 * h)
  }, numeric(1))
  expect_equal(as.numeric(o$obj$gr(par)), num, tolerance = 1e-5)
})
