test_that("prior constructors store the right codes and parameters", {
  p <- prior_pc(u = 0.5, prob = 0.05)
  expect_s3_class(p, "brtps_prior")
  expect_equal(p$code, 1L)
  expect_equal(p$par[1], -log(0.05) / 0.5)
  expect_equal(prior_exponential(2)$par[1], 2)
  expect_equal(prior_lognormal(1, 2)$par[1:2], c(1, 2))
  expect_equal(prior_gamma(2, 3)$code, 6L)
  expect_equal(prior_flat()$code, 0L)
  expect_match(format(prior_normal(0, 10)), "normal\\(0, 10\\)")
})

test_that("prior constructors validate their arguments", {
  expect_error(prior_pc(u = -1), "positive")
  expect_error(prior_pc(prob = 1), "strictly between")
  expect_error(prior_normal(sd = 0), "positive")
  expect_error(prior_lognormal(meanlog = NA), "finite")
  expect_error(brtps_priors(prior_pc()), "named")
  expect_error(brtps_priors(sigma = 1), "Not a prior")
})

test_that("priors are resolved with defaults and support checks", {
  fam <- brtps:::resolve_family(gaussian())
  model_priors <- list(alpha = prior_lognormal(0, 1))
  r <- brtps:::resolve_priors(brtps_priors(), fam, c("(Intercept)", "x1"), model_priors)
  expect_equal(r$tmb$prior_beta_code, c(0L, 0L))
  expect_equal(r$tmb$prior_alpha_code, 5L)
  expect_equal(r$tmb$prior_phi_code, c(sigma = 1L))
  r2 <- brtps:::resolve_priors(brtps_priors(intercept = prior_normal(0, 5),
                                            beta = prior_normal(0, 1)),
                               fam, c("(Intercept)", "x1"), model_priors)
  expect_equal(r2$tmb$prior_beta_par[, 2], c(5, 1))
  expect_error(brtps:::resolve_priors(brtps_priors(sigma = prior_normal()), fam, "x1",
                                      model_priors), "positive parameter")
  expect_error(brtps:::resolve_priors(brtps_priors(beta = prior_gamma()), fam, "x1",
                                      model_priors), "real-valued")
  r0 <- brtps:::resolve_priors(brtps_priors(), fam, character(0), model_priors)
  expect_equal(dim(r0$tmb$prior_beta_par), c(0L, 3L))
})
