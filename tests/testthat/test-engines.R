# MCMC engines. These need optional packages and take a few seconds, so they
# are skipped on CRAN and when the packages are missing.

test_that("the SparseNUTS engine returns a valid fit", {
  skip_on_cran()
  skip_if_not_installed("SparseNUTS")
  fit <- suppressMessages(
    brtps(y ~ x1 + x2, data = brtps_sim, model = regtps_s(s1, s2, k_s = 40),
          chains = 2, num_samples = 300, cores = 1, seed = 11)
  )
  expect_s3_class(fit$fit, "tmbfit")
  expect_equal(posterior::nchains(fit$draws), 2)
  expect_equal(posterior::niterations(fit$draws), 300)
  expect_true("lp__" %in% posterior::variables(fit$draws))
  s <- summary(fit)
  expect_true(all(c("rhat", "ess_bulk") %in% names(s$hyper)))
  expect_equal(unname(coef(fit)[2]), 0.5, tolerance = 0.2)
})

test_that("the tmbstan engine returns a valid fit", {
  skip_on_cran()
  skip_if_not_installed("tmbstan")
  fit <- suppressWarnings(suppressMessages(
    brtps(y ~ x1 + x2, data = brtps_sim, model = regtps_s(s1, s2, k_s = 40),
          engine = "tmbstan", chains = 2, num_samples = 200, cores = 1, seed = 11,
          refresh = 0)
  ))
  expect_s4_class(fit$fit, "stanfit")
  expect_equal(posterior::ndraws(fit$draws), 400)
  expect_equal(unname(coef(fit)[2]), 0.5, tolerance = 0.2)
})

test_that("a missing SparseNUTS gives installation instructions", {
  skip_if(requireNamespace("SparseNUTS", quietly = TRUE))
  expect_error(brtps(y ~ x1, data = brtps_sim, model = regtps_s(s1, s2, k_s = 40)),
               "install_github")
})
