fit_lap <- function(formula = y ~ x1 + x2, data = brtps_sim, ...) {
  brtps(formula, data = data, model = regtps_s(s1, s2, k_s = 60),
        engine = "laplace", num_samples = 400, seed = 1, ...)
}

test_that("the Laplace engine recovers the simulated parameters", {
  fit <- fit_lap()
  expect_s3_class(fit, "brtps")
  b <- coef(fit)
  expect_named(b, c("(Intercept)", "x1", "x2"))
  # with a centred field the intercept estimates 1 + mean(u)
  expect_equal(unname(b[1]), 1 + mean(brtps_sim$u), tolerance = 0.1)
  expect_equal(unname(b[2]), 0.5, tolerance = 0.2)
  expect_equal(unname(b[3]), -0.8, tolerance = 0.25)
  s <- summary(fit)
  expect_true(s$hyper["sigma", "mean"] > 0.2 && s$hyper["sigma", "mean"] < 0.4)
  field <- predict(fit, type = "spatial")
  expect_gt(cor(field$mean, brtps_sim$u), 0.9)
})

test_that("draws have the documented variables", {
  fit <- fit_lap()
  vars <- posterior::variables(fit$draws)
  expect_true(all(c("(Intercept)", "x1", "x2", "sigma", "alpha", "z_tilde[1]") %in% vars))
  expect_equal(posterior::ndraws(fit$draws), 400)
  d <- posterior::as_draws_df(fit, variable = c("x1", "alpha"))
  expect_equal(posterior::variables(d), c("x1", "alpha"))
})

test_that("predictions are consistent", {
  fit <- fit_lap()
  p_all <- predict(fit)
  expect_equal(dim(p_all), c(nrow(brtps_sim), 5L))
  expect_named(p_all, c("mean", "sd", "q2.5", "q50", "q97.5"))
  p_new <- predict(fit, newdata = brtps_sim[1:7, ])
  expect_equal(p_new$mean, p_all$mean[1:7], tolerance = 1e-8)
  expect_equal(fitted(fit), p_all$mean)
  expect_equal(residuals(fit), brtps_sim$y - p_all$mean)
  grid <- expand.grid(s1 = c(0.2, 0.5), s2 = c(0.3, 0.7))
  expect_equal(nrow(predict(fit, newdata = grid, type = "spatial")), 4)
  expect_error(predict(fit, newdata = data.frame(a = 1)), "coordinate column")
  draws <- predict(fit, type = "predictive", summary = FALSE, ndraws = 25, seed = 3)
  expect_equal(dim(draws), c(25L, nrow(brtps_sim)))
  ll <- brtps_loglik(fit)
  expect_equal(dim(ll), c(400L, nrow(brtps_sim)))
  expect_equal(dim(brtps_loglik(fit, array = TRUE)), c(400L, 1L, nrow(brtps_sim)))
})

test_that("intercept-free formulas use the paper's null space", {
  fit <- fit_lap(y ~ 0 + x1)
  expect_equal(fit$model_info$basis$d_null, 3L)
  expect_false(fit$model_info$basis$center)
  fit_c <- fit_lap(y ~ x1)
  expect_true(fit_c$model_info$basis$center)
})

test_that("offsets and missing values are handled", {
  dat <- brtps_sim
  dat$off <- 0.25
  fit <- fit_lap()
  fit_off <- fit_lap(y ~ x1 + x2 + offset(off), data = dat)
  expect_equal(unname(coef(fit)[1] - coef(fit_off)[1]), 0.25, tolerance = 1e-3)
  dat$x1[3] <- NA
  dat$s2[8] <- NA
  expect_message(fit_na <- fit_lap(data = dat), "Dropping 2 row")
  expect_equal(length(fit_na$data$y), nrow(dat) - 2)
})

test_that("coordinates can be given as strings or vectors", {
  fit_chr <- brtps(y ~ x1, data = brtps_sim, model = regtps_s("s1", "s2", k_s = 40),
                   engine = "laplace", num_samples = 50, seed = 1)
  fit_sym <- brtps(y ~ x1, data = brtps_sim, model = regtps_s(s1, s2, k_s = 40),
                   engine = "laplace", num_samples = 50, seed = 1)
  expect_equal(coef(fit_chr), coef(fit_sym))
  lon <- brtps_sim$s1
  lat <- brtps_sim$s2
  fit_vec <- brtps(y ~ x1, data = brtps_sim, model = regtps_s(lon, lat, k_s = 40),
                   engine = "laplace", num_samples = 50, seed = 1)
  expect_equal(coef(fit_vec), coef(fit_sym))
})

test_that("brtps() validates its inputs", {
  expect_error(brtps(y ~ x1, data = brtps_sim), "regtps_s")
  expect_error(brtps(y ~ x1, data = as.list(brtps_sim), model = regtps_s(s1, s2)),
               "data frame")
  expect_error(brtps(~ x1, data = brtps_sim, model = regtps_s(s1, s2)), "response")
  expect_error(brtps(y ~ x1, data = brtps_sim, model = regtps_s(s1, s2), family = poisson(),
                     engine = "laplace"), "not available yet")
  expect_error(brtps(y ~ x1, data = brtps_sim, model = regtps_s(s1, s2),
                     priors = brtps_priors(tau = prior_pc()), engine = "laplace"),
               "Unknown prior")
  dat <- brtps_sim
  dat$x3 <- 2 * dat$x1
  expect_error(brtps(y ~ x1 + x3, data = dat, model = regtps_s(s1, s2), engine = "laplace"),
               "rank deficient")
})

test_that("a large-scale response triggers a note about the sigma prior", {
  dat <- brtps_sim
  dat$y <- dat$y * 20
  expect_message(fit_lap(data = dat), "unit scale")
  expect_no_message(fit_lap(data = dat, priors = brtps_priors(sigma = prior_pc(20, 0.05))))
})
