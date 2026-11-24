# The paper's simulation locations (Appendix D): set.seed(1234 + i),
# N = 50 * i points on the unit square, k = floor(0.99 * N).
paper_locations <- function(i) {
  set.seed(1234 + i)
  n <- 50 * i
  data.frame(s1 = runif(n), s2 = runif(n))
}

test_that("truncation levels reproduce Table 1 of the paper", {
  expected <- c(45, 90, 131, 177)
  for (i in 1:4) {
    loc <- paper_locations(i)
    b <- regtps_basis(loc$s1, loc$s2, k = floor(0.99 * nrow(loc)),
                      m_trunc = 0.99, center = FALSE)
    expect_equal(b$M, expected[i])
    expect_equal(b$d_null, 3L)
  }
})

test_that("truncation levels reproduce Table 3 of the paper (scenario 3)", {
  loc <- paper_locations(3)
  M <- vapply(c(0.85, 0.90, 0.95, 0.99), function(g) {
    regtps_basis(loc$s1, loc$s2, k = 148, m_trunc = g, center = FALSE)$M
  }, numeric(1))
  expect_equal(M, c(83, 93, 107, 131))
})

test_that("basis pieces are consistent", {
  loc <- paper_locations(1)
  b <- regtps_basis(loc$s1, loc$s2, k = 40)
  expect_equal(dim(b$Phi), c(50L, b$M))
  expect_length(b$v, b$M)
  expect_true(all(b$v[1:3] == 0))
  expect_true(all(b$v[-(1:3)] > 0))
  expect_true(all(diff(b$v_all) >= 0))
  expect_equal(b$cumvar[length(b$cumvar)], 1)
  # Gamma has orthonormal columns
  expect_equal(crossprod(b$Gamma), diag(b$M), tolerance = 1e-10)
  # initial alpha follows eq. A22, clipped to [-3, 3] on the log scale
  expect_equal(b$log_alpha_init, min(max(log(1 / median(b$v_all[b$v_all > 0])), -3), 3))
})

test_that("an integer truncation keeps that many components", {
  loc <- paper_locations(1)
  b <- regtps_basis(loc$s1, loc$s2, k = 40, m_trunc = 25)
  expect_equal(b$M, 25)
  expect_equal(ncol(b$Phi), 25)
  expect_equal(regtps_basis(loc$s1, loc$s2, k = 40, m_trunc = 1)$M, 40)
})

test_that("the centring constraint removes the constant from the null space", {
  loc <- paper_locations(1)
  b <- regtps_basis(loc$s1, loc$s2, k = 40, center = TRUE)
  expect_equal(b$d_null, 2L)
  expect_equal(b$ncol_basis, 39L)
  expect_equal(colSums(b$Phi), rep(0, b$M), tolerance = 1e-8)
})

test_that("invalid inputs give clear errors", {
  loc <- paper_locations(1)
  expect_error(regtps_basis(loc$s1, loc$s2, k = 60), "larger than the number of distinct")
  expect_error(regtps_basis(loc$s1, loc$s2, k = 40, m_trunc = 2), "fewer components")
  expect_error(regtps_basis(loc$s1, loc$s2, k = 40, m_trunc = 50), "larger than the basis")
  expect_error(regtps_basis(loc$s1, loc$s2[-1]), "same length")
  expect_error(regtps_s(s1), "both coordinates")
  expect_error(regtps_s(s1, s2, k_s = 2.5), "whole number")
  expect_error(regtps_s(s1, s2, m_s_trunc = 1.5), "proportion")
  expect_error(regtps_s(s1, s2, center = "yes"), "NULL, TRUE or FALSE")
})
