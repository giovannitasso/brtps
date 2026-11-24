# Simulated data set shipped with brtps (data/brtps_sim.rda).
#
# 150 locations on the unit square, a Gaussian random field with Matern
# covariance (nu = 1, range rho = 0.3, sigma_u = 1) as in Simulation 2 of
# Cavieres and Krumscheid (2026), two covariates and Gaussian noise
# (sigma_e = 0.3):
#   y = 1 + 0.5 * x1 - 0.8 * x2 + u(s) + e
# Run from the package root with: Rscript data-raw/brtps_sim.R

matern_cov <- function(coords, nu, rho, sigma2) {
  D <- as.matrix(dist(coords))
  D[D == 0] <- 1e-10
  scaling <- (sqrt(2 * nu) * D) / rho
  C <- (2^(1 - nu)) / gamma(nu) * scaling^nu * besselK(scaling, nu)
  diag(C) <- 1
  sigma2 * C
}

set.seed(2026)
n <- 150
s1 <- runif(n)
s2 <- runif(n)
C <- matern_cov(cbind(s1, s2), nu = 1, rho = 0.3, sigma2 = 1)
u <- drop(t(chol(C)) %*% rnorm(n))
x1 <- rnorm(n)
x2 <- runif(n, -1, 1)
y <- 1 + 0.5 * x1 - 0.8 * x2 + u + rnorm(n, 0, 0.3)

brtps_sim <- data.frame(s1 = s1, s2 = s2, x1 = x1, x2 = x2, y = y, u = u)
save(brtps_sim, file = "data/brtps_sim.rda", compress = "xz")
