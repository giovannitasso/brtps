#' Spatial regTPS component
#'
#' Specifies the spatial random field of a [brtps()] model. The field is
#' approximated by the regularized thin plate spline (regTPS) construction of
#' Cavieres and Krumscheid (2026): a thin plate regression spline basis of
#' dimension `k_s` is built with \pkg{mgcv}, the penalty matrix \eqn{S} is
#' eigendecomposed, and the expansion is truncated to the first \eqn{M}
#' eigenfunctions.
#'
#' @param s1,s2 Spatial coordinates. Usually the (unquoted) names of two
#'   columns of `data` in [brtps()]; a character string or a numeric vector
#'   also works. Use column names if you want to call [predict.brtps()] on new
#'   locations.
#' @param k_s Number of thin plate spline basis functions (`K` in the paper,
#'   `k` in [mgcv::s()]). It must not exceed the number of distinct locations.
#'   The default, `NULL`, uses `min(floor(0.99 * n), 200)`, where `n` is the
#'   number of distinct locations (the paper uses `floor(0.99 * n)` for small
#'   data sets).
#' @param m_s_trunc Truncation of the expansion. A number in `(0, 1]` is the
#'   share of the total prior variance to retain (\eqn{\gamma} in eq. 39 of the
#'   paper), and an integer of 2 or more is the number of retained components
#'   \eqn{M} directly. Default `0.99`.
#' @param center Identifiability between the formula intercept and the
#'   constant mode of the spline null space. `NULL` (default) applies a
#'   sum-to-zero constraint to the field when the formula has an intercept,
#'   as \pkg{mgcv} does, and none otherwise; then the null space is
#'   \{1, s1, s2\} as in the paper. `TRUE` or `FALSE` force either choice.
#'
#' @details
#' With `m_s_trunc` as a proportion, \eqn{M} is the smallest number of
#' components whose prior variances \eqn{\lambda_k = 1 / (1 + \alpha v_k)}
#' add up to that share of the total, evaluated at
#' \eqn{\alpha = 1 / \mathrm{median}(v_k > 0)} (eq. A22). With this rule and
#' `center = FALSE` the package reproduces the truncation levels reported in
#' Tables 1 and 3 of the paper.
#'
#' @return An object of class `regtps_s` to pass to the `model` argument of
#'   [brtps()].
#' @seealso [regtps_basis()] to inspect the basis and its spectrum.
#' @export
#' @examples
#' regtps_s(s1, s2, k_s = 100, m_s_trunc = 0.99)
regtps_s <- function(s1, s2, k_s = NULL, m_s_trunc = 0.99, center = NULL) {
  if (missing(s1) || missing(s2)) {
    stop("regtps_s() needs both coordinates, 's1' and 's2'.", call. = FALSE)
  }
  check_k(k_s)
  check_trunc(m_s_trunc)
  if (!is.null(center) && !(is.logical(center) && length(center) == 1 && !is.na(center))) {
    stop("'center' must be NULL, TRUE or FALSE.", call. = FALSE)
  }
  structure(
    list(
      s1 = substitute(s1),
      s2 = substitute(s2),
      env = parent.frame(),
      k_s = k_s,
      m_s_trunc = m_s_trunc,
      center = center
    ),
    class = c("regtps_s", "brtps_model")
  )
}

#' @export
print.regtps_s <- function(x, ...) {
  cat("Spatial regTPS component\n")
  cat("  coordinates: ", deparse(x$s1), ", ", deparse(x$s2), "\n", sep = "")
  cat("  k_s:         ", if (is.null(x$k_s)) "default" else x$k_s, "\n", sep = "")
  cat("  m_s_trunc:   ", x$m_s_trunc, "\n", sep = "")
  cat("  center:      ", if (is.null(x$center)) "automatic" else x$center, "\n", sep = "")
  invisible(x)
}

#' regTPS basis and spectrum
#'
#' Builds the truncated regTPS basis used by [regtps_s()]: the thin plate
#' regression spline basis \eqn{\tilde\Phi} and penalty \eqn{S} from
#' \pkg{mgcv}, the eigendecomposition \eqn{S = \Gamma V \Gamma^\top} (eq. 33)
#' and the design matrix \eqn{\Phi = \tilde\Phi \Gamma_M} of the retained
#' components (eq. 53).
#'
#' @param s1,s2 Numeric vectors of coordinates.
#' @param k Number of basis functions (see `k_s` in [regtps_s()]).
#' @param m_trunc Truncation (see `m_s_trunc` in [regtps_s()]).
#' @param center Logical; apply a sum-to-zero constraint to the field.
#' @param alpha_ref Reference value \eqn{\alpha_{ref}} in eq. A22.
#'
#' @return An object of class `regtps_basis` with, among others, `Phi` (the
#'   \eqn{N \times M} design matrix), `v` (the \eqn{M} retained penalty
#'   eigenvalues, zero for the null space), `M`, `k`, `d_null`, `cumvar` (the
#'   cumulative share of prior variance of all components) and the
#'   \pkg{mgcv} smooth object used for prediction. `plot()` shows the
#'   cumulative variance and the eigenvalue spectrum (Figures 3 and 4 of the
#'   paper).
#' @export
#' @examples
#' set.seed(1)
#' b <- regtps_basis(runif(100), runif(100), k = 60)
#' b
#' plot(b)
regtps_basis <- function(s1, s2, k = NULL, m_trunc = 0.99, center = FALSE,
                         alpha_ref = 1) {
  if (!is.numeric(s1) || !is.numeric(s2) || length(s1) != length(s2)) {
    stop("'s1' and 's2' must be numeric vectors of the same length.", call. = FALSE)
  }
  if (anyNA(s1) || anyNA(s2)) stop("Coordinates contain missing values.", call. = FALSE)
  check_k(k)
  check_trunc(m_trunc)
  n_unique <- nrow(unique(cbind(s1, s2)))
  if (is.null(k)) k <- max(10L, min(floor(0.99 * n_unique), 200L))
  k <- as.integer(k)
  if (k > n_unique) {
    stop("k = ", k, " is larger than the number of distinct locations (",
         n_unique, ").", call. = FALSE)
  }

  dat <- data.frame(s1 = s1, s2 = s2)
  sm <- mgcv::smoothCon(mgcv::s(s1, s2, bs = "tp", k = k), data = dat,
                        knots = NULL, absorb.cons = center,
                        scale.penalty = TRUE)[[1]]
  Xtilde <- sm$X
  S <- sm$S[[1]]

  eig <- eigen(S, symmetric = TRUE)
  ord <- order(eig$values)                     # ascending, eq. 34
  v_all <- eig$values[ord]
  Gamma <- eig$vectors[, ord, drop = FALSE]

  # Null space: thin plate splines of order 2 in 2-D have {1, s1, s2}; the
  # sum-to-zero constraint removes the constant.
  d_null <- if (isTRUE(center)) 2L else 3L
  v_all[seq_len(d_null)] <- 0

  alpha_trunc <- alpha_ref / median(v_all[v_all > 0])   # eq. A22
  lambda <- 1 / (1 + alpha_trunc * v_all)                # eq. 38
  cumvar <- cumsum(lambda) / sum(lambda)

  if (m_trunc <= 1) {
    M <- if (m_trunc >= 1) length(v_all) else which(cumvar >= m_trunc)[1]  # eq. 39
  } else {
    M <- as.integer(m_trunc)
    if (M > length(v_all)) {
      stop("m_trunc = ", M, " is larger than the basis dimension (",
           length(v_all), ").", call. = FALSE)
    }
  }
  if (M < d_null) {
    stop("The truncation keeps fewer components (", M, ") than the null space (",
         d_null, ").", call. = FALSE)
  }
  keep <- seq_len(M)

  structure(
    list(
      Phi = Xtilde %*% Gamma[, keep, drop = FALSE],
      v = v_all[keep],
      Gamma = Gamma[, keep, drop = FALSE],
      v_all = v_all,
      lambda = lambda,
      cumvar = cumvar,
      M = M,
      k = k,
      ncol_basis = ncol(Xtilde),
      d_null = d_null,
      center = isTRUE(center),
      m_trunc = m_trunc,
      alpha_ref = alpha_ref,
      alpha_trunc = alpha_trunc,
      log_alpha_init = min(max(log(alpha_trunc), -3), 3),
      n_unique = n_unique,
      smooth = sm
    ),
    class = "regtps_basis"
  )
}

#' @export
print.regtps_basis <- function(x, ...) {
  cat("regTPS basis\n")
  cat("  basis functions (k):   ", x$k, "\n", sep = "")
  cat("  retained components (M): ", x$M, " (", round(100 * x$cumvar[x$M], 1),
      "% of the prior variance)\n", sep = "")
  cat("  null space:            ", x$d_null,
      if (x$center) " (sum-to-zero constraint)" else "", "\n", sep = "")
  invisible(x)
}

#' @export
plot.regtps_basis <- function(x, ...) {
  op <- graphics::par(mfrow = c(1, 2))
  on.exit(graphics::par(op))
  K <- length(x$cumvar)
  graphics::plot(seq_len(K), x$cumvar, type = "l", xlab = "Number of components",
                 ylab = "Cumulative prior variance", main = "Variance explained", ...)
  graphics::abline(v = x$M, lty = 2)
  if (x$m_trunc <= 1) graphics::abline(h = x$m_trunc, lty = 3)
  pos <- x$v_all > 0
  graphics::plot(which(pos), x$v_all[pos], log = "y", xlab = "Component k",
                 ylab = expression(v[k]), main = "Penalty eigenvalues", pch = 20, ...)
  graphics::abline(v = x$M, lty = 2)
  invisible(x)
}

# Builds the model-specific pieces for brtps(). New model types (for example a
# spatio-temporal regtps_st) add a build_model() and a finish_model() method.
build_model <- function(model, data, has_intercept) UseMethod("build_model")

finish_model <- function(model, built, rows) UseMethod("finish_model")

#' @export
build_model.regtps_s <- function(model, data, has_intercept) {
  s1 <- eval_coord(model$s1, data, model$env, "s1")
  s2 <- eval_coord(model$s2, data, model$env, "s2")
  if (length(s1) != nrow(data) || length(s2) != nrow(data)) {
    stop("The coordinates must have one value per row of 'data'.", call. = FALSE)
  }
  list(
    tmb_model = "regtps_s",
    coords = data.frame(s1 = s1, s2 = s2),
    coord_names = c(coord_name(model$s1), coord_name(model$s2)),
    center = if (is.null(model$center)) has_intercept else model$center
  )
}

# Builds the basis once the rows used in the fit are known.
#' @export
finish_model.regtps_s <- function(model, built, rows) {
  coords <- built$coords[rows, , drop = FALSE]
  basis <- regtps_basis(coords$s1, coords$s2, k = model$k_s,
                        m_trunc = model$m_s_trunc, center = built$center)
  built$coords <- coords
  built$basis <- basis
  built$tmb_data <- list(Phi = basis$Phi, v = basis$v)
  built
}

# Default priors of the model's own hyperparameters.
model_default_priors <- function(model) UseMethod("model_default_priors")

#' @export
model_default_priors.regtps_s <- function(model) {
  list(alpha = prior_lognormal(0, 1))
}

# Design matrix of the field at new locations (or at the observed ones).
field_design <- function(model, built, newdata) UseMethod("field_design")

#' @export
field_design.regtps_s <- function(model, built, newdata) {
  if (is.null(newdata)) return(built$basis$Phi)
  coords <- lapply(1:2, function(i) {
    name <- built$coord_names[i]
    if (is.na(name)) name <- c("s1", "s2")[i]
    if (!name %in% names(newdata)) {
      stop("'newdata' must contain the coordinate column '", name, "'.", call. = FALSE)
    }
    as.numeric(newdata[[name]])
  })
  if (anyNA(coords[[1]]) || anyNA(coords[[2]])) {
    stop("The coordinates in 'newdata' contain missing values.", call. = FALSE)
  }
  Xt <- mgcv::PredictMat(built$basis$smooth,
                         data.frame(s1 = coords[[1]], s2 = coords[[2]]))
  Xt %*% built$basis$Gamma
}

eval_coord <- function(expr, data, env, what) {
  value <- eval(expr, data, env)
  if (is.character(value) && length(value) == 1) {
    if (!value %in% names(data)) {
      stop("Column '", value, "' (", what, ") was not found in 'data'.", call. = FALSE)
    }
    value <- data[[value]]
  }
  if (!is.numeric(value)) stop("Coordinate '", what, "' must be numeric.", call. = FALSE)
  as.numeric(value)
}

coord_name <- function(expr) {
  if (is.symbol(expr)) return(as.character(expr))
  if (is.character(expr) && length(expr) == 1) return(expr)
  NA_character_
}

check_k <- function(k) {
  if (!is.null(k) && !(is.numeric(k) && length(k) == 1 && !is.na(k) && k >= 4 &&
                       k == round(k))) {
    stop("'k_s' must be NULL or a whole number of at least 4.", call. = FALSE)
  }
}

check_trunc <- function(m) {
  ok <- is.numeric(m) && length(m) == 1 && !is.na(m) && m > 0 &&
    (m <= 1 || m == round(m))
  if (!ok) {
    stop("'m_s_trunc' must be a proportion in (0, 1] or a whole number of ",
         "components.", call. = FALSE)
  }
}
