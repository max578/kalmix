# linalg-internal.R -- small numerically-careful matrix helpers.
#
# The filtering and smoothing recursions repeatedly Cholesky-factor a
# covariance and re-symmetrise the result of a matrix product. These helpers
# isolate that arithmetic with a clear error when a covariance has drifted
# non-positive-definite, and a log-Gaussian-density evaluator the mixture and
# HMM filters share.

# Matrix guards ---------------------------------------------------------------

#' Cholesky factor of a covariance, with a clear error on failure
#'
#' Returns the upper-triangular Cholesky factor of a symmetric positive-definite
#' matrix. A tiny jitter is added to the diagonal first so that a covariance
#' that is positive semi-definite only up to rounding still factors; a genuine
#' breakdown raises an informative error rather than the opaque base message.
#'
#' @param x A symmetric positive (semi-)definite matrix.
#' @param what Character scalar — what the matrix represents, for the error.
#'
#' @returns The upper-triangular Cholesky factor.
#' @noRd
#' @keywords internal
.safe_chol <- function(x, what = "covariance") {
  jitter <- diag(.Machine$double.eps^0.5 * (mean(diag(x)) + 1), nrow = nrow(x))
  out <- tryCatch(chol(x + jitter), error = function(e) NULL)
  if (is.null(out)) {
    cli::cli_abort(c(
      "The {what} is not positive definite and cannot be factored.",
      "i" = "Check the model's covariance matrices are positive definite."
    ))
  }
  out
}

#' Force a matrix to be exactly symmetric
#'
#' Averages a matrix with its transpose to remove the asymmetry that
#' accumulates from finite-precision products in the covariance updates.
#'
#' @param x A numeric matrix expected to be symmetric.
#'
#' @returns The symmetrised matrix.
#' @noRd
#' @keywords internal
.symmetrise <- function(x) {
  (x + t(x)) / 2
}

# Densities and weights -------------------------------------------------------

#' Log-density of a multivariate Gaussian
#'
#' Evaluates \eqn{\log N(x; \mu, \Sigma)} through the Cholesky factor of
#' \eqn{\Sigma}, avoiding an explicit inverse. Shared by the mixture filter
#' (component evidence) and the hidden Markov model filter (emission density).
#'
#' @param x Numeric vector — the point to evaluate.
#' @param mean Numeric vector — the mean \eqn{\mu}.
#' @param sigma Numeric matrix — the covariance \eqn{\Sigma}.
#'
#' @returns Numeric scalar — the log-density.
#' @noRd
#' @keywords internal
.dmvnorm_log <- function(x, mean, sigma) {
  d <- length(x)
  chol_sigma <- .safe_chol(sigma, "observation covariance")
  resid <- backsolve(chol_sigma, x - mean, transpose = TRUE)
  log_det <- 2 * sum(log(diag(chol_sigma)))
  -0.5 * (d * log(2 * pi) + log_det + sum(resid^2))
}

#' Normalise a vector of log-weights to probabilities
#'
#' Exponentiates and normalises log-weights with the standard log-sum-exp shift
#' for numerical stability, returning both the probabilities and the
#' log-normaliser (which the mixture and HMM filters accumulate into the
#' log-likelihood).
#'
#' @param log_w Numeric vector of unnormalised log-weights.
#'
#' @returns A list with `prob` (the normalised probabilities) and `log_norm`
#'   (the log of the summed unnormalised weights).
#' @noRd
#' @keywords internal
.softmax_with_norm <- function(log_w) {
  shift <- max(log_w)
  unnorm <- exp(log_w - shift)
  total <- sum(unnorm)
  list(prob = unnorm / total, log_norm = shift + log(total))
}
