# student-t-internal.R -- the heavy-tailed (Student-t) observation model.
#
# A Student-t observation noise is the exact Gaussian scale-mixture
#   v_t | lambda_t ~ N(0, R / lambda_t),  lambda_t ~ Gamma(nu/2, nu/2),
# whose marginal v_t ~ t_nu(0, R) (Bishop 2006, sec. 2.3.7). Conditional on the
# latent scale lambda_t the update is an ordinary Kalman step with an inflated
# observation scale, so an outlier (large innovation) drives lambda_t down and
# is automatically down-weighted. The per-step scale and the state are coupled,
# so they are solved by the variational fixed point of Agamennoni, Nieto and
# Nebot (2012) / Sarkka and Nummenmaa (2009): alternate the Kalman update and the
# scale update until lambda_t converges. The recursion is deterministic -- no
# particle approximation -- and collapses to the Gaussian filter as nu -> Inf.
#
# References:
#   Agamennoni, G., Nieto, J. I. and Nebot, E. M. (2012). Approximate inference
#     in state-space models with heavy-tailed noise. IEEE Trans. Signal Process.
#     60(10), 5024-5037.
#   Sarkka, S. and Nummenmaa, A. (2009). Recursive noise adaptive Kalman
#     filtering by variational Bayesian approximations. IEEE Trans. Automat.
#     Control 54(3), 596-600.
#   Roth, M., Ozkan, E. and Gustafsson, F. (2013). A Student's t filter for heavy
#     tailed process and measurement noise. ICASSP, 5770-5774.

#' Multivariate Student-t log-density at a zero-mean point
#'
#' The log density of a `d`-variate Student-t with `df` degrees of freedom,
#' scale matrix `scale`, evaluated at `e` with location zero. Used as the
#' one-step-ahead predictive likelihood of the robust filter, with
#' `scale = B P_{t|t-1} B^T + R` the predictive t-scale.
#'
#' @param e A numeric vector of length `d`.
#' @param scale A `d` by `d` symmetric positive-definite scale matrix.
#' @param df The degrees of freedom (a finite scalar greater than zero).
#'
#' @returns A numeric scalar log-density.
#' @noRd
#' @keywords internal
.mvt_logpdf <- function(e, scale, df) {
  d <- length(e)
  sc <- .safe_chol(scale, "t scale")
  log_det <- 2 * sum(log(diag(sc)))
  quad <- sum(backsolve(sc, e, transpose = TRUE)^2)
  lgamma((df + d) / 2) - lgamma(df / 2) - 0.5 * d * log(df * pi) -
    0.5 * log_det - 0.5 * (df + d) * log1p(quad / df)
}

#' One robust (Student-t) Kalman update step
#'
#' Solves the variational fixed point for the latent scale `lambda` and returns
#' the robust filtered state. Each iteration runs an ordinary Kalman update with
#' observation scale `R / lambda`, then updates `lambda` from the expected
#' squared Mahalanobis residual under the t-scale `R` (the term that down-weights
#' outliers). The stored innovation covariance is the predictive t-scale
#' `S0 = B P_pred B^T + R` (`lambda = 1`), so the adequacy diagnostics whiten by
#' the right reference, and the returned log-likelihood is the Student-t
#' predictive density of the innovation.
#'
#' @param x_pred,p_pred The predicted state mean and covariance.
#' @param b The observation matrix at this step.
#' @param r The observation t-scale matrix \eqn{R}.
#' @param e The one-step-ahead innovation \eqn{y_t - B x_{t|t-1}}.
#' @param df The Student-t degrees of freedom.
#' @param tol,max_iter Fixed-point convergence tolerance and iteration cap.
#'
#' @returns A list with `x` (filtered mean), `p` (filtered covariance), `s` (the
#'   predictive t-scale innovation covariance) and `loglik` (the t predictive
#'   log-density of `e`).
#' @noRd
#' @keywords internal
.kalman_update_t <- function(x_pred, p_pred, b, r, e, df,
                             tol = 1e-8, max_iter = 50L) {
  d <- length(e)
  r_inv <- chol2inv(.safe_chol(r, "observation scale"))
  bpb <- b %*% p_pred %*% t(b)
  lambda <- 1
  x <- x_pred
  p <- p_pred
  for (it in seq_len(max_iter)) {
    s <- bpb + r / lambda
    s_inv <- chol2inv(.safe_chol(s, "innovation covariance"))
    gain <- p_pred %*% t(b) %*% s_inv
    x <- x_pred + as.numeric(gain %*% e)
    p <- .symmetrise(p_pred - gain %*% b %*% p_pred)
    # Expected squared Mahalanobis residual under the t-scale R: the data term
    # plus the state-uncertainty trace term (Agamennoni 2012, eq. 19).
    resid <- e - as.numeric(b %*% (x - x_pred))
    delta2 <- sum(resid * as.numeric(r_inv %*% resid)) +
      sum(diag(r_inv %*% b %*% p %*% t(b)))
    lambda_new <- (df + d) / (df + delta2)
    if (abs(lambda_new - lambda) < tol) {
      lambda <- lambda_new
      break
    }
    lambda <- lambda_new
  }
  # One final update at the converged scale so (x, p) match the last lambda.
  s <- bpb + r / lambda
  s_inv <- chol2inv(.safe_chol(s, "innovation covariance"))
  gain <- p_pred %*% t(b) %*% s_inv
  x <- x_pred + as.numeric(gain %*% e)
  p <- .symmetrise(p_pred - gain %*% b %*% p_pred)

  s0 <- bpb + r
  list(x = x, p = p, s = s0, loglik = .mvt_logpdf(e, s0, df))
}

#' Scale and family-fit diagnostics under a Student-t observation model
#'
#' The Gaussian innovation diagnostics test the standardised innovations against
#' a standard normal; under a Student-t observation family the correct reference
#' is the t itself, so this reads the two non-whiteness tests off the
#' probability integral transform (PIT). The per-step squared Mahalanobis
#' \eqn{m_t = w_t^\top w_t} satisfies \eqn{m_t / d \sim F(d, \nu)}, so its PIT is
#' uniform under correct specification: a shifted PIT mean signals a mis-stated
#' scale, and non-uniformity (Jarque-Bera on the inverse-normal PIT) signals a
#' wrong tail family. For a univariate series the signed innovation
#' \eqn{w_t \sim t_\nu} is used for the family-fit so skewness is caught too.
#'
#' @param w An `n` by `d` matrix of standardised innovations.
#' @param df The Student-t degrees of freedom.
#'
#' @returns A list with `scale_p` and `family_p`.
#' @noRd
#' @keywords internal
.t_scale_family_tests <- function(w, df) {
  n <- nrow(w)
  d <- ncol(w)
  clamp <- function(u) pmin(pmax(u, .Machine$double.eps), 1 - .Machine$double.eps)
  maha <- rowSums(w^2)
  u_rad <- clamp(stats::pf(maha / d, df1 = d, df2 = df))

  # Scale: the radial PIT has mean 0.5 under correct specification; a one-sample
  # t-test of that mean is valid for any df (no moment condition on the t).
  sd_u <- stats::sd(u_rad)
  scale_p <- if (is.finite(sd_u) && sd_u > 0) {
    tstat <- (mean(u_rad) - 0.5) / (sd_u / sqrt(n))
    2 * stats::pt(-abs(tstat), df = n - 1L)
  } else {
    NA_real_
  }

  # Family-fit: inverse-normal PIT should be standard normal. For a univariate
  # series the signed innovation is used (catches skew); otherwise the radial PIT.
  g <- if (d == 1L) {
    stats::qnorm(clamp(stats::pt(w[, 1L], df = df)))
  } else {
    stats::qnorm(u_rad)
  }
  list(scale_p = scale_p, family_p = .jarque_bera_p(g))
}

#' Clone an ssm with a different observation family
#'
#' Rebuilds an [ssm] reusing its system matrices and prior but with the given
#' observation family and degrees of freedom. Used by [estimate_obs_df()] to
#' evaluate the likelihood under each candidate `df`.
#'
#' @param model An [ssm].
#' @param family `"gaussian"` or `"student_t"`.
#' @param df The Student-t degrees of freedom (`Inf` for Gaussian).
#'
#' @returns An [ssm].
#' @noRd
#' @keywords internal
.with_family <- function(model, family, df) {
  ssm(
    transition = model@transition,
    observation = model@observation,
    state_cov = model@state_cov,
    obs_cov = model@obs_cov,
    init_state = model@init_state,
    init_cov = model@init_cov,
    obs_family = family,
    obs_df = df
  )
}
