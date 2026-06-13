# kalman.R -- the linear-Gaussian Kalman filter and RTS smoother.
#
# kalman_filter() runs the closed-form forward recursion of Kalman (1960) over
# an ssm: a predict step that propagates the state mean and covariance through
# the transition, then an update step that corrects them with each observation
# through the Kalman gain, accumulating the prediction-error-decomposition
# log-likelihood. rts_smoother() runs the Rauch-Tung-Striebel (1965) backward
# recursion over a completed filter pass to obtain the smoothed states given
# the whole series. Both are exact for a linear-Gaussian model.
#
# Assimilative Causal Inference (aci.R) reads off this filter/smoother pair: the
# time-resolved causal information is KL( smoother_t || filter_t ) per step
# (Andreou, Chen and Bollt 2026, eq. 7), and the objective causal information
# rate gives a decision lead-time (their eqs. 8-9). The Gaussian KL is a closed
# form implemented natively in aci.R; when kernR is installed its
# relative_entropy() is used as an independent oracle in the test suite.

#' A Kalman filter pass
#'
#' An S7 object holding the output of a forward Kalman pass over an [ssm], as
#' returned by [kalman_filter()]. The filtered quantities condition on
#' observations up to and including each step; the predicted quantities
#' condition on observations strictly before it. The innovations and their
#' covariances are retained because the log-likelihood and the
#' interrupted-time-series contrast are read off them.
#'
#' @param filtered_mean An `n` by `m` matrix of filtered state means, row `t`
#'   being \eqn{E[x_t \mid y_{1:t}]}.
#' @param filtered_cov A list of `n` filtered state covariances.
#' @param predicted_mean An `n` by `m` matrix of one-step-ahead predicted state
#'   means, row `t` being \eqn{E[x_t \mid y_{1:t-1}]}.
#' @param predicted_cov A list of `n` one-step-ahead predicted state
#'   covariances.
#' @param innovation An `n` by `d` matrix of innovations (prediction errors).
#' @param innovation_cov A list of `n` innovation covariances.
#' @param log_lik Numeric scalar — the model log-likelihood from the prediction
#'   error decomposition.
#' @param model The [ssm] that was filtered.
#'
#' @returns An S7 object of class `kalman_fit`.
#' @family state-space
#' @export
kalman_fit <- S7::new_class(
  name = "kalman_fit",
  package = "kalmix",
  properties = list(
    filtered_mean = S7::class_double,
    filtered_cov = S7::class_list,
    predicted_mean = S7::class_double,
    predicted_cov = S7::class_list,
    innovation = S7::class_double,
    innovation_cov = S7::class_list,
    log_lik = S7::class_double,
    model = ssm
  )
)

#' A Rauch-Tung-Striebel smoother pass
#'
#' An S7 object holding the output of a backward RTS smoothing pass, as
#' returned by [rts_smoother()]. Each smoothed quantity conditions on the whole
#' observed series, \eqn{E[x_t \mid y_{1:n}]} and its covariance.
#'
#' @param smoothed_mean An `n` by `m` matrix of smoothed state means.
#' @param smoothed_cov A list of `n` smoothed state covariances.
#' @param log_lik Numeric scalar — the model log-likelihood (carried through
#'   from the filter pass for convenience).
#' @param model The [ssm] that was smoothed.
#'
#' @returns An S7 object of class `rts_fit`.
#' @family state-space
#' @export
rts_fit <- S7::new_class(
  name = "rts_fit",
  package = "kalmix",
  properties = list(
    smoothed_mean = S7::class_double,
    smoothed_cov = S7::class_list,
    log_lik = S7::class_double,
    model = ssm
  )
)

#' Kalman filter for a linear-Gaussian state-space model
#'
#' Runs the closed-form forward recursion of Kalman (1960) over an [ssm].
#' Starting from the model's prior, each step predicts the state through the
#' transition and then updates it with the observation through the Kalman gain
#' \deqn{K_t = P_{t|t-1} B_t^\top S_t^{-1}, \qquad
#'       S_t = B_t P_{t|t-1} B_t^\top + R_t,}
#' where \eqn{S_t} is the innovation covariance. The model log-likelihood is
#' accumulated from the prediction error decomposition,
#' \deqn{\log p(y_{1:n}) = -\tfrac12 \sum_t \left( d \log 2\pi +
#'       \log\lvert S_t\rvert + e_t^\top S_t^{-1} e_t \right),}
#' with innovation \eqn{e_t = y_t - B_t x_{t|t-1}}.
#'
#' The recursion is exact for a linear-Gaussian model. Time-varying system
#' matrices are honoured: a model whose `transition`, `observation`,
#' `state_cov` or `obs_cov` was supplied as a list of `n` matrices uses the
#' step-`t` matrix at step `t`.
#'
#' @section Assimilative causal inference:
#' The filtering distribution \eqn{N(x_{t|t}, P_{t|t})} returned here is one of
#' the two distributions [aci()] compares step by step: the time-resolved causal
#' information is its relative entropy from the smoothing distribution of
#' [rts_smoother()] (Andreou, Chen and Bollt 2026, eq. 7). See [aci()] and
#' [causal_information_rate()].
#'
#' @param model An [ssm].
#' @param y A numeric vector (univariate) or `n` by `d` matrix (multivariate)
#'   of observations.
#'
#' @returns A [kalman_fit].
#' @family state-space
#' @seealso [rts_smoother()] for the backward smoothing pass; [ssm()] for the
#'   model.
#' @references
#' Kalman, R. E. (1960). A new approach to linear filtering and prediction
#' problems. *Journal of Basic Engineering*, 82(1), 35--45.
#' @export
#' @examples
#' ## Recover a smooth signal from a noisy local-level series.
#' model <- ssm(
#'   transition = 1, observation = 1,
#'   state_cov = 0.01, obs_cov = 1,
#'   init_state = 0, init_cov = 10
#' )
#' set.seed(1)
#' level <- cumsum(rnorm(100, sd = 0.1))
#' y <- level + rnorm(100)
#' fit <- kalman_filter(model, y)
#' fit@log_lik
kalman_filter <- function(model, y) {
  if (!S7::S7_inherits(model, ssm)) {
    cli::cli_abort("`model` must be an {.cls ssm}.")
  }
  y <- .check_observations(y, model)

  n <- nrow(y)
  m <- model@state_dim
  d <- model@obs_dim

  filtered_mean <- matrix(0, nrow = n, ncol = m)
  predicted_mean <- matrix(0, nrow = n, ncol = m)
  innovation <- matrix(0, nrow = n, ncol = d)
  filtered_cov <- vector("list", n)
  predicted_cov <- vector("list", n)
  innovation_cov <- vector("list", n)

  x <- model@init_state
  p <- model@init_cov
  log_lik <- 0
  const <- d * log(2 * pi)

  for (t in seq_len(n)) {
    a <- .at_step(model@transition, t)
    b <- .at_step(model@observation, t)
    q <- .at_step(model@state_cov, t)
    r <- .at_step(model@obs_cov, t)

    # Predict ----------------------------------------------------------------

    x_pred <- as.numeric(a %*% x)
    p_pred <- a %*% p %*% t(a) + q
    predicted_mean[t, ] <- x_pred
    predicted_cov[[t]] <- p_pred

    # Update -----------------------------------------------------------------

    e <- y[t, ] - as.numeric(b %*% x_pred)
    s <- b %*% p_pred %*% t(b) + r
    s_chol <- .safe_chol(s, "innovation covariance")
    s_inv <- chol2inv(s_chol)
    gain <- p_pred %*% t(b) %*% s_inv

    x <- x_pred + as.numeric(gain %*% e)
    p <- p_pred - gain %*% b %*% p_pred
    p <- .symmetrise(p)

    filtered_mean[t, ] <- x
    filtered_cov[[t]] <- p
    innovation[t, ] <- e
    innovation_cov[[t]] <- s

    # Prediction error decomposition -----------------------------------------

    log_det <- 2 * sum(log(diag(s_chol)))
    quad <- sum(e * as.numeric(s_inv %*% e))
    log_lik <- log_lik - 0.5 * (const + log_det + quad)
  }

  kalman_fit(
    filtered_mean = filtered_mean,
    filtered_cov = filtered_cov,
    predicted_mean = predicted_mean,
    predicted_cov = predicted_cov,
    innovation = innovation,
    innovation_cov = innovation_cov,
    log_lik = log_lik,
    model = model
  )
}

#' Rauch-Tung-Striebel smoother for a linear-Gaussian state-space model
#'
#' Runs the backward recursion of Rauch, Tung and Striebel (1965) over a
#' completed Kalman pass to obtain the state estimates conditioned on the whole
#' series, \eqn{E[x_t \mid y_{1:n}]}. From the last filtered state backwards,
#' each step blends the filtered estimate with the smoothed estimate of the
#' next state through the smoother gain
#' \deqn{C_t = P_{t|t} A_{t+1}^\top P_{t+1|t}^{-1},}
#' giving \eqn{x_{t|n} = x_{t|t} + C_t (x_{t+1|n} - x_{t+1|t})} and the matching
#' covariance update. The result is exact for a linear-Gaussian model.
#'
#' Either a fitted [kalman_fit] or an [ssm] plus observations may be supplied;
#' in the latter case the filter pass is run first.
#'
#' @section Assimilative causal inference:
#' This smoother supplies the second of the two Gaussian distributions [aci()]
#' compares: the time-resolved causal information is the relative entropy of the
#' smoothing distribution \eqn{N(x_{t|n}, P_{t|n})} from the filtering
#' distribution (Andreou, Chen and Bollt 2026, eq. 7). See [aci()].
#'
#' @param object A [kalman_fit] from [kalman_filter()], or an [ssm].
#' @param y A numeric vector or matrix of observations. Required when `object`
#'   is an [ssm]; ignored when `object` is already a [kalman_fit].
#'
#' @returns An [rts_fit].
#' @family state-space
#' @seealso [kalman_filter()].
#' @references
#' Rauch, H. E., Tung, F. and Striebel, C. T. (1965). Maximum likelihood
#' estimates of linear dynamic systems. *AIAA Journal*, 3(8), 1445--1450.
#' @export
#' @examples
#' model <- ssm(
#'   transition = 1, observation = 1,
#'   state_cov = 0.01, obs_cov = 1,
#'   init_state = 0, init_cov = 10
#' )
#' set.seed(1)
#' y <- cumsum(rnorm(100, sd = 0.1)) + rnorm(100)
#' smoothed <- rts_smoother(model, y)
#' head(smoothed@smoothed_mean)
rts_smoother <- function(object, y = NULL) {
  fit <- if (S7::S7_inherits(object, kalman_fit)) {
    object
  } else if (S7::S7_inherits(object, ssm)) {
    if (is.null(y)) {
      cli::cli_abort("`y` is required when `object` is an {.cls ssm}.")
    }
    kalman_filter(object, y)
  } else {
    cli::cli_abort("`object` must be a {.cls kalman_fit} or an {.cls ssm}.")
  }

  model <- fit@model
  n <- nrow(fit@filtered_mean)
  m <- model@state_dim

  smoothed_mean <- fit@filtered_mean
  smoothed_cov <- fit@filtered_cov

  # Backward recursion from the penultimate step ------------------------------

  for (t in rev(seq_len(n - 1L))) {
    a_next <- .at_step(model@transition, t + 1L)
    p_pred_next <- fit@predicted_cov[[t + 1L]]
    p_filt <- fit@filtered_cov[[t]]

    gain <- p_filt %*% t(a_next) %*% chol2inv(
      .safe_chol(p_pred_next, "predicted covariance")
    )

    mean_gap <- smoothed_mean[t + 1L, ] - fit@predicted_mean[t + 1L, ]
    smoothed_mean[t, ] <- fit@filtered_mean[t, ] +
      as.numeric(gain %*% mean_gap)

    cov_gap <- smoothed_cov[[t + 1L]] - p_pred_next
    smoothed_cov[[t]] <- .symmetrise(
      p_filt + gain %*% cov_gap %*% t(gain)
    )
  }

  rts_fit(
    smoothed_mean = smoothed_mean,
    smoothed_cov = smoothed_cov,
    log_lik = fit@log_lik,
    model = model
  )
}

#' @export
S7::method(print, kalman_fit) <- function(x, ...) {
  n <- nrow(x@filtered_mean)
  cat(sprintf(
    "<kalman_fit>: %d step%s, state dim %d, log-likelihood %.3f\n",
    n, if (n == 1L) "" else "s", ncol(x@filtered_mean), x@log_lik
  ))
  invisible(x)
}

#' @export
S7::method(print, rts_fit) <- function(x, ...) {
  n <- nrow(x@smoothed_mean)
  cat(sprintf(
    "<rts_fit>: %d step%s, state dim %d (smoothed)\n",
    n, if (n == 1L) "" else "s", ncol(x@smoothed_mean)
  ))
  invisible(x)
}
