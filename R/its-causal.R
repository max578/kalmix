# its-causal.R -- single-case interrupted-time-series causal contrast.
#
# its_causal() is the minimal state-space causal lead: it asks what would have
# happened to a single series had an intervention at a known time not occurred.
# A local-level (or local-linear-trend) state-space model is fitted to the
# pre-intervention segment, projected forward to form the counterfactual the
# series would have followed, and the observed post-intervention path is
# contrasted with that counterfactual. The pointwise and cumulative effects
# carry predictive intervals from the forecast covariance, so a null effect is
# reported honestly rather than asserted. This is the single-case
# interrupted-time-series design done on the Kalman forecast.

#' A single-case interrupted-time-series causal contrast
#'
#' An S7 object holding the output of [its_causal()]: the observed
#' post-intervention series, the model's counterfactual forecast and its
#' predictive interval, the pointwise and cumulative effects, and the average
#' post-intervention effect with an interval.
#'
#' @param time A numeric or integer vector of the post-intervention time
#'   indices.
#' @param observed A numeric vector of the observed post-intervention values.
#' @param counterfactual A numeric vector of the forecast counterfactual.
#' @param counterfactual_lower,counterfactual_upper Numeric vectors — the
#'   predictive-interval bounds of the counterfactual.
#' @param pointwise_effect A numeric vector — observed minus counterfactual.
#' @param cumulative_effect A numeric vector — the running sum of the pointwise
#'   effect.
#' @param avg_effect Numeric scalar — the average post-intervention effect.
#' @param avg_effect_lower,avg_effect_upper Numeric scalars — the interval on
#'   the average effect.
#' @param level Numeric scalar — the interval coverage used.
#'
#' @returns An S7 object of class `its_fit`.
#' @family causal
#' @export
its_fit <- S7::new_class(
  name = "its_fit",
  package = "kalmix",
  properties = list(
    time = S7::class_double,
    observed = S7::class_double,
    counterfactual = S7::class_double,
    counterfactual_lower = S7::class_double,
    counterfactual_upper = S7::class_double,
    pointwise_effect = S7::class_double,
    cumulative_effect = S7::class_double,
    avg_effect = S7::class_double,
    avg_effect_lower = S7::class_double,
    avg_effect_upper = S7::class_double,
    level = S7::class_double
  )
)

#' Single-case interrupted-time-series causal effect
#'
#' Estimates the causal effect of an intervention on a single time series by
#' contrasting the observed post-intervention path with the counterfactual a
#' state-space model forecasts from the pre-intervention segment. A local-level
#' or local-linear-trend [ssm] is fitted to `x[1:(intervention - 1)]`, its
#' state is projected forward across the post-intervention window, and the
#' observed values are differenced against that forecast. The forecast
#' covariance gives a predictive interval, so the pointwise, cumulative and
#' average effects each carry an honest uncertainty band -- a genuinely null
#' effect shows an interval that straddles zero.
#'
#' The design is the single-case (N-of-1) interrupted time series; it is causal
#' only under the maintained assumption that, absent the intervention, the
#' pre-period dynamics would have continued. That assumption is the user's to
#' defend; the function quantifies the contrast it implies.
#'
#' @param x A numeric vector — the full series, pre- and post-intervention.
#' @param intervention Integer scalar — the index of the first
#'   post-intervention observation. The pre-period is `x[1:(intervention - 1)]`.
#' @param trend Logical scalar — whether to use a local-linear-trend model
#'   (`TRUE`) rather than a local level (`FALSE`). A trend lets the
#'   counterfactual continue a pre-period slope. Defaults to `FALSE`.
#' @param level Numeric scalar in `(0, 1)` — the predictive-interval coverage.
#'   Defaults to `0.95`.
#'
#' @returns An [its_fit].
#' @family causal
#' @seealso [ssm()], [kalman_filter()].
#' @references
#' Brodersen, K. H., Gallusser, F., Koehler, J., Remy, N. and Scott, S. L.
#' (2015). Inferring causal impact using Bayesian structural time-series
#' models. *Annals of Applied Statistics*, 9(1), 247--274.
#' @export
#' @examples
#' ## A series with a genuine post-intervention lift of about +3.
#' set.seed(1)
#' pre <- cumsum(rnorm(60, sd = 0.2))
#' post <- pre[60] + cumsum(rnorm(40, sd = 0.2)) + 3
#' x <- c(pre, post)
#' fit <- its_causal(x, intervention = 61L)
#' fit@avg_effect
its_causal <- function(x, intervention, trend = FALSE, level = 0.95) {
  .check_its_input(x, intervention, level)
  intervention <- as.integer(intervention)
  pre <- x[seq_len(intervention - 1L)]
  post <- x[intervention:length(x)]
  n_post <- length(post)

  # Fit the pre-period state-space model and filter it ------------------------

  model <- .its_pre_model(pre, trend)
  fit <- kalman_filter(model, pre)
  last <- nrow(fit@filtered_mean)
  x_state <- fit@filtered_mean[last, ]
  p_state <- fit@filtered_cov[[last]]

  # Project the counterfactual forward across the post-period -----------------

  a <- model@transition[[1L]]
  b <- model@observation[[1L]]
  q <- model@state_cov[[1L]]
  r <- model@obs_cov[[1L]]

  # Forecast the counterfactual, tracking the joint covariance of the state and
  # its running sum so the cumulative- and average-effect intervals reflect the
  # correlation of the forecast errors across the post-period, not an
  # independence assumption that would understate them for a drifting series.

  # The augmented state (x_t, c_t) stacks the latent state with its running sum
  # c_t = c_{t-1} + x_t; its covariance recursion gives var(sum_t B x_t) in the
  # lower-right block, capturing the across-time error correlation exactly. The
  # same process-noise draw w_t enters both blocks, so all four blocks of the
  # augmented noise covariance equal Q -- a zero off-diagonal block would drop
  # the state/running-sum noise correlation and understate the interval.
  m <- model@state_dim
  cf <- numeric(n_post)
  cf_var <- numeric(n_post)
  # Seed the augmented covariance with the filtered state uncertainty in the
  # state block; the running-sum block starts at zero.
  joint_cov <- matrix(0, nrow = 2L * m, ncol = 2L * m)
  joint_cov[seq_len(m), seq_len(m)] <- p_state
  obs_var_sum <- 0
  big_a <- rbind(cbind(a, matrix(0, m, m)), cbind(a, diag(m)))
  big_q <- rbind(cbind(q, q), cbind(q, q))

  for (t in seq_len(n_post)) {
    x_state <- as.numeric(a %*% x_state)
    p_state <- a %*% p_state %*% t(a) + q
    cf[t] <- as.numeric(b %*% x_state)
    cf_var[t] <- as.numeric(b %*% p_state %*% t(b) + r)

    joint_cov <- big_a %*% joint_cov %*% t(big_a) + big_q
    obs_var_sum <- obs_var_sum + as.numeric(r)
  }

  # Contrast and propagate the interval --------------------------------------

  z <- stats::qnorm(1 - (1 - level) / 2)
  cf_sd <- sqrt(cf_var)
  pointwise <- post - cf
  cumulative <- cumsum(pointwise)

  # Variance of the summed counterfactual: the running-sum block of the joint
  # state covariance, mapped through the observation matrix, plus the summed
  # (independent) observation noise.
  b_sum <- cbind(matrix(0, nrow(b), m), b)
  cum_state_var <- as.numeric(b_sum %*% joint_cov %*% t(b_sum)) + obs_var_sum
  avg_effect <- mean(pointwise)
  avg_sd <- sqrt(cum_state_var) / n_post

  its_fit(
    time = as.numeric(intervention:length(x)),
    observed = post,
    counterfactual = cf,
    counterfactual_lower = cf - z * cf_sd,
    counterfactual_upper = cf + z * cf_sd,
    pointwise_effect = pointwise,
    cumulative_effect = cumulative,
    avg_effect = avg_effect,
    avg_effect_lower = avg_effect - z * avg_sd,
    avg_effect_upper = avg_effect + z * avg_sd,
    level = level
  )
}

#' Fit the pre-period state-space model by maximum likelihood
#'
#' Constructs a local-level or local-linear-trend [ssm] for the
#' interrupted-time-series counterfactual and calibrates its noise scale by
#' maximising the Kalman log-likelihood over the pre-period. A single parameter
#' -- the log signal-to-noise ratio, the share of variance the state carries
#' relative to the observation -- is optimised by [stats::optimize]; the
#' overall scale is then matched to the pre-period variance. Estimating this
#' split rather than fixing it heuristically is what makes the forecast
#' interval well-calibrated: a near-random-walk series is allowed to put almost
#' all its variance in the state, widening the counterfactual band correctly.
#'
#' @param pre The pre-intervention numeric series.
#' @param trend Logical scalar — local-linear-trend if `TRUE`, else local
#'   level.
#'
#' @returns An [ssm] with the maximum-likelihood noise scale.
#' @noRd
#' @keywords internal
.its_pre_model <- function(pre, trend) {
  v <- stats::var(diff(pre))
  if (!is.finite(v) || v <= 0) {
    v <- stats::var(pre)
  }
  if (!is.finite(v) || v <= 0) {
    v <- 1
  }

  # Maximise the Kalman log-likelihood over the log signal-to-noise ratio.
  # ratio = state-noise / observation-noise; everything else scales with v.
  neg_log_lik <- function(log_ratio) {
    ratio <- exp(log_ratio)
    model <- .its_build_model(pre, trend, v, ratio)
    -kalman_filter(model, pre)@log_lik
  }
  opt <- tryCatch(
    stats::optimize(neg_log_lik, interval = c(-8, 8)),
    error = function(e) NULL
  )
  ratio <- if (is.null(opt)) 0.5 else exp(opt$minimum)

  .its_build_model(pre, trend, v, ratio)
}

#' Assemble a local-level or local-linear-trend model at a given noise ratio
#'
#' Builds the [ssm] the interrupted-time-series counterfactual uses for one
#' candidate signal-to-noise `ratio` (state-noise relative to observation
#' noise) and overall scale `v`. Factored out so the likelihood optimisation in
#' `.its_pre_model()` and the final model build share one construction.
#'
#' @param pre The pre-intervention numeric series.
#' @param trend Logical scalar — local-linear-trend if `TRUE`, else local level.
#' @param v Numeric scalar — the overall variance scale.
#' @param ratio Numeric scalar — the state-to-observation noise ratio.
#'
#' @returns An [ssm].
#' @noRd
#' @keywords internal
.its_build_model <- function(pre, trend, v, ratio) {
  obs_var <- v / (1 + ratio)
  state_var <- ratio * obs_var
  if (isTRUE(trend)) {
    a <- matrix(c(1, 0, 1, 1), nrow = 2L)
    b <- matrix(c(1, 0), nrow = 1L)
    q <- diag(c(state_var, state_var / 5))
    ssm(
      transition = a, observation = b, state_cov = q, obs_cov = obs_var,
      init_state = c(pre[1L], 0), init_cov = diag(c(v * 10, v))
    )
  } else {
    ssm(
      transition = 1, observation = 1, state_cov = state_var, obs_cov = obs_var,
      init_state = pre[1L], init_cov = v * 10
    )
  }
}

#' Validate the input to its_causal()
#'
#' @param x The full series.
#' @param intervention The intervention index.
#' @param level The interval coverage.
#'
#' @returns `NULL`, invisibly; called for its error side effect.
#' @noRd
#' @keywords internal
.check_its_input <- function(x, intervention, level) {
  if (!is.numeric(x) || !is.null(dim(x)) || anyNA(x)) {
    cli::cli_abort("`x` must be a numeric vector with no missing values.")
  }
  if (length(intervention) != 1L || intervention != as.integer(intervention)) {
    cli::cli_abort("`intervention` must be a single integer index.")
  }
  if (intervention < 11L) {
    cli::cli_abort(c(
      "`intervention` must leave at least ten pre-intervention observations.",
      "i" = "The counterfactual model is fitted on the pre-period."
    ))
  }
  if (intervention > length(x)) {
    cli::cli_abort("`intervention` must not exceed `length(x)`.")
  }
  if (length(level) != 1L || level <= 0 || level >= 1) {
    cli::cli_abort("`level` must be a single number in (0, 1).")
  }
  invisible(NULL)
}

#' @export
S7::method(print, its_fit) <- function(x, ...) {
  pct <- round(100 * x@level)
  cat(sprintf(
    "<its_fit>: single-case interrupted time series (%d post-period steps)\n",
    length(x@observed)
  ))
  cat(sprintf(
    "  average effect : %.4f  [%.4f, %.4f]  (%d%% interval)\n",
    x@avg_effect, x@avg_effect_lower, x@avg_effect_upper, pct
  ))
  cat(sprintf(
    "  cumulative effect : %.4f\n",
    x@cumulative_effect[length(x@cumulative_effect)]
  ))
  spans_zero <- x@avg_effect_lower <= 0 && x@avg_effect_upper >= 0
  cat(sprintf(
    "  interval %s zero\n",
    if (spans_zero) "straddles" else "excludes"
  ))
  invisible(x)
}
