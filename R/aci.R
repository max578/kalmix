# aci.R -- the time-resolved Assimilative Causal Inference (ACI) read-out.
#
# Assimilative Causal Inference (Andreou, Chen and Bollt 2026) casts causality
# as a Bayesian data-assimilation inverse problem: a hidden state is a cause of
# an observed effect at time t exactly when incorporating the future of the
# effect sharpens the estimate of the state, that is when the smoothing
# distribution of the state departs from its filtering distribution. The causal
# information a measurement carries about the hidden state at each step is the
# relative entropy (Kullback-Leibler divergence) between the smoother and the
# filter state distributions (their eq. 7); time-aggregated, the divergence
# profile against an expanding future window yields a threshold-free causal
# information rate and a decision lead-time (their eq. 9).
#
# kalmix is the orchestra's home for the closed-form Gaussian filter/smoother
# pair and this read-out. The divergence is a standard closed form for two
# Gaussians and is implemented natively here in a handful of lines; when kernR
# is installed its `relative_entropy()` is used as an independent oracle in the
# test suite (it is not a runtime dependency). aci() runs the Kalman filter and
# the Rauch-Tung-Striebel smoother of kalman.R and reads the causal information
# series off the two passes; causal_information_rate() forms the objective CIR
# and lead-time from an expanding-window divergence profile.

#' An Assimilative Causal Inference read-out
#'
#' An S7 object holding the output of [aci()]: the per-step causal information
#' series, its time aggregate (the mean causal information), the objective
#' causal information rate and the implied decision lead-time, alongside the
#' filter and smoother passes the read-out was computed from.
#'
#' The causal information at step `t` is the relative entropy of the smoothing
#' distribution from the filtering distribution of the latent state,
#' \deqn{I_t = D_{\mathrm{KL}}\!\left(
#'   N(x_{t\mid n}, P_{t\mid n}) \,\Vert\, N(x_{t\mid t}, P_{t\mid t}) \right),}
#' the information the future of the series adds to the present estimate of the
#' state. It is non-negative and zero only where the future is uninformative
#' about the state at `t`.
#'
#' @param causal_information A numeric vector of length `n` -- the per-step
#'   causal information series \eqn{I_t}, the smoother-from-filter relative
#'   entropy at each step.
#' @param mean_causal_information Numeric scalar -- the time-average of
#'   `causal_information`, a single summary of how much the future informs the
#'   state over the series.
#' @param causal_information_rate Numeric scalar -- the objective,
#'   threshold-free causal information rate (a lead-time in the time units of
#'   the series), or `NA_real_` if it was not requested.
#' @param lead_time Numeric scalar -- the decision lead-time, an alias of the
#'   causal information rate in time units, or `NA_real_` if not requested.
#' @param dt Numeric scalar -- the sampling interval the lead-time is expressed
#'   in. Defaults to `1`, so the lead-time is in steps.
#' @param filter The [kalman_fit] the read-out was computed from.
#' @param smoother The [rts_fit] the read-out was computed from.
#' @param grounding Character scalar -- the grounding token of the causal
#'   verdict, `"grounded"` or `"[unverified]"` (provenance_vocabulary.md §1).
#' @param grounding_reason Character scalar -- a short label for why the verdict
#'   carries that token (`"grounded"`, `"mechanism_unverified"`,
#'   `"model_inadequate"` or `"insufficient_data"`).
#' @param adequacy The [innov_diag] model-adequacy diagnostics the grounding was
#'   read from.
#'
#' @returns An S7 object of class `aci_fit`.
#' @family causal
#' @export
aci_fit <- S7::new_class(
  name = "aci_fit",
  package = "kalmix",
  properties = list(
    causal_information = S7::class_double,
    mean_causal_information = S7::class_double,
    causal_information_rate = S7::new_property(
      class = S7::class_double,
      default = NA_real_
    ),
    lead_time = S7::new_property(
      class = S7::class_double,
      default = NA_real_
    ),
    dt = S7::new_property(class = S7::class_double, default = 1),
    filter = kalman_fit,
    smoother = rts_fit,
    grounding = S7::new_property(
      class = S7::class_character,
      default = "[unverified]"
    ),
    grounding_reason = S7::new_property(
      class = S7::class_character,
      default = NA_character_
    ),
    adequacy = S7::class_any
  ),
  validator = function(self) {
    if (anyNA(self@causal_information) ||
        any(self@causal_information < -1e-8)) {
      return("`causal_information` must be non-negative and free of NA")
    }
    if (length(self@mean_causal_information) != 1L) {
      return("`mean_causal_information` must be a scalar")
    }
    if (length(self@dt) != 1L || !is.finite(self@dt) || self@dt <= 0) {
      return("`dt` must be a finite positive scalar")
    }
    NULL
  }
)

#' Relative entropy between two Gaussians (the ACI causal metric)
#'
#' Computes the Kullback-Leibler divergence
#' \eqn{D_{\mathrm{KL}}(p \Vert q)} between two `k`-variate Gaussians
#' \eqn{p = N(\mu_p, \Sigma_p)} and \eqn{q = N(\mu_q, \Sigma_q)} in closed form,
#' \deqn{D_{\mathrm{KL}}(p \Vert q) = \tfrac12 \left[
#'   \operatorname{tr}(\Sigma_q^{-1}\Sigma_p)
#'   + (\mu_q - \mu_p)^\top \Sigma_q^{-1} (\mu_q - \mu_p) - k
#'   + \log\frac{\lvert\Sigma_q\rvert}{\lvert\Sigma_p\rvert} \right].}
#' This is the operational statistic of Assimilative Causal Inference (Andreou,
#' Chen and Bollt 2026): with `p` the smoother posterior of the latent state and
#' `q` the filter posterior, a non-zero value identifies the future of the
#' observed series as informative about -- a cause of -- the state at that step.
#'
#' The divergence is non-negative, zero exactly when the two Gaussians coincide,
#' and asymmetric in its arguments. This is canonical Gaussian arithmetic, not a
#' reimplementation of a divergence API: the package ships it natively so the
#' read-out is self-contained, and cross-checks the value against
#' `kernR::relative_entropy()` in its test suite when kernR is installed.
#'
#' @param mu_p Numeric vector -- the mean of `p` (the smoother posterior).
#' @param sigma_p A `k` by `k` covariance matrix of `p`.
#' @param mu_q Numeric vector -- the mean of `q` (the filter posterior).
#' @param sigma_q A `k` by `k` covariance matrix of `q`, positive definite (it
#'   is inverted through its Cholesky factor).
#'
#' @returns Numeric scalar -- the relative entropy, floored at zero to absorb
#'   rounding noise.
#' @references
#' Kullback, S. and Leibler, R. A. (1951). On information and sufficiency.
#' *Annals of Mathematical Statistics*, 22(1), 79--86.
#'
#' Andreou, M., Chen, N. and Bollt, E. (2026). Assimilative causal inference.
#' *Nature Communications*, 17, 1854.
#' @noRd
#' @keywords internal
.gaussian_relative_entropy <- function(mu_p, sigma_p, mu_q, sigma_q) {
  mu_p <- as.numeric(mu_p)
  mu_q <- as.numeric(mu_q)
  k <- length(mu_p)
  sigma_p <- as.matrix(sigma_p)
  sigma_q <- as.matrix(sigma_q)

  q_chol <- .safe_chol(sigma_q, "filter covariance")
  q_inv <- chol2inv(q_chol)
  delta <- mu_q - mu_p

  tr_term <- sum(q_inv * sigma_p)
  quad <- sum(delta * as.numeric(q_inv %*% delta))
  log_det_q <- 2 * sum(log(diag(q_chol)))
  log_det_p <- as.numeric(determinant(sigma_p, logarithm = TRUE)$modulus)

  val <- 0.5 * (tr_term + quad - k + (log_det_q - log_det_p))
  max(val, 0)
}

#' Assimilative Causal Inference read-out for a state-space model
#'
#' Runs the closed-form Kalman filter and Rauch-Tung-Striebel smoother of
#' [kalman_filter()] and [rts_smoother()] over an [ssm], then reads the
#' time-resolved causal information series off the two passes following Andreou,
#' Chen and Bollt (2026). The causal information at each step is the relative
#' entropy of the smoothing distribution of the latent state from its filtering
#' distribution,
#' \deqn{I_t = D_{\mathrm{KL}}\!\left(
#'   N(x_{t\mid n}, P_{t\mid n}) \,\Vert\, N(x_{t\mid t}, P_{t\mid t}) \right),}
#' the information the future of the series adds to the present estimate of the
#' state (their eq. 7). The series is summarised by its time-average, the mean
#' causal information.
#'
#' By default the objective causal information rate and the decision lead-time
#' are also computed (`lead_time = TRUE`), through [causal_information_rate()]:
#' the divergence between an expanding-future-window smoother and the complete
#' smoother is integrated into a single threshold-free lead-time in the time
#' units of the series (their eqs. 8--9). Set `lead_time = FALSE` to return the
#' causal information series alone, which is cheaper for long series.
#'
#' The recovered state is interpretable as a cause of the observed series only
#' under the maintained assumption that the [ssm] is the data-generating
#' mechanism. That assumption is not taken on trust: the read-out is grounded
#' through [innovation_diagnostics()], and a model whose standardised
#' innovations are not white, not correctly scaled, or not Gaussian is declared
#' inadequate, so the verdict abstains with a `[unverified]` grounding token.
#' Passing the diagnostics establishes self-consistency only; the verdict is
#' labelled `"grounded"` solely when an adequate model is accompanied by
#' declared, dated `mechanism` provenance, and is the honest `"[unverified]"`
#' otherwise.
#'
#' @param model An [ssm].
#' @param y A numeric vector or `n` by `d` matrix of observations.
#' @param lead_time Logical scalar -- whether to compute the objective causal
#'   information rate and lead-time as well as the causal information series.
#'   Defaults to `TRUE`.
#' @param dt Numeric scalar -- the sampling interval, so the lead-time is in
#'   meaningful time units (for example years for an annual ENSO index).
#'   Defaults to `1`, giving the lead-time in steps.
#' @param eval_points An optional integer vector of step indices at which to
#'   anchor the lead-time computation; passed to [causal_information_rate()].
#'   Defaults to `NULL`, which spreads a handful of anchors across the interior
#'   of the series.
#' @param max_lag Integer scalar -- the largest future-window lag, in steps,
#'   used for the lead-time. Defaults to one fifth of the series length.
#' @param mechanism Optional declared provenance for the model as the
#'   data-generating mechanism: `NULL` (the default, an unverified
#'   mechanism), or a list carrying a `verified_on` `Date` on which the model
#'   was checked against the external authority that owns it. Only a verified
#'   mechanism can raise an adequate read-out to a `"grounded"` verdict.
#' @param alpha Numeric scalar -- the significance level the model-adequacy
#'   diagnostics abstain at. Defaults to `0.05`.
#'
#' @returns An [aci_fit].
#' @family causal
#' @seealso [causal_information_rate()] for the lead-time alone;
#'   [kalman_filter()] and [rts_smoother()] for the underlying passes.
#' @references
#' Andreou, M., Chen, N. and Bollt, E. (2026). Assimilative causal inference.
#' *Nature Communications*, 17, 1854.
#'
#' Kalman, R. E. (1960). A new approach to linear filtering and prediction
#' problems. *Journal of Basic Engineering*, 82(1), 35--45.
#'
#' Rauch, H. E., Tung, F. and Striebel, C. T. (1965). Maximum likelihood
#' estimates of linear dynamic systems. *AIAA Journal*, 3(8), 1445--1450.
#' @export
#' @examples
#' ## A latent oscillation seen through a slow, noisy observation channel: the
#' ## future of the series is informative about the hidden state, so the causal
#' ## information is positive.
#' model <- ssm(
#'   transition = matrix(c(0.6, -0.5, 0.5, 0.6), nrow = 2),
#'   observation = matrix(c(1, 0), nrow = 1),
#'   state_cov = diag(c(0.2, 0.2)),
#'   obs_cov = 1,
#'   init_state = c(0, 0),
#'   init_cov = diag(c(5, 5))
#' )
#' set.seed(1)
#' x <- numeric(200)
#' s <- c(0, 0)
#' for (t in seq_len(200)) {
#'   s <- as.numeric(model@transition[[1]] %*% s) + rnorm(2, sd = 0.45)
#'   x[t] <- s[1] + rnorm(1)
#' }
#' fit <- aci(model, x, lead_time = FALSE)
#' fit@mean_causal_information
aci <- function(model,
                y,
                lead_time = TRUE,
                dt = 1,
                eval_points = NULL,
                max_lag = NULL,
                mechanism = NULL,
                alpha = 0.05) {
  if (!S7::S7_inherits(model, ssm)) {
    cli::cli_abort("`model` must be an {.cls ssm}.")
  }
  if (length(dt) != 1L || !is.finite(dt) || dt <= 0) {
    cli::cli_abort("`dt` must be a finite positive scalar.")
  }

  filter <- kalman_filter(model, y)
  smoother <- rts_smoother(filter)
  info <- .causal_information_series(filter, smoother)

  rate <- NA_real_
  lead <- NA_real_
  if (isTRUE(lead_time)) {
    rate <- causal_information_rate(
      model, y,
      dt = dt, eval_points = eval_points, max_lag = max_lag
    )
    lead <- rate
  }

  adequacy <- innovation_diagnostics(filter, alpha = alpha)
  grounding <- .aci_grounding(adequacy, mechanism)

  aci_fit(
    causal_information = info,
    mean_causal_information = mean(info),
    causal_information_rate = rate,
    lead_time = lead,
    dt = dt,
    filter = filter,
    smoother = smoother,
    grounding = grounding$grounding,
    grounding_reason = grounding$reason,
    adequacy = adequacy
  )
}

#' The per-step causal information series
#'
#' Reads the smoother-from-filter relative entropy off a matched [kalman_fit]
#' and [rts_fit] at every step, the time-resolved causal information of Andreou,
#' Chen and Bollt (2026, eq. 7). Factored out so [aci()] and the lead-time
#' anchors share one definition.
#'
#' @param filter A [kalman_fit].
#' @param smoother The [rts_fit] from the same model and series.
#'
#' @returns A numeric vector of length `n`.
#' @noRd
#' @keywords internal
.causal_information_series <- function(filter, smoother) {
  n <- nrow(filter@filtered_mean)
  vapply(
    seq_len(n),
    function(t) {
      .gaussian_relative_entropy(
        smoother@smoothed_mean[t, ], smoother@smoothed_cov[[t]],
        filter@filtered_mean[t, ], filter@filtered_cov[[t]]
      )
    },
    numeric(1L)
  )
}

#' Objective causal information rate and decision lead-time
#'
#' Computes the objective, threshold-free causal information rate of Andreou,
#' Chen and Bollt (2026, eqs. 8--9) for a state-space model: the effective
#' horizon over which the future of the observed series keeps informing the
#' estimate of the latent state, returned as a decision lead-time in the time
#' units of the series.
#'
#' At an anchor step \eqn{t_0} the smoother is recomputed on an expanding future
#' window: the relative entropy \eqn{D(L)} between the smoother that has seen the
#' series only up to \eqn{t_0 + L} and the complete smoother falls from its
#' filter value at \eqn{L = 0} towards zero as the window grows. The subjective
#' range at tolerance \eqn{\varepsilon} is
#' \eqn{\tau_\varepsilon = \inf\{L : D(L) \le \varepsilon\}}; integrating the
#' tolerance out gives the threshold-free range
#' \deqn{\tau(t_0) = \frac1M \int_0^M \tau_\varepsilon\, d\varepsilon
#'   = \frac1M \int_0^{L_{\max}} D(L)\, dL,}
#' with \eqn{M = D(0)} the filter value, the second equality by integration by
#' parts of the decreasing profile. The rate is the average of \eqn{\tau(t_0)}
#' over the anchors, scaled to time units by `dt`.
#'
#' Two engines compute the same divergence profile. The default `"online"`
#' engine grows the future window through a fixed-point smoother recursion: the
#' smoothed estimate of the anchor state is updated in place as each later
#' observation is assimilated, by accumulating the Rauch-Tung-Striebel gains, so
#' the whole profile costs one filter-smoother pass plus `length(eval_points)`
#' light forward recursions and no re-smoothing. It is exact and applies to
#' time-varying models. The `"expanding"` engine instead reruns the filter and
#' smoother on each truncated series -- the direct construction of the published
#' equations, retained as an independent cross-check and restricted to
#' time-invariant models. The two agree to numerical precision.
#'
#' @param model An [ssm].
#' @param y A numeric vector or matrix of observations.
#' @param dt Numeric scalar -- the sampling interval; the returned lead-time is
#'   in these time units. Defaults to `1`.
#' @param eval_points An optional integer vector of anchor step indices.
#'   Defaults to `NULL`, which places anchors across the interior of the series
#'   (avoiding the ends, where the future window is short).
#' @param max_lag Integer scalar -- the largest future-window lag in steps.
#'   Defaults to one fifth of the series length, capped so an anchor still has
#'   future to integrate over.
#' @param n_lag Integer scalar -- the number of window lengths the divergence
#'   profile is evaluated at, including zero. Defaults to `12`.
#' @param engine Character -- `"online"` (default) computes the lead-time with a
#'   fixed-point smoother recursion that avoids re-smoothing and handles
#'   time-varying models; `"expanding"` reruns the smoother on each truncated
#'   series and requires a time-invariant model. The two agree to numerical
#'   precision and cross-check each other.
#'
#' @returns Numeric scalar -- the objective causal information rate (a lead-time
#'   in the time units set by `dt`), `0` where no future information is
#'   recoverable.
#' @family causal
#' @seealso [aci()].
#' @references
#' Andreou, M., Chen, N. and Bollt, E. (2026). Assimilative causal inference.
#' *Nature Communications*, 17, 1854.
#' @export
#' @examples
#' model <- ssm(
#'   transition = matrix(c(0.6, -0.5, 0.5, 0.6), nrow = 2),
#'   observation = matrix(c(1, 0), nrow = 1),
#'   state_cov = diag(c(0.2, 0.2)),
#'   obs_cov = 1,
#'   init_state = c(0, 0),
#'   init_cov = diag(c(5, 5))
#' )
#' set.seed(1)
#' x <- numeric(200)
#' s <- c(0, 0)
#' for (t in seq_len(200)) {
#'   s <- as.numeric(model@transition[[1]] %*% s) + rnorm(2, sd = 0.45)
#'   x[t] <- s[1] + rnorm(1)
#' }
#' causal_information_rate(model, x)
causal_information_rate <- function(model,
                                    y,
                                    dt = 1,
                                    eval_points = NULL,
                                    max_lag = NULL,
                                    n_lag = 12L,
                                    engine = c("online", "expanding")) {
  engine <- match.arg(engine)
  if (!S7::S7_inherits(model, ssm)) {
    cli::cli_abort("`model` must be an {.cls ssm}.")
  }
  if (length(dt) != 1L || !is.finite(dt) || dt <= 0) {
    cli::cli_abort("`dt` must be a finite positive scalar.")
  }
  y <- .check_observations(y, model)
  n <- nrow(y)

  # The expanding engine reruns the filter and smoother on truncated series,
  # which a time-varying model cannot honour (its per-step matrices would no
  # longer match the truncated length). The online engine assimilates forward in
  # place and has no such restriction, so only the expanding engine needs a
  # static model.
  if (identical(engine, "expanding")) {
    static <- all(vapply(
      c("transition", "observation", "state_cov", "obs_cov"),
      function(nm) length(S7::prop(model, nm)) == 1L,
      logical(1L)
    ))
    if (!static) {
      cli::cli_abort(c(
        "The expanding engine needs a time-invariant {.cls ssm}.",
        "i" = "Use {.code engine = \"online\"} for a time-varying model."
      ))
    }
  }

  if (is.null(max_lag)) {
    max_lag <- max(5L, as.integer(n / 5L))
  }
  max_lag <- as.integer(max_lag)
  if (max_lag < 1L) {
    cli::cli_abort("`max_lag` must be at least one step.")
  }

  eval_points <- .aci_eval_points(eval_points, n, max_lag)
  lags <- unique(as.integer(round(seq(0, max_lag, length.out = n_lag))))

  # One filter-smoother pass; the complete smoother is the reference every
  # window is compared with, and the filter feeds the online recursion.
  filter <- kalman_filter(model, y)
  complete <- rts_smoother(filter)

  per_anchor <- vapply(
    eval_points,
    function(t0) {
      if (identical(engine, "online")) {
        .aci_objective_online(filter, complete, t0, lags)
      } else {
        .aci_objective_at(model, y, t0, lags, complete)
      }
    },
    numeric(1L)
  )
  mean(per_anchor) * dt
}

#' Objective causal information range at one anchor step
#'
#' Forms the expanding-future-window divergence profile \eqn{D(L)} at anchor
#' `t0` and reduces it to the threshold-free range
#' \eqn{(1/M)\int_0^{L_{\max}} D(L)\, dL} by the trapezoidal rule. \eqn{D(0)} is
#' the filter-from-complete-smoother relative entropy; \eqn{D(L)} for `L > 0`
#' reruns the smoother on the series truncated at `t0 + L`. The range is
#' returned in steps; the caller scales by `dt`.
#'
#' @param model An [ssm].
#' @param y The `n` by `d` observation matrix.
#' @param t0 Integer scalar -- the anchor step.
#' @param lags An increasing integer vector of window lengths starting at zero.
#' @param complete The [rts_fit] over the whole series.
#'
#' @returns Numeric scalar -- the objective range at `t0`, in steps.
#' @noRd
#' @keywords internal
.aci_objective_at <- function(model, y, t0, lags, complete) {
  n <- nrow(y)
  mc <- complete@smoothed_mean[t0, ]
  pc <- complete@smoothed_cov[[t0]]

  divergence <- vapply(
    lags,
    function(lag) {
      end <- min(t0 + lag, n)
      window <- if (lag == 0L || end <= t0) {
        # No future incorporated: the filter posterior at t0 is the lag-0 state.
        filter <- kalman_filter(model, y[seq_len(t0), , drop = FALSE])
        list(
          mean = filter@filtered_mean[t0, ],
          cov = filter@filtered_cov[[t0]]
        )
      } else {
        smoothed <- rts_smoother(
          kalman_filter(model, y[seq_len(end), , drop = FALSE])
        )
        list(
          mean = smoothed@smoothed_mean[t0, ],
          cov = smoothed@smoothed_cov[[t0]]
        )
      }
      .gaussian_relative_entropy(window$mean, window$cov, mc, pc)
    },
    numeric(1L)
  )

  .cir_trapezoid(as.numeric(lags), divergence)
}

#' Objective causal information range at one anchor by the online smoother
#'
#' The fixed-point-smoother equivalent of [.aci_objective_at()]: forms the same
#' expanding-future-window divergence profile \eqn{D(L)} at anchor `t0` without
#' re-smoothing. The smoothed estimate of the anchor state under data to
#' \eqn{t_0 + L} is grown forward from the filtered estimate by accumulating
#' the Rauch-Tung-Striebel gains
#' \eqn{C_t = P_{t|t} A_{t+1}^\top P_{t+1|t}^{-1}} into a product \eqn{\Phi},
#' telescoping the smoother recursion so that
#' \eqn{x_{t_0 \mid t_0 + L} = x_{t_0 \mid t_0} +
#'   \sum_{t=t_0+1}^{t_0+L} \Phi_{t_0,t}\,(x_{t \mid t} - x_{t \mid t-1})}
#' and the covariance analogously. This is exact and honours time-varying
#' models. The profile is recorded at the requested `lags` and reduced by the
#' same trapezoidal rule.
#'
#' @param filter The [kalman_fit] over the whole series.
#' @param complete The [rts_fit] over the whole series.
#' @param t0 Integer scalar -- the anchor step.
#' @param lags An increasing integer vector of window lengths starting at zero.
#'
#' @returns Numeric scalar -- the objective range at `t0`, in steps.
#' @noRd
#' @keywords internal
.aci_objective_online <- function(filter, complete, t0, lags) {
  model <- filter@model
  n <- nrow(filter@filtered_mean)
  m <- ncol(filter@filtered_mean)
  mc <- complete@smoothed_mean[t0, ]
  pc <- complete@smoothed_cov[[t0]]

  want <- as.integer(lags)
  max_lag <- max(want)
  divergence <- rep(NA_real_, length(want))

  # L = 0: the filter posterior at the anchor against the complete smoother.
  phi <- diag(1, m)
  x_jm <- filter@filtered_mean[t0, ]
  p_jm <- filter@filtered_cov[[t0]]
  divergence[want == 0L] <- .gaussian_relative_entropy(x_jm, p_jm, mc, pc)

  for (lag in seq_len(max_lag)) {
    t <- t0 + lag
    if (t > n) {
      break
    }
    a_t <- .at_step(model@transition, t)
    c_prev <- filter@filtered_cov[[t - 1L]] %*% t(a_t) %*% chol2inv(
      .safe_chol(filter@predicted_cov[[t]], "predicted covariance")
    )
    phi <- phi %*% c_prev
    x_jm <- x_jm + as.numeric(
      phi %*% (filter@filtered_mean[t, ] - filter@predicted_mean[t, ])
    )
    p_jm <- .symmetrise(
      p_jm + phi %*%
        (filter@filtered_cov[[t]] - filter@predicted_cov[[t]]) %*% t(phi)
    )
    if (any(want == lag)) {
      divergence[want == lag] <- .gaussian_relative_entropy(x_jm, p_jm, mc, pc)
    }
  }

  # If the series ended before max_lag, the window cannot extend further; carry
  # the last computed divergence forward so the profile stays well-formed.
  if (anyNA(divergence)) {
    last <- max(which(!is.na(divergence)))
    divergence[is.na(divergence)] <- divergence[last]
  }

  .cir_trapezoid(as.numeric(want), divergence)
}

#' Threshold-free range from a divergence-versus-lag profile
#'
#' Reduces a decreasing divergence profile \eqn{D(L)} to the objective causal
#' influence range \eqn{(1/M)\int_0^{L_{\max}} D(L)\, dL} by the trapezoidal
#' rule, with \eqn{M = D(0)} the normalising filter value (Andreou, Chen and
#' Bollt 2026, eq. 9). Returns zero when there is no recoverable future
#' information: a near-zero normaliser would otherwise turn a flat,
#' numerically-noisy profile into a spurious large ratio, so `M` is floored
#' against a small absolute tolerance below which the anchor is treated as
#' carrying no causal information.
#'
#' @param lag An increasing numeric vector of lags, including zero.
#' @param divergence The matching divergence values; `divergence[1]` is `M`.
#'
#' @returns Numeric scalar -- the range, in the units of `lag`.
#' @noRd
#' @keywords internal
.cir_trapezoid <- function(lag, divergence) {
  m <- divergence[1L]
  if (!is.finite(m) || m <= 1e-6) {
    return(0)
  }
  k <- length(lag)
  area <- sum(diff(lag) * (divergence[-k] + divergence[-1L]) / 2)
  max(area / m, 0)
}

#' Default anchor steps for the lead-time computation
#'
#' Validates a user-supplied set of anchor steps, or places a default spread of
#' anchors across the interior of the series so each has a usable future window
#' for the expanding-window divergence profile.
#'
#' @param eval_points The user value (an integer vector or `NULL`).
#' @param n Integer scalar -- the series length.
#' @param max_lag Integer scalar -- the largest future-window lag.
#'
#' @returns An integer vector of anchor steps in `seq_len(n)`.
#' @noRd
#' @keywords internal
.aci_eval_points <- function(eval_points, n, max_lag) {
  if (!is.null(eval_points)) {
    eval_points <- as.integer(eval_points)
    if (anyNA(eval_points) || any(eval_points < 1L) ||
        any(eval_points > n)) {
      cli::cli_abort("`eval_points` must be step indices in `seq_len(n)`.")
    }
    return(unique(eval_points))
  }
  # Spread anchors over the interior, leaving room for the future window and a
  # short lead-in for the filter to settle.
  lo <- max(2L, as.integer(0.1 * n))
  hi <- max(lo, n - max_lag)
  if (hi <= lo) {
    return(as.integer(round(seq(2L, n - 1L, length.out = min(6L, n - 2L)))))
  }
  unique(as.integer(round(seq(lo, hi, length.out = 6L))))
}

#' @export
S7::method(print, aci_fit) <- function(x, ...) {
  n <- length(x@causal_information)
  cat(sprintf(
    "<aci_fit>: assimilative causal inference (%d step%s)\n",
    n, if (n == 1L) "" else "s"
  ))
  cat(sprintf(
    "  mean causal information : %.4f nats\n",
    x@mean_causal_information
  ))
  if (is.finite(x@lead_time)) {
    unit <- if (isTRUE(all.equal(x@dt, 1))) "steps" else "time units"
    cat(sprintf(
      "  decision lead-time      : %.3f %s\n",
      x@lead_time, unit
    ))
  } else {
    cat("  decision lead-time      : not computed\n")
  }
  reason <- if (is.na(x@grounding_reason)) {
    ""
  } else {
    sprintf(" (%s)", x@grounding_reason)
  }
  cat(sprintf(
    "  grounding               : %s%s\n", x@grounding, reason
  ))
  invisible(x)
}
