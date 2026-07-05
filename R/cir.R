# cir.R -- the objective causal information range (CIR) and lead-time.
#
# causal_information_rate() reduces the expanding-future-window divergence
# profile of Andreou, Chen and Bollt (2026, eqs. 8-9) to a single
# threshold-free decision lead-time. Two engines compute the same profile: the
# default "online" engine grows each anchor's smoothed estimate forward in
# place by accumulating Rauch-Tung-Striebel gains (the discrete linear-Gaussian
# specialisation of the adaptive online smoother of Andreou, Chen and Li,
# arXiv:2411.05870), and the "expanding" engine re-runs the smoother on
# truncated series as an independent cross-check. Split out of aci.R so the
# read-out verb and the range machinery each stay one-screen navigable.

# The rate verb -----------------------------------------------------------

#' Objective causal information rate and decision lead-time
#'
#' Computes the objective, threshold-free causal information rate of Andreou,
#' Chen and Bollt (2026, eqs. 8--9) for a state-space model: the effective
#' horizon over which the future of the observed series keeps informing the
#' estimate of the latent state, returned as a decision lead-time in the time
#' units of the series.
#'
#' At an anchor step \eqn{t_0} the smoother is recomputed on an expanding future
#' window: the relative entropy \eqn{D(L)} of the complete smoother from the
#' smoother that has seen the series only up to \eqn{t_0 + L} (their eq. 8,
#' with the complete smoother as the integrating density) falls from its
#' filter value at \eqn{L = 0} -- the per-step causal information of eq. 7 --
#' towards zero as the window grows. The subjective
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
#'   precision and cross-check each other. A regime-switching model always
#'   uses the expanding construction (the online recursion is
#'   linear-Gaussian only) and `engine` is ignored.
#' @param transition,init_prob The regime-transition matrix and initial
#'   regime distribution when `model` is a list of regimes; passed to
#'   [mixture_filter()] and [mixture_smoother()]. Ignored for an [ssm].
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
                                    engine = c("online", "expanding"),
                                    transition = NULL,
                                    init_prob = NULL) {
  engine <- match.arg(engine)
  is_regime <- is.list(model) && !S7::S7_inherits(model, ssm)
  if (!is_regime && !S7::S7_inherits(model, ssm)) {
    cli::cli_abort(
      "`model` must be an {.cls ssm} or a list of {.cls ssm} regimes."
    )
  }
  if (length(dt) != 1L || !is.finite(dt) || dt <= 0) {
    cli::cli_abort("`dt` must be a finite positive scalar.")
  }
  ref_model <- if (is_regime) model[[1L]] else model
  y <- .check_observations(y, ref_model)
  n <- nrow(y)

  # The expanding engine reruns the filter and smoother on truncated series,
  # which a time-varying model cannot honour (its per-step matrices would no
  # longer match the truncated length). The online engine assimilates forward in
  # place and has no such restriction, so only the expanding engine needs a
  # static model. Regime models always use the expanding construction.
  if (is_regime || identical(engine, "expanding")) {
    check_static <- function(mod) {
      all(vapply(
        c("transition", "observation", "state_cov", "obs_cov",
          "state_intercept", "obs_intercept"),
        function(nm) length(S7::prop(mod, nm)) == 1L,
        logical(1L)
      ))
    }
    static <- if (is_regime) {
      all(vapply(model, check_static, logical(1L)))
    } else {
      check_static(model)
    }
    if (!static) {
      cli::cli_abort(c(
        "The expanding engine needs time-invariant models.",
        "i" = "Use {.code engine = \"online\"} for a time-varying {.cls ssm}."
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

  if (is_regime) {
    filter <- mixture_filter(model, y, transition = transition,
                             init_prob = init_prob)
    complete <- mixture_smoother(model, y, transition = transition,
                                 init_prob = init_prob)
    per_anchor <- vapply(
      eval_points,
      function(t0) {
        .aci_objective_regime(model, y, transition, init_prob, t0, lags,
                              complete, filter)
      },
      numeric(1L)
    )
    return(mean(per_anchor) * dt)
  }

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

#' Objective causal information range at one anchor for a regime model
#'
#' The regime-switching counterpart of `.aci_objective_at()`: the divergence
#' profile compares the complete Kim smoother's collapsed posterior at the
#' anchor with the collapsed estimate under data to `t0 + L` (the GPB1
#' filter at `L = 0`, a truncated Kim smoother beyond). Both densities are
#' the model class's own Gaussian collapses, so the profile inherits that
#' approximation.
#'
#' @param models The list of regime [ssm] objects.
#' @param y The `n` by `d` observation matrix.
#' @param transition,init_prob The regime chain, as given to the verbs.
#' @param t0 Integer scalar -- the anchor step.
#' @param lags An increasing integer vector of window lengths starting at zero.
#' @param complete The [regime_smooth] over the whole series.
#' @param filter The [regime_fit] over the whole series.
#'
#' @returns Numeric scalar -- the objective range at `t0`, in steps.
#' @noRd
#' @keywords internal
.aci_objective_regime <- function(models, y, transition, init_prob,
                                  t0, lags, complete, filter) {
  n <- nrow(y)
  mc <- complete@smoothed_mean[t0, ]
  pc <- complete@smoothed_cov[[t0]]

  divergence <- vapply(
    lags,
    function(lag) {
      end <- min(t0 + lag, n)
      window <- if (lag == 0L || end <= t0) {
        list(
          mean = filter@state_mean[t0, ],
          cov = filter@state_cov[[t0]]
        )
      } else {
        sm <- mixture_smoother(
          models, y[seq_len(end), , drop = FALSE],
          transition = transition, init_prob = init_prob
        )
        list(
          mean = sm@smoothed_mean[t0, ],
          cov = sm@smoothed_cov[[t0]]
        )
      }
      .gaussian_relative_entropy(mc, pc, window$mean, window$cov)
    },
    numeric(1L)
  )

  .cir_trapezoid(as.numeric(lags), divergence)
}

# Per-anchor range engines ------------------------------------------------

#' Objective causal information range at one anchor step
#'
#' Forms the expanding-future-window divergence profile \eqn{D(L)} at anchor
#' `t0` and reduces it to the threshold-free range
#' \eqn{(1/M)\int_0^{L_{\max}} D(L)\, dL} by the trapezoidal rule. \eqn{D(0)} is
#' the complete-smoother-from-filter relative entropy (the per-step causal
#' information of eq. 7); \eqn{D(L)} for `L > 0` reruns the smoother on the
#' series truncated at `t0 + L` and measures the complete smoother from that
#' lagged estimate (eq. 8). The range is returned in steps; the caller scales
#' by `dt`.
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
      .gaussian_relative_entropy(mc, pc, window$mean, window$cov)
    },
    numeric(1L)
  )

  .cir_trapezoid(as.numeric(lags), divergence)
}

#' Objective causal information range at one anchor by the online smoother
#'
#' The fixed-point-smoother equivalent of `.aci_objective_at()`: forms the same
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
  divergence[want == 0L] <- .gaussian_relative_entropy(mc, pc, x_jm, p_jm)

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
      divergence[want == lag] <- .gaussian_relative_entropy(mc, pc, x_jm, p_jm)
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

# Profile reduction and anchor placement -----------------------------------

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

