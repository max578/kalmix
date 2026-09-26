# cir.R -- the objective causal influence range (CIR) and lead-time.
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
#
# The profile builders return the profile itself -- lag, divergence, and which
# lags the record could actually resolve -- and the reduction to a range is one
# function, so the quadrature grid, the choice of functional and the censoring
# contract are all visible in one place rather than buried in three engines.

# The range verb ----------------------------------------------------------

#' Objective causal influence range and decision lead-time
#'
#' Computes the objective, threshold-free causal influence range of Andreou,
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
#' towards zero as the window grows.
#'
#' Two functionals reduce that profile to a range, and they are different
#' functionals rather than two quadratures of one. The subjective range at
#' tolerance \eqn{\varepsilon} is the *last* lag at which the profile still
#' exceeds \eqn{\varepsilon},
#' \eqn{\tau_\varepsilon = \sup\{L : D(L) > \varepsilon\}}, and
#' `functional = "objective_exact"` averages it over the tolerance,
#' \eqn{M^{-1}\int_0^M \tau_\varepsilon \, d\varepsilon} with
#' \eqn{M = \max_L D(L)}. The default `functional = "objective"` is the
#' computationally efficient underestimate the source paper gives,
#' \deqn{\tau(t_0) = \frac1M \int_0^{L_{\max}} D(L)\, dL.}
#' The two coincide when the profile decreases with lag, by the layer-cake
#' identity, and the integral form is strictly the smaller as soon as it does
#' not; a profile that rises before it decays is the ordinary case near the
#' start of a record. Monotonicity is assumed for convenience in the source
#' paper and is not required by the theory, so `"objective"` is reported as a
#' lower bound and `attr(x, "monotone")` records whether the two agree. The
#' range is averaged over the anchors and scaled to time units by `dt`.
#'
#' Two engines compute the same divergence profile. The default `"online"`
#' engine grows the future window through a fixed-point smoother recursion: the
#' smoothed estimate of the anchor state is updated in place as each later
#' observation is assimilated, by accumulating the Rauch-Tung-Striebel gains, so
#' the whole profile costs one filter-smoother pass plus `length(eval_points)`
#' light forward recursions and no re-smoothing. It is exact and applies to
#' time-varying models. The `"expanding"` engine instead reruns the filter and
#' smoother on each truncated series -- the direct construction of the published
#' equations, retained as a cross-check on the smoother recursion and restricted
#' to time-invariant models. The two share the quadrature and the reduction, so
#' their agreement grades the recursion and not the range functional; the
#' independent grading of the functional is the conformance suite against
#' `acir`.
#'
#' An anchor whose future window runs past the end of the record cannot be
#' resolved. The profile is then integrated over the lags the record does
#' support and no further -- the value is a lower bound, not an estimate -- and
#' the anchor is marked censored. So is an anchor whose profile has not decayed
#' to within `margin` of its peak by the last lag it could be evaluated at,
#' since the window was too short to see the influence expire. `attr(x,
#' "censored")` is `TRUE` when any anchor is censored and `attr(x,
#' "censored_fraction")` gives the proportion that are.
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
#' @param n_lag Integer scalar or `NULL` -- the number of window lengths the
#'   divergence profile is evaluated at, including zero. Defaults to `NULL`,
#'   which evaluates every lag in the window, the finest grid the record
#'   supports, capped at `401` points. A thinned grid over-estimates the
#'   integral of a sharply decaying profile, so whenever the grid is thinned
#'   the quadrature is compared with the same reduction on every other point
#'   and a warning is issued when the two disagree by more than `tol` in
#'   relative terms.
#' @param engine Character -- `"online"` (default) computes the lead-time with a
#'   fixed-point smoother recursion that avoids re-smoothing and handles
#'   time-varying models; `"expanding"` reruns the smoother on each truncated
#'   series and requires a time-invariant model. The two agree to numerical
#'   precision and cross-check the smoother recursion. A regime-switching model
#'   always uses the expanding construction (the online recursion is
#'   linear-Gaussian only) and `engine` is ignored.
#' @param functional Character -- `"objective"` (default) for the efficient
#'   integral form, `"objective_exact"` for the threshold-averaged form. See
#'   the details.
#' @param epsilon Numeric vector -- the tolerance grid the exact functional
#'   averages the subjective range over, in nats. Defaults to the 129-point
#'   logarithmic grid `10^seq(-6, 0.5, length.out = 129)`, matching `acir`'s
#'   former default. Ignored by `functional = "objective"`.
#' @param margin Numeric scalar in `(0, 1)` -- the fraction of its peak the
#'   profile must have decayed below by its last evaluated lag for the anchor
#'   to count as resolved. Defaults to `0.1`.
#' @param tol Numeric scalar -- the relative disagreement between the quadrature
#'   and the same quadrature on half the grid points above which a warning is
#'   issued. Defaults to `0.02`. Set to `Inf` to silence the check.
#' @param transition,init_prob The regime-transition matrix and initial
#'   regime distribution when `model` is a list of regimes; passed to
#'   [mixture_filter()] and [mixture_smoother()]. Ignored for an [ssm].
#'
#' @returns Numeric scalar -- the objective causal influence range (a lead-time
#'   in the time units set by `dt`), `0` where no future information is
#'   recoverable, carrying the attributes `censored` (logical),
#'   `censored_fraction` (numeric), `monotone` (logical, whether every anchor's
#'   profile decreases with lag), `n_lag` (the grid size used) and `converged`
#'   (logical, whether a thinned grid agreed with this one to within `tol`;
#'   `NA` when the check does not apply, for example a fully resolved grid).
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
                                   n_lag = NULL,
                                   engine = c("online", "expanding"),
                                   functional = c("objective",
                                                  "objective_exact"),
                                   epsilon = .cir_epsilon_grid(),
                                   margin = 0.1,
                                   tol = 0.02,
                                   transition = NULL,
                                   init_prob = NULL) {
  engine <- match.arg(engine)
  functional <- match.arg(functional)
  is_regime <- is.list(model) && !S7::S7_inherits(model, ssm)
  if (!is_regime && !S7::S7_inherits(model, ssm)) {
    cli::cli_abort(
      "`model` must be an {.cls ssm} or a list of {.cls ssm} regimes."
    )
  }
  if (length(dt) != 1L || !is.finite(dt) || dt <= 0) {
    cli::cli_abort("`dt` must be a finite positive scalar.")
  }
  if (length(margin) != 1L || !is.finite(margin) ||
      margin <= 0 || margin >= 1) {
    cli::cli_abort("`margin` must be a scalar in `(0, 1)`.")
  }
  if (length(tol) != 1L || is.na(tol) || tol <= 0) {
    cli::cli_abort("`tol` must be a positive scalar, or `Inf`.")
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
  lags <- .cir_lag_grid(max_lag, n_lag)

  if (is_regime) {
    filter <- mixture_filter(model, y, transition = transition,
                             init_prob = init_prob)
    complete <- mixture_smoother(model, y, transition = transition,
                                 init_prob = init_prob)
    profiles <- lapply(
      eval_points,
      function(t0) {
        .cir_profile_regime(model, y, transition, init_prob, t0, lags,
                            complete, filter)
      }
    )
  } else {
    # One filter-smoother pass; the complete smoother is the reference every
    # window is compared with, and the filter feeds the online recursion.
    filter <- kalman_filter(model, y)
    complete <- rts_smoother(filter)
    profiles <- lapply(
      eval_points,
      function(t0) {
        if (identical(engine, "online")) {
          .cir_profile_online(filter, complete, t0, lags)
        } else {
          .cir_profile_expanding(model, y, t0, lags, complete)
        }
      }
    )
  }

  .cir_aggregate(profiles, dt, functional, epsilon, margin, tol,
                 length(lags), max_lag)
}

# Per-anchor profile engines ----------------------------------------------

#' The tolerance grid the exact functional averages over
#'
#' The 129-point logarithmic grid from `1e-6` to `10^0.5` nats, matching
#' `acir`'s former `aci_cir()` default (now `aci_range()`). The reference
#' implementation spans the same range with 513 points.
#'
#' @returns A numeric vector.
#' @noRd
#' @keywords internal
.cir_epsilon_grid <- function() {
  10^seq(-6, 0.5, length.out = 129L)
}

#' The lag grid the divergence profile is evaluated at
#'
#' Every integer lag up to `max_lag` by default, that being the finest grid the
#' record supports and therefore the exact quadrature for the profile actually
#' available; a thinned grid over-estimates the integral of a sharply decaying
#' profile. The default caps at 401 points so a very long window cannot make
#' the re-smoothing engines quadratic without the caller asking for it.
#'
#' @param max_lag Integer scalar -- the largest window lag.
#' @param n_lag Integer scalar or `NULL` -- the requested grid size.
#'
#' @returns An increasing integer vector starting at zero.
#' @noRd
#' @keywords internal
.cir_lag_grid <- function(max_lag, n_lag) {
  if (is.null(n_lag)) {
    n_lag <- min(as.integer(max_lag) + 1L, 401L)
  }
  n_lag <- as.integer(n_lag)
  if (length(n_lag) != 1L || is.na(n_lag) || n_lag < 2L) {
    cli::cli_abort("`n_lag` must be a single integer of at least two.")
  }
  unique(as.integer(round(seq(0, max_lag, length.out = n_lag))))
}

#' Divergence profile at one anchor for a regime model
#'
#' The regime-switching counterpart of `.cir_profile_expanding()`: the profile
#' compares the complete Kim smoother's collapsed posterior at the anchor with
#' the collapsed estimate under data to `t0 + L` (the GPB1 filter at `L = 0`, a
#' truncated Kim smoother beyond). Both densities are the model class's own
#' Gaussian collapses, so the profile inherits that approximation.
#'
#' @param models The list of regime [ssm] objects.
#' @param y The `n` by `d` observation matrix.
#' @param transition,init_prob The regime chain, as given to the verbs.
#' @param t0 Integer scalar -- the anchor step.
#' @param lags An increasing integer vector of window lengths starting at zero.
#' @param complete The [regime_smooth] over the whole series.
#' @param filter The [regime_fit] over the whole series.
#'
#' @returns A list with `lag`, `divergence` and `resolved`; see
#'   `.cir_profile_online()`.
#' @noRd
#' @keywords internal
.cir_profile_regime <- function(models, y, transition, init_prob,
                                t0, lags, complete, filter) {
  n <- nrow(y)
  mc <- complete@smoothed_mean[t0, ]
  pc <- complete@smoothed_cov[[t0]]
  resolved <- t0 + lags <= n

  divergence <- vapply(
    seq_along(lags),
    function(i) {
      if (!resolved[i]) {
        return(NA_real_)
      }
      lag <- lags[i]
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

  list(lag = as.numeric(lags), divergence = divergence, resolved = resolved)
}

#' Divergence profile at one anchor by re-smoothing truncated series
#'
#' Forms the expanding-future-window divergence profile \eqn{D(L)} at anchor
#' `t0`. \eqn{D(0)} is the complete-smoother-from-filter relative entropy (the
#' per-step causal information of eq. 7); \eqn{D(L)} for `L > 0` reruns the
#' smoother on the series truncated at `t0 + L` and measures the complete
#' smoother from that lagged estimate (eq. 8). Lags the record cannot reach are
#' returned as `NA` and marked unresolved, never extrapolated.
#'
#' @param model An [ssm].
#' @param y The `n` by `d` observation matrix.
#' @param t0 Integer scalar -- the anchor step.
#' @param lags An increasing integer vector of window lengths starting at zero.
#' @param complete The [rts_fit] over the whole series.
#'
#' @returns A list with `lag`, `divergence` and `resolved`.
#' @noRd
#' @keywords internal
.cir_profile_expanding <- function(model, y, t0, lags, complete) {
  n <- nrow(y)
  mc <- complete@smoothed_mean[t0, ]
  pc <- complete@smoothed_cov[[t0]]
  resolved <- t0 + lags <= n

  divergence <- vapply(
    seq_along(lags),
    function(i) {
      if (!resolved[i]) {
        return(NA_real_)
      }
      lag <- lags[i]
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

  list(lag = as.numeric(lags), divergence = divergence, resolved = resolved)
}

#' Divergence profile at one anchor by the online smoother
#'
#' The fixed-point-smoother equivalent of `.cir_profile_expanding()`: forms the
#' same expanding-future-window divergence profile \eqn{D(L)} at anchor `t0`
#' without re-smoothing. The smoothed estimate of the anchor state under data to
#' \eqn{t_0 + L} is grown forward from the filtered estimate by accumulating
#' the Rauch-Tung-Striebel gains
#' \eqn{C_t = P_{t|t} A_{t+1}^\top P_{t+1|t}^{-1}} into a product \eqn{\Phi},
#' telescoping the smoother recursion so that
#' \eqn{x_{t_0 \mid t_0 + L} = x_{t_0 \mid t_0} +
#'   \sum_{t=t_0+1}^{t_0+L} \Phi_{t_0,t}\,(x_{t \mid t} - x_{t \mid t-1})}
#' and the covariance analogously. This is exact and honours time-varying
#' models. Lags beyond the end of the record are returned as `NA` and marked
#' unresolved.
#'
#' @param filter The [kalman_fit] over the whole series.
#' @param complete The [rts_fit] over the whole series.
#' @param t0 Integer scalar -- the anchor step.
#' @param lags An increasing integer vector of window lengths starting at zero.
#'
#' @returns A list with `lag` (numeric), `divergence` (numeric, `NA` at
#'   unresolved lags) and `resolved` (logical).
#' @noRd
#' @keywords internal
.cir_profile_online <- function(filter, complete, t0, lags) {
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

  list(
    lag = as.numeric(want),
    divergence = divergence,
    resolved = t0 + want <= n
  )
}

# Profile reduction and anchor placement -----------------------------------

#' Reduce a set of per-anchor profiles to one lead-time
#'
#' Reduces every anchor's profile with `.cir_reduce()`, averages the resolved
#' ranges, scales to time units, checks the quadrature against the same
#' reduction on half the grid points, and attaches the censoring and
#' monotonicity record.
#'
#' @param profiles A list of profiles from the per-anchor engines.
#' @param dt Numeric scalar -- the sampling interval.
#' @param functional,epsilon,margin,tol As in [causal_information_rate()].
#' @param n_lag Integer scalar -- the grid size used, recorded on the result.
#' @param max_lag Integer scalar -- the largest window lag. The convergence
#'   check is skipped when the grid already carries every lag up to it, since
#'   no finer grid exists to compare against.
#'
#' @returns A numeric scalar with the documented attributes.
#' @noRd
#' @keywords internal
.cir_aggregate <- function(profiles, dt, functional, epsilon, margin, tol,
                           n_lag, max_lag) {
  reduce_all <- function(thin) {
    lapply(profiles, function(p) {
      keep <- if (thin) seq(1L, length(p$lag), by = 2L) else seq_along(p$lag)
      ok <- keep[p$resolved[keep]]
      .cir_reduce(p$lag[ok], p$divergence[ok], functional = functional,
                  epsilon = epsilon, margin = margin)
    })
  }
  full <- reduce_all(FALSE)
  value <- mean(vapply(full, function(r) r$value, numeric(1L))) * dt

  censored <- vapply(full, function(r) r$censored, logical(1L))
  # A profile the record could not evaluate to the requested lag is censored
  # whatever its shape: the value is an integral over a shorter window.
  truncated <- vapply(profiles, function(p) !all(p$resolved), logical(1L))
  censored <- censored | truncated

  # `converged` is NA when the check does not apply (a resolved grid has
  # nothing finer to compare against, or the check is switched off), and
  # otherwise records whether the thinned-grid quadrature agreed with the
  # full one to within `tol` -- the same test the warning below is raised
  # from, so a caller reading the attribute sees exactly what the warning
  # reported.
  resolved_grid <- n_lag >= as.integer(max_lag) + 1L
  converged <- NA
  if (is.finite(tol) && n_lag > 3L && !resolved_grid) {
    half <- reduce_all(TRUE)
    coarse <- mean(vapply(half, function(r) r$value, numeric(1L))) * dt
    if (is.finite(coarse) && is.finite(value) && abs(value) > 0) {
      shift <- abs(coarse - value) / abs(value)
      converged <- shift <= tol
      if (!converged) {
        cli::cli_warn(c(
          paste(
            "The causal influence range has not converged on this",
            "quadrature grid."
          ),
          "i" = paste0(
            "Dropping every other lag moves it by ",
            sprintf("%.1f%%", 100 * shift), " (`tol` is ",
            sprintf("%.1f%%", 100 * tol), ")."
          ),
          "i" = "Raise {.arg n_lag} or {.arg max_lag}."
        ))
      }
    }
  }

  structure(
    value,
    censored = any(censored),
    censored_fraction = mean(censored),
    monotone = all(vapply(full, function(r) isTRUE(r$monotone), logical(1L))),
    n_lag = n_lag,
    converged = converged
  )
}

#' Threshold-free range from a divergence-versus-lag profile
#'
#' Reduces a divergence profile \eqn{D(L)} to a causal influence range by one of
#' two functionals (Andreou, Chen and Bollt 2026, eqs. 8--9). `"objective"` is
#' the efficient integral form \eqn{M^{-1}\int D(L)\,dL} with
#' \eqn{M = \max_L D(L)} the peak; `"objective_exact"` averages the subjective
#' range over the tolerance grid. Returns zero when there is no recoverable
#' future information: a near-zero peak would otherwise turn a flat,
#' numerically-noisy profile into a spurious large ratio, so `M` is floored
#' against a small absolute tolerance below which the anchor is treated as
#' carrying no causal information.
#'
#' @param lag An increasing numeric vector of lags, including zero.
#' @param divergence The matching divergence values, free of `NA`.
#' @param functional Character -- `"objective"` or `"objective_exact"`.
#' @param epsilon The tolerance grid for the exact functional.
#' @param margin Numeric scalar -- the fraction of its peak the profile must
#'   have fallen below by its last lag to count as resolved.
#'
#' @returns A list with `value`, `peak`, `monotone` and `censored`.
#' @noRd
#' @keywords internal
.cir_reduce <- function(lag, divergence,
                        functional = c("objective", "objective_exact"),
                        epsilon = .cir_epsilon_grid(),
                        margin = 0.1) {
  functional <- match.arg(functional)
  k <- length(lag)
  peak <- if (k == 0L) NA_real_ else max(divergence)
  if (k < 2L || !is.finite(peak) || peak <= 1e-6) {
    return(list(value = 0, peak = peak, monotone = TRUE, censored = FALSE))
  }
  monotone <- all(diff(divergence) <= 1e-10)
  censored <- divergence[k] > margin * peak

  value <- if (identical(functional, "objective")) {
    area <- sum(diff(lag) * (divergence[-k] + divergence[-1L]) / 2)
    max(area / peak, 0)
  } else {
    tau <- .cir_subjective_range(lag, divergence, epsilon)
    inside <- epsilon < peak
    # tau is a decreasing step function of the tolerance; integrate it over
    # (0, peak] with the endpoints the definition fixes -- the whole evaluated
    # window as the tolerance vanishes, zero at the peak itself.
    nodes <- c(0, epsilon[inside], peak)
    heights <- c(max(tau), tau[inside], 0)
    ord <- order(nodes)
    nodes <- nodes[ord]
    heights <- heights[ord]
    j <- length(nodes)
    max(
      sum(diff(nodes) * (heights[-j] + heights[-1L]) / 2) / peak,
      0
    )
  }

  list(value = value, peak = peak, monotone = monotone, censored = censored)
}

#' Objective range from a divergence-versus-lag profile
#'
#' The efficient integral functional alone, kept as a named entry point because
#' it is the quantity the conformance suite grades against
#' `acir::aci_range(method = "l1_linf")$tau`.
#'
#' @param lag An increasing numeric vector of lags, including zero.
#' @param divergence The matching divergence values, free of `NA`.
#'
#' @returns Numeric scalar -- the range, in the units of `lag`.
#' @noRd
#' @keywords internal
.cir_trapezoid <- function(lag, divergence) {
  .cir_reduce(lag, divergence, functional = "objective")$value
}

#' The subjective causal influence range at each tolerance
#'
#' The subjective range at tolerance \eqn{\varepsilon} is the elapsed lag after
#' which the divergence profile stays below \eqn{\varepsilon} -- the *last* lag
#' at which it exceeds the tolerance, advanced to the next grid point (Andreou,
#' Chen and Bollt 2026, eq. 8). A profile that never exceeds the tolerance has
#' range zero; one still above it at the final lag returns that lag, which is a
#' lower bound.
#'
#' @param lag An increasing numeric vector of lags, including zero.
#' @param divergence The matching divergence values.
#' @param epsilon A numeric vector of tolerances, in nats.
#'
#' @returns A numeric vector the length of `epsilon`.
#' @noRd
#' @keywords internal
.cir_subjective_range <- function(lag, divergence, epsilon) {
  k <- length(lag)
  vapply(
    epsilon,
    function(e) {
      above <- which(divergence > e)
      if (!length(above)) {
        return(0)
      }
      lag[min(max(above) + 1L, k)]
    },
    numeric(1L)
  )
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
