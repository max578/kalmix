# filter.R -- the pairs-trading regime filter (kalmix_filter).
#
# kalmix_filter() is the pairs-trading application of the general
# regime-switching machinery: it turns a fitted mean-reverting spread into a
# set of linear-Gaussian regimes -- one per volatility state of the spread --
# and filters the observed spread through the native mixture_filter(). Each
# regime is the discrete-time state-space form of the Ornstein-Uhlenbeck spread
# of Elliott, van der Hoek and Malcolm (2005); the regimes differ in their
# observation noise, so the filter tracks whether the spread is in a calm or a
# turbulent state -- the regime read-out a mean-reversion trade conditions on.
#
# The engine is the package's own native filter and never requires a sibling
# package. An optional proxymix path is honoured if that package and its
# mixture-Kalman primitive are installed, but it is never on the critical path.

#' Filter a regime-switching spread state (mixture Kalman)
#'
#' Filters the latent mean-reverting state of a spread through a finite
#' Gaussian-mixture Kalman recursion, the regime-switching generalisation of
#' the scalar filter of Elliott, van der Hoek and Malcolm (2005). The fitted
#' [spread_model] is expanded into `n_regimes` linear-Gaussian regimes that
#' share the spread's mean-reversion dynamics but differ in observation noise,
#' from calm to turbulent; the observed spread is filtered through
#' [mixture_filter()], and the returned [regime_fit] carries both the
#' filtered spread state and the posterior probability of each volatility
#' regime at every step.
#'
#' The recursion runs on the package's native engine and needs no sibling
#' package. If `proxymix` and its mixture-Kalman primitive are installed and
#' `engine = "proxymix"` is requested, that engine is used instead; otherwise
#' the native filter is used and remains the default.
#'
#' The regime models operate on deviations from the fitted long-run mean: the
#' discrete-time Ornstein-Uhlenbeck dynamics are the plain autoregressive
#' contraction of \eqn{x_t - \theta}, so the observed spread is de-meaned by
#' `model@theta` before filtering and the returned state estimates are
#' reported back on the observed scale. A spread whose long-run mean is far
#' from zero is therefore tracked around \eqn{\theta}, not pulled toward zero.
#'
#' @param model A [spread_model] giving the mean-reversion dynamics of the
#'   spread, as returned by [ou_fit()].
#' @param y A numeric vector — the observed spread series to filter.
#' @param n_regimes Integer scalar — the number of volatility regimes to track.
#'   Defaults to `2`.
#' @param vol_multipliers A numeric vector of length `n_regimes` giving each
#'   regime's observation-noise multiple of the fitted spread variance, or
#'   `NULL` for a geometric spread from calm to turbulent. Defaults to `NULL`.
#' @param engine Character scalar — `"native"` (default) or `"proxymix"`.
#'
#' @returns A [regime_fit] over the spread's volatility regimes. When
#'   `n_regimes = 1` the regime machinery is skipped and the [kalman_fit] of a
#'   plain Kalman pass is returned instead.
#' @family filter
#' @seealso [mixture_filter()] for the general engine; [ou_fit()] for the model.
#' @references
#' Elliott, R. J., van der Hoek, J. and Malcolm, W. P. (2005). Pairs trading.
#' *Quantitative Finance*, 5(3), 271--276.
#' @export
#' @examples
#' x <- as.numeric(stats::arima.sim(list(ar = 0.85), n = 400))
#' model <- ou_fit(x)
#' fit <- kalmix_filter(model, x)
#' tail(round(fit@regime_prob, 2))
kalmix_filter <- function(model,
                          y,
                          n_regimes = 2L,
                          vol_multipliers = NULL,
                          engine = c("native", "proxymix")) {
  engine <- match.arg(engine)
  if (!S7::S7_inherits(model, spread_model)) {
    cli::cli_abort("`model` must be a {.cls spread_model}.")
  }
  if (!is.numeric(y) || length(y) < 2L) {
    cli::cli_abort("`y` must be a numeric vector of at least two observations.")
  }
  n_regimes <- as.integer(n_regimes)
  if (length(n_regimes) != 1L || n_regimes < 1L) {
    cli::cli_abort("`n_regimes` must be a positive integer scalar.")
  }

  multipliers <- .spread_vol_multipliers(vol_multipliers, n_regimes)
  models <- .spread_regime_models(model, multipliers)

  # The regimes are the OU deviation form, so the filter runs on the de-meaned
  # spread and the state estimates are shifted back to the observed scale.
  # Innovations, covariances and the log-likelihood are invariant to the shift.
  y_dev <- y - model@theta

  if (n_regimes == 1L) {
    fit <- kalman_filter(models[[1L]], y_dev)
    fit@filtered_mean <- fit@filtered_mean + model@theta
    fit@predicted_mean <- fit@predicted_mean + model@theta
    return(fit)
  }
  fit <- mixture_filter(models, y_dev, engine = engine)
  fit@state_mean <- fit@state_mean + model@theta
  fit
}

#' Per-regime observation-noise multipliers for the spread filter
#'
#' Returns the validated vector of observation-noise multiples, one per regime.
#' A `NULL` argument yields a geometric ladder from `0.5` (calm) upward, so two
#' regimes contrast a calm and a turbulent spread and more regimes interpolate.
#'
#' @param vol_multipliers The user argument, or `NULL`.
#' @param n_regimes Integer scalar — the number of regimes.
#'
#' @returns A numeric vector of length `n_regimes`.
#' @noRd
#' @keywords internal
.spread_vol_multipliers <- function(vol_multipliers, n_regimes) {
  if (is.null(vol_multipliers)) {
    return(0.5 * 2^(seq_len(n_regimes) - 1L))
  }
  vol_multipliers <- as.numeric(vol_multipliers)
  if (length(vol_multipliers) != n_regimes || any(vol_multipliers <= 0)) {
    cli::cli_abort(
      "`vol_multipliers` must be {n_regimes} positive value{?s}."
    )
  }
  vol_multipliers
}

#' Build the per-regime state-space models for a spread
#'
#' Expands a fitted [spread_model] into one linear-Gaussian [ssm] per volatility
#' regime, built on the spread's deviation from its fitted long-run mean. Every
#' regime shares the discrete-time mean-reversion dynamics implied by the
#' Ornstein-Uhlenbeck fit -- the autoregressive contraction \eqn{e^{-\kappa
#' \Delta t}} of the deviation toward zero -- and differs only in its
#' observation-noise scale, set by the regime's volatility multiplier. The
#' caller is responsible for de-meaning the observations to match.
#'
#' @param model A [spread_model].
#' @param multipliers A numeric vector of per-regime observation-noise
#'   multipliers.
#'
#' @returns A list of [ssm] objects, one per regime.
#' @noRd
#' @keywords internal
.spread_regime_models <- function(model, multipliers) {
  phi <- exp(-model@kappa * model@dt)
  # Stationary one-step innovation variance of the discrete OU process.
  innov_var <- model@sigma^2 * (1 - phi^2) / (2 * model@kappa)
  base_obs <- max(innov_var, .Machine$double.eps)

  lapply(multipliers, function(mult) {
    # Deviation form: the state is x_t - theta, so the transition is the plain
    # contraction phi with no intercept term and the prior sits at zero with
    # the stationary OU variance. An intercept-free transition applied to the
    # raw spread would revert the state to zero, not theta.
    ssm(
      transition = phi,
      observation = 1,
      state_cov = innov_var,
      obs_cov = base_obs * mult,
      init_state = 0,
      init_cov = model@sigma^2 / (2 * model@kappa)
    )
  })
}

#' Whether the proxymix mixture-Kalman primitive is available
#'
#' @returns Logical scalar.
#' @noRd
#' @keywords internal
.proxymix_has_filter <- function() {
  requireNamespace("proxymix", quietly = TRUE) &&
    exists("gmm_filter", envir = asNamespace("proxymix"))
}
