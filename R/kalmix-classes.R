# kalmix-classes.R -- the mean-reverting-spread application class.
#
# spread_model is the S7 class for the pairs-trading application of the
# package: the Ornstein-Uhlenbeck parameters of a cointegration spread --
# mean-reversion speed, long-run mean, volatility -- plus the half-life and the
# fitting metadata. It feeds kalmix_filter(), which expands a fitted spread into
# regime-switching state-space models for the general mixture_filter() engine.
# The general state-space machinery itself lives in ssm-class.R, kalman.R,
# mixture.R, hmm.R, changepoint.R and its-causal.R.

# The spread_model class ------------------------------------------------------

#' A fitted mean-reverting spread model
#'
#' An S7 object holding the Ornstein-Uhlenbeck (OU) parameters of a
#' cointegration spread, as returned by [ou_fit()]. The continuous-time
#' model is
#' \deqn{dX_t = \kappa (\theta - X_t)\, dt + \sigma\, dW_t,}
#' so `kappa` is the speed of mean reversion, `theta` the long-run mean,
#' and `sigma` the instantaneous volatility. The mean-reversion half-life
#' is \eqn{\log 2 / \kappa}.
#'
#' @param kappa Numeric scalar — the mean-reversion speed (per unit time).
#'   Strictly positive for a mean-reverting spread.
#' @param theta Numeric scalar — the long-run mean the spread reverts to.
#' @param sigma Numeric scalar — the instantaneous volatility. Strictly
#'   positive.
#' @param dt Numeric scalar — the sampling interval the model was fitted at,
#'   in the same time units as `kappa`.
#' @param n_obs Integer scalar — the number of spread observations used.
#' @param hedge_ratio Numeric scalar — the cointegration hedge ratio the
#'   spread was built with, or `NA_real_` if the spread was supplied
#'   directly.
#' @param call The matched call to [ou_fit()].
#'
#' @returns An S7 object of class `spread_model`.
#' @family classes
#' @export
#' @examples
#' x <- as.numeric(stats::arima.sim(list(ar = 0.8), n = 500))
#' ou_fit(x)
spread_model <- S7::new_class(
  name = "spread_model",
  package = "kalmix",
  properties = list(
    kappa = S7::class_double,
    theta = S7::class_double,
    sigma = S7::class_double,
    dt = S7::new_property(class = S7::class_double, default = 1),
    n_obs = S7::new_property(class = S7::class_integer, default = NA_integer_),
    hedge_ratio = S7::new_property(
      class = S7::class_double,
      default = NA_real_
    ),
    call = S7::class_any
  ),
  validator = function(self) {
    if (length(self@kappa) != 1L || !is.finite(self@kappa) ||
        self@kappa <= 0) {
      return("`kappa` must be a finite positive scalar")
    }
    if (length(self@theta) != 1L || !is.finite(self@theta)) {
      return("`theta` must be a finite scalar")
    }
    if (length(self@sigma) != 1L || !is.finite(self@sigma) ||
        self@sigma <= 0) {
      return("`sigma` must be a finite positive scalar")
    }
    NULL
  }
)

# Spread-model summaries ------------------------------------------------------

#' Mean-reversion half-life of a spread model
#'
#' Returns the time for a deviation from the long-run mean to halve in
#' expectation, \eqn{\log 2 / \kappa}, in the time units of the fit.
#'
#' @param model A [spread_model].
#'
#' @returns Numeric scalar.
#' @family classes
#' @export
#' @examples
#' x <- as.numeric(stats::arima.sim(list(ar = 0.8), n = 500))
#' half_life(ou_fit(x))
half_life <- function(model) {
  if (!S7::S7_inherits(model, spread_model)) {
    cli::cli_abort("`model` must be a {.cls spread_model}.")
  }
  log(2) / model@kappa
}

#' @export
S7::method(print, spread_model) <- function(x, ...) {
  cat(sprintf("<spread_model>: Ornstein-Uhlenbeck (n = %s, dt = %s)\n",
              if (is.na(x@n_obs)) "NA" else as.character(x@n_obs),
              format(x@dt)))
  cat(sprintf("  kappa (reversion) : %.5f\n", x@kappa))
  cat(sprintf("  theta (long mean) : %.5f\n", x@theta))
  cat(sprintf("  sigma (vol)       : %.5f\n", x@sigma))
  cat(sprintf("  half-life         : %.3f\n", half_life(x)))
  if (!is.na(x@hedge_ratio)) {
    cat(sprintf("  hedge ratio       : %.5f\n", x@hedge_ratio))
  }
  invisible(x)
}
