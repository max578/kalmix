## Standalone spread-construction and mean-reversion estimation.
##
## Neither function depends on a sibling package: kalmix builds and fits a
## cointegration spread on its own. This keeps kalmix standalone-functional.

#' Build a cointegration spread from two price series
#'
#' Forms the spread `a - alpha - beta * b` whose stationarity underlies a
#' pairs trade, where `beta` is the ordinary-least-squares hedge ratio of
#' `a` on `b` and `alpha` the intercept. The hedge ratio and intercept are
#' attached to the result as attributes.
#'
#' @param a Numeric vector — the first asset's price (or log-price) series.
#' @param b Numeric vector — the second asset's series, the same length as
#'   `a`.
#' @param intercept Logical scalar — whether to include an intercept in the
#'   hedge regression. Defaults to `TRUE`.
#'
#' @returns A numeric vector of the same length as `a`, the spread, carrying
#'   `hedge_ratio` and `intercept` attributes.
#' @family spread
#' @export
#' @examples
#' set.seed(1)
#' b <- cumsum(stats::rnorm(300))
#' a <- 2 * b + stats::rnorm(300)
#' s <- spread_series(a, b)
#' attr(s, "hedge_ratio")
spread_series <- function(a, b, intercept = TRUE) {
  .check_price_pair(a, b)
  fit <- if (isTRUE(intercept)) {
    stats::lm(a ~ b)
  } else {
    stats::lm(a ~ b - 1)
  }
  coefs <- stats::coef(fit)
  beta <- unname(coefs[["b"]])
  alpha <- if (isTRUE(intercept)) unname(coefs[["(Intercept)"]]) else 0
  spread <- as.numeric(stats::residuals(fit))
  attr(spread, "hedge_ratio") <- beta
  attr(spread, "intercept") <- alpha
  spread
}

#' Fit Ornstein-Uhlenbeck mean reversion to a spread
#'
#' Estimates the mean-reverting (Ornstein-Uhlenbeck) dynamics of a spread by
#' fitting the exact discrete-time autoregression
#' \deqn{X_t = a + b\, X_{t-1} + \varepsilon_t,}
#' then mapping the autoregressive coefficients back to the continuous-time
#' parameters: \eqn{\kappa = -\log(b)/\Delta t}, \eqn{\theta = a/(1-b)}, and
#' \eqn{\sigma^2 = 2\kappa\, \mathrm{Var}(\varepsilon) / (1 - b^2)}.
#'
#' A spread that is not mean-reverting has an autoregressive coefficient
#' outside `(0, 1)`; the function refuses such a series rather than report a
#' meaningless negative half-life.
#'
#' @param x Numeric vector — the spread series (for example from
#'   [spread_series()]).
#' @param dt Numeric scalar — the sampling interval, in the time units the
#'   reported `kappa` should use. Defaults to `1`.
#'
#' @returns A [spread_model].
#' @family spread
#' @export
#' @examples
#' x <- as.numeric(stats::arima.sim(list(ar = 0.85), n = 800))
#' ou_fit(x)
ou_fit <- function(x, dt = 1) {
  .check_spread(x, dt)

  # Exact discrete-time AR(1) regression -------------------------------------

  n <- length(x)
  y <- x[-1L]
  x_lag <- x[-n]
  fit <- stats::lm(y ~ x_lag)
  a_hat <- unname(stats::coef(fit)[["(Intercept)"]])
  b_hat <- unname(stats::coef(fit)[["x_lag"]])
  if (b_hat <= 0 || b_hat >= 1) {
    # A typed refusal: `x` is not compatible with an Ornstein-Uhlenbeck fit at
    # all, so no result is returned. The extra classes let cross-member
    # integration code recognise this as a structured decline (see
    # ORCHESTRA_dev/integration/refusal_contract.R) rather than a bare error.
    cli::cli_abort(
      c(
        "`x` is not consistent with an Ornstein-Uhlenbeck (mean-reverting) process.",
        "i" = "The fitted AR(1) coefficient is {round(b_hat, 4)}; it must lie in (0, 1)."
      ),
      class = c("orchestra_refusal", "kalmix_refusal")
    )
  }

  # Map back to continuous-time Ornstein-Uhlenbeck parameters -----------------

  resid_var <- stats::var(stats::residuals(fit))
  kappa <- -log(b_hat) / dt
  theta <- a_hat / (1 - b_hat)
  sigma <- sqrt(resid_var * 2 * kappa / (1 - b_hat^2))

  spread_model(
    kappa = kappa,
    theta = theta,
    sigma = sigma,
    dt = dt,
    n_obs = as.integer(n),
    hedge_ratio = attr(x, "hedge_ratio") %||% NA_real_,
    call = match.call()
  )
}

#' Validate a pair of price series
#'
#' @param a,b The two series.
#'
#' @returns `NULL`, invisibly; called for its error side effect.
#' @noRd
#' @keywords internal
.check_price_pair <- function(a, b) {
  if (!is.numeric(a) || !is.numeric(b)) {
    cli::cli_abort("`a` and `b` must both be numeric vectors.")
  }
  if (length(a) != length(b)) {
    cli::cli_abort("`a` and `b` must have the same length.")
  }
  if (length(a) < 3L) {
    cli::cli_abort("at least three observations are needed.")
  }
  if (anyNA(a) || anyNA(b)) {
    cli::cli_abort("`a` and `b` must not contain missing values.")
  }
  invisible(NULL)
}

#' Validate a spread series and sampling interval
#'
#' @param x The spread series.
#' @param dt The sampling interval.
#'
#' @returns `NULL`, invisibly; called for its error side effect.
#' @noRd
#' @keywords internal
.check_spread <- function(x, dt) {
  if (!is.numeric(x) || length(x) < 10L) {
    cli::cli_abort("`x` must be a numeric vector of at least ten observations.")
  }
  if (anyNA(x)) {
    cli::cli_abort("`x` must not contain missing values.")
  }
  if (length(dt) != 1L || !is.numeric(dt) || dt <= 0) {
    cli::cli_abort("`dt` must be a positive numeric scalar.")
  }
  invisible(NULL)
}

`%||%` <- function(a, b) if (is.null(a)) b else a
