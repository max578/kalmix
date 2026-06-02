## Dormant mixture-Kalman adapter.
##
## kalmix's headline capability -- filtering a regime-switching spread state
## through a finite Gaussian-mixture Kalman recursion -- is supplied by the
## proxymix mixture-Kalman engine (`proxymix::gmm_filter()`, horizon B). That
## primitive is commissioned but not yet shipped, so this adapter is
## dormant-but-correct (constellation Invariant 6): the signature is stable
## and it errors with a clear commissioning notice. kalmix is proxymix's
## first consumer of `gmm_filter()` and therefore reviews that API from the
## applied side before it sets.

#' Filter a regime-switching spread state (mixture Kalman)
#'
#' Filters the latent mean-reverting state of a spread through a finite
#' Gaussian-mixture Kalman recursion, the regime-switching generalisation of
#' the scalar filter of Elliott, van der Hoek and Malcolm (2005): each
#' mixture component is one cointegration regime, and the recursion is the
#' component-wise Kalman update plus an evidence reweighting.
#'
#' The recursion is supplied by the optional `proxymix` package
#' (`proxymix::gmm_filter()`). Until that primitive ships, this function is
#' dormant: it validates its arguments and then errors with a commissioning
#' notice. The signature is stable, so downstream code and tests can be
#' written against it now.
#'
#' @param model A [spread_model] giving the per-regime mean-reversion
#'   dynamics.
#' @param y Numeric vector — the observed spread series to filter.
#' @param n_regimes Integer scalar — the number of mixture components
#'   (cointegration regimes) to track. Defaults to `2`.
#' @param ... Reserved for forward compatibility with
#'   `proxymix::gmm_filter()`.
#'
#' @returns A filtered-state object (once the proxymix engine is available).
#' @family filter
#' @export
#' @examples
#' x <- as.numeric(stats::arima.sim(list(ar = 0.85), n = 400))
#' model <- ou_fit(x)
#' try(kalmix_filter(model, x))
kalmix_filter <- function(model, y, n_regimes = 2L, ...) {
  if (!S7::S7_inherits(model, spread_model)) {
    cli::cli_abort("`model` must be a {.cls spread_model}.")
  }
  if (!is.numeric(y) || length(y) < 2L) {
    cli::cli_abort("`y` must be a numeric vector of at least two observations.")
  }
  if (length(n_regimes) != 1L || n_regimes < 1L) {
    cli::cli_abort("`n_regimes` must be a positive integer scalar.")
  }
  .require_proxymix_filter()

  # Reached only once the proxymix mixture-Kalman primitive is available. The
  # call is resolved dynamically (rather than via `proxymix::gmm_filter`) so
  # the package checks clean while the commissioned primitive is still
  # pending; `.require_proxymix_filter()` has already guaranteed it exists.
  gmm_filter <- get("gmm_filter", envir = asNamespace("proxymix")) # nocov
  gmm_filter(model = model, y = y, n_regimes = as.integer(n_regimes), ...) # nocov
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

#' Require the proxymix mixture-Kalman primitive, or error clearly
#'
#' @returns `TRUE`, invisibly; errors if the primitive is unavailable.
#' @noRd
#' @keywords internal
.require_proxymix_filter <- function() {
  if (!requireNamespace("proxymix", quietly = TRUE)) {
    cli::cli_abort(c(
      "`kalmix_filter()` needs the {.pkg proxymix} mixture-Kalman engine.",
      "i" = "Install {.pkg proxymix}, then re-run."
    ))
  }
  if (!.proxymix_has_filter()) {
    cli::cli_abort(c(
      "The proxymix mixture-Kalman primitive {.fn proxymix::gmm_filter} is not available yet.",
      "i" = "kalmix has commissioned it (horizon B); it is pending in proxymix.",
      "i" = "See the {.emph kalmix roadmap} vignette."
    ))
  }
  invisible(TRUE)
}
