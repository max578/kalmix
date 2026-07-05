# mixture-internal.R -- validation for the regime-switching mixture filter.
#
# These helpers check the regime models share a dimension, build a sensible
# default regime-transition matrix and initial distribution when the caller
# omits them, and validate any supplied ones. They keep mixture_filter()
# reading as its contract.

# Regime-model validation -----------------------------------------------------

#' Validate a list of regime models
#'
#' Checks that `models` is a list of at least two [ssm] objects that share a
#' common state and observation dimension, so they can be filtered against the
#' same observation series and collapsed into one Gaussian.
#'
#' @param models The `models` argument to [mixture_filter()].
#'
#' @returns `NULL`, invisibly; called for its error side effect.
#' @noRd
#' @keywords internal
.check_regime_models <- function(models) {
  if (!is.list(models) || length(models) < 2L) {
    cli::cli_abort("`models` must be a list of at least two {.cls ssm} objects.")
  }
  is_ssm <- vapply(models, function(m) S7::S7_inherits(m, ssm), logical(1L))
  if (!all(is_ssm)) {
    cli::cli_abort("every element of `models` must be an {.cls ssm}.")
  }
  state_dims <- vapply(models, function(m) m@state_dim, integer(1L))
  obs_dims <- vapply(models, function(m) m@obs_dim, integer(1L))
  if (length(unique(state_dims)) != 1L) {
    cli::cli_abort("all regime models must share one state dimension.")
  }
  if (length(unique(obs_dims)) != 1L) {
    cli::cli_abort("all regime models must share one observation dimension.")
  }
  invisible(NULL)
}

# Transition and prior defaults -----------------------------------------------

#' Validate or default the regime-transition matrix
#'
#' Returns a validated `k` by `k` row-stochastic transition matrix. A `NULL`
#' argument yields a persistent default that stays in the current regime with
#' probability `0.95` and spreads the rest uniformly over the others.
#'
#' @param transition The `transition` argument, or `NULL`.
#' @param k Integer scalar — the number of regimes.
#'
#' @returns A `k` by `k` numeric matrix with rows summing to one.
#' @noRd
#' @keywords internal
.check_regime_transition <- function(transition, k) {
  if (is.null(transition)) {
    stay <- 0.95
    off <- (1 - stay) / (k - 1L)
    transition <- matrix(off, nrow = k, ncol = k)
    diag(transition) <- stay
    return(transition)
  }
  if (!is.matrix(transition) || nrow(transition) != k ||
      ncol(transition) != k) {
    cli::cli_abort("`transition` must be a {k} by {k} matrix.")
  }
  if (any(transition < 0) || anyNA(transition)) {
    cli::cli_abort("`transition` entries must be non-negative and non-missing.")
  }
  if (!isTRUE(all.equal(unname(rowSums(transition)), rep(1, k)))) {
    cli::cli_abort("each row of `transition` must sum to one.")
  }
  transition
}

#' Validate or default the initial regime distribution
#'
#' Returns a validated length-`k` probability vector. A `NULL` argument yields a
#' uniform prior over the regimes.
#'
#' @param init_prob The `init_prob` argument, or `NULL`.
#' @param k Integer scalar — the number of regimes.
#'
#' @returns A length-`k` numeric vector summing to one.
#' @noRd
#' @keywords internal
.check_regime_init <- function(init_prob, k) {
  if (is.null(init_prob)) {
    return(rep(1 / k, k))
  }
  init_prob <- as.numeric(init_prob)
  if (length(init_prob) != k || any(init_prob < 0) || anyNA(init_prob)) {
    cli::cli_abort("`init_prob` must be {k} non-negative, non-missing values.")
  }
  if (!isTRUE(all.equal(sum(init_prob), 1))) {
    cli::cli_abort("`init_prob` must sum to one.")
  }
  init_prob
}
