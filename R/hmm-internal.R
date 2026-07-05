# hmm-internal.R -- emission, log-sum-exp and validation for the HMM verbs.
#
# These helpers supply the per-state log-emission matrix the forward-backward
# and Viterbi recursions share, a numerically-stable log-sum-exp, and the
# validators that keep hmm(), hmm_filter() and hmm_viterbi() reading as their
# contracts.

# Numerical helpers -----------------------------------------------------------

#' Log-sum-exp of a numeric vector
#'
#' Computes \eqn{\log \sum_i \exp(x_i)} with the standard max-shift so the sum
#' does not overflow or underflow. The recursions accumulate probabilities in
#' the log domain and combine them with this.
#'
#' @param x A numeric vector of log-values.
#'
#' @returns Numeric scalar.
#' @noRd
#' @keywords internal
.logsumexp <- function(x) {
  shift <- max(x)
  if (!is.finite(shift)) {
    return(shift)
  }
  shift + log(sum(exp(x - shift)))
}

#' Per-state Gaussian log-emission matrix
#'
#' Builds the `n` by `k` matrix whose `(t, j)` entry is the log-density of
#' observation `t` under the Gaussian emission of state `j`. Shared by
#' [hmm_filter()] and [hmm_viterbi()].
#'
#' @param model An [hmm].
#' @param y The numeric observation vector.
#'
#' @returns An `n` by `k` numeric matrix of log-densities.
#' @noRd
#' @keywords internal
.hmm_log_emission <- function(model, y) {
  k <- model@n_states
  vapply(seq_len(k), function(j) {
    stats::dnorm(
      y,
      mean = model@emission_mean[j],
      sd = model@emission_sd[j],
      log = TRUE
    )
  }, numeric(length(y)))
}

# Validation helpers ----------------------------------------------------------

#' Validate the assembled properties of an hmm object
#'
#' Cross-checks that the initial distribution, transition matrix and emission
#' parameters all share one state count, that the probabilities are well-formed,
#' and that the emission standard deviations are strictly positive. Returns an
#' error string (S7 validator contract) or `NULL`.
#'
#' @param self An `hmm` object under construction.
#'
#' @returns `NULL` if valid, otherwise a character scalar describing the fault.
#' @noRd
#' @keywords internal
.validate_hmm_properties <- function(self) {
  k <- length(self@init_prob)
  if (k < 1L) {
    return("`init_prob` must have at least one state")
  }
  if (any(self@init_prob < 0) || anyNA(self@init_prob) ||
      !isTRUE(all.equal(sum(self@init_prob), 1))) {
    return("`init_prob` must be non-negative and sum to one")
  }
  if (nrow(self@transition) != k || ncol(self@transition) != k) {
    return(sprintf("`transition` must be a %d by %d matrix", k, k))
  }
  if (any(self@transition < 0) || anyNA(self@transition) ||
      !isTRUE(all.equal(unname(rowSums(self@transition)), rep(1, k)))) {
    return("each row of `transition` must be non-negative and sum to one")
  }
  if (length(self@emission_mean) != k) {
    return(sprintf("`emission_mean` must have %d entries", k))
  }
  if (length(self@emission_sd) != k || any(self@emission_sd <= 0) ||
      anyNA(self@emission_sd)) {
    return(sprintf("`emission_sd` must be %d positive values", k))
  }
  NULL
}

#' Validate an observation vector for the HMM verbs
#'
#' Checks that the observations form a finite numeric vector of at least two
#' points. The HMM verbs are univariate in this release.
#'
#' @param y The `y` argument to [hmm_filter()] or [hmm_viterbi()].
#'
#' @returns The observations as a plain numeric vector.
#' @noRd
#' @keywords internal
.check_hmm_observations <- function(y) {
  if (!is.numeric(y) || !is.null(dim(y))) {
    cli::cli_abort("`y` must be a numeric vector.")
  }
  if (length(y) < 2L) {
    cli::cli_abort("`y` must have at least two observations.")
  }
  if (anyNA(y)) {
    cli::cli_abort("`y` must not contain missing values.")
  }
  as.numeric(y)
}
