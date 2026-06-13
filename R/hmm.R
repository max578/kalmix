# hmm.R -- discrete-state hidden Markov model filtering and Viterbi decoding.
#
# A hidden Markov model is the discrete-state cousin of the linear-Gaussian
# state-space model: a finite Markov chain over hidden states emits an
# observation at each step. hmm() builds the specification (initial
# distribution, transition matrix, per-state Gaussian emission), hmm_filter()
# runs the forward-backward recursion for the filtered and smoothed state
# posteriors and the log-likelihood, and hmm_viterbi() returns the single most
# likely state path. All three work in the log domain for numerical stability.

#' A hidden Markov model specification
#'
#' An S7 object specifying a discrete-state hidden Markov model with univariate
#' Gaussian emissions, as built by the constructor of the same name. The hidden
#' state follows a Markov chain with the given initial distribution and
#' transition matrix; in state `j` the observation is Gaussian with mean
#' `emission_mean[j]` and standard deviation `emission_sd[j]`.
#'
#' @param init_prob A length-`k` vector of initial state probabilities, summing
#'   to one.
#' @param transition A `k` by `k` row-stochastic transition matrix.
#' @param emission_mean A length-`k` vector of per-state emission means.
#' @param emission_sd A length-`k` vector of per-state emission standard
#'   deviations, all strictly positive.
#'
#' @returns An S7 object of class `hmm`.
#' @family hidden-markov
#' @seealso [hmm_filter()], [hmm_viterbi()].
#' @export
#' @examples
#' ## A two-state mean-shift model: a low state and a high state.
#' model <- hmm(
#'   init_prob = c(0.5, 0.5),
#'   transition = matrix(c(0.95, 0.05, 0.05, 0.95), nrow = 2, byrow = TRUE),
#'   emission_mean = c(0, 3),
#'   emission_sd = c(1, 1)
#' )
#' model
hmm <- S7::new_class(
  name = "hmm",
  package = "kalmix",
  properties = list(
    init_prob = S7::class_double,
    transition = S7::class_double,
    emission_mean = S7::class_double,
    emission_sd = S7::class_double,
    n_states = S7::new_property(
      class = S7::class_integer,
      getter = function(self) length(self@init_prob)
    )
  ),
  validator = function(self) {
    .validate_hmm_properties(self)
  }
)

#' Hidden Markov model filtering and smoothing
#'
#' Runs the forward-backward recursion over an [hmm]. The forward pass gives
#' the filtered state posteriors \eqn{P(s_t = j \mid y_{1:t})} and the model
#' log-likelihood; the backward pass combines with it to give the smoothed
#' posteriors \eqn{P(s_t = j \mid y_{1:n})}. The recursion is carried in the
#' log domain and renormalised each step, so it is stable over long series.
#'
#' @param model An [hmm].
#' @param y A numeric vector of observations.
#'
#' @returns An [hmm_fit].
#' @family hidden-markov
#' @seealso [hmm_viterbi()] for the most likely path.
#' @references
#' Rabiner, L. R. (1989). A tutorial on hidden Markov models and selected
#' applications in speech recognition. *Proceedings of the IEEE*, 77(2),
#' 257--286.
#' @export
#' @examples
#' model <- hmm(
#'   init_prob = c(0.5, 0.5),
#'   transition = matrix(c(0.95, 0.05, 0.05, 0.95), nrow = 2, byrow = TRUE),
#'   emission_mean = c(0, 3), emission_sd = c(1, 1)
#' )
#' set.seed(1)
#' y <- c(rnorm(40, 0), rnorm(40, 3))
#' fit <- hmm_filter(model, y)
#' fit@log_lik
hmm_filter <- function(model, y) {
  if (!S7::S7_inherits(model, hmm)) {
    cli::cli_abort("`model` must be an {.cls hmm}.")
  }
  y <- .check_hmm_observations(y)
  n <- length(y)
  k <- model@n_states

  log_emit <- .hmm_log_emission(model, y)
  log_trans <- log(model@transition + .Machine$double.xmin)
  log_init <- log(model@init_prob + .Machine$double.xmin)

  # Forward pass (filtered posteriors and log-likelihood) --------------------

  log_alpha <- matrix(0, nrow = n, ncol = k)
  filtered <- matrix(0, nrow = n, ncol = k)
  log_lik <- 0

  log_alpha[1L, ] <- log_init + log_emit[1L, ]
  norm <- .softmax_with_norm(log_alpha[1L, ])
  filtered[1L, ] <- norm$prob
  log_lik <- log_lik + norm$log_norm
  log_alpha[1L, ] <- log(norm$prob + .Machine$double.xmin)

  for (t in 2:n) {
    for (j in seq_len(k)) {
      log_alpha[t, j] <- .logsumexp(log_alpha[t - 1L, ] + log_trans[, j]) +
        log_emit[t, j]
    }
    norm <- .softmax_with_norm(log_alpha[t, ])
    filtered[t, ] <- norm$prob
    log_lik <- log_lik + norm$log_norm
    log_alpha[t, ] <- log(norm$prob + .Machine$double.xmin)
  }

  # Backward pass and smoothing ----------------------------------------------

  log_beta <- matrix(0, nrow = n, ncol = k)
  for (t in rev(seq_len(n - 1L))) {
    for (i in seq_len(k)) {
      log_beta[t, i] <- .logsumexp(
        log_trans[i, ] + log_emit[t + 1L, ] + log_beta[t + 1L, ]
      )
    }
    log_beta[t, ] <- log_beta[t, ] - max(log_beta[t, ])
  }

  smoothed <- matrix(0, nrow = n, ncol = k)
  for (t in seq_len(n)) {
    smoothed[t, ] <- .softmax_with_norm(log_alpha[t, ] + log_beta[t, ])$prob
  }

  hmm_fit(
    filtered = filtered,
    smoothed = smoothed,
    state_path = max.col(smoothed, ties.method = "first"),
    log_lik = log_lik,
    model = model
  )
}

#' Most likely hidden-state path (Viterbi)
#'
#' Returns the single most likely sequence of hidden states given the whole
#' observed series, by the Viterbi (1967) dynamic-programming recursion in the
#' log domain. Unlike the per-step smoothed mode from [hmm_filter()], the
#' Viterbi path is jointly optimal: it is the state sequence that maximises the
#' joint posterior, and so respects the transition structure across steps.
#'
#' @param model An [hmm].
#' @param y A numeric vector of observations.
#'
#' @returns An integer vector of length `length(y)`, the most likely state
#'   index at each step.
#' @family hidden-markov
#' @seealso [hmm_filter()].
#' @references
#' Viterbi, A. (1967). Error bounds for convolutional codes and an
#' asymptotically optimum decoding algorithm. *IEEE Transactions on Information
#' Theory*, 13(2), 260--269.
#' @export
#' @examples
#' model <- hmm(
#'   init_prob = c(0.5, 0.5),
#'   transition = matrix(c(0.95, 0.05, 0.05, 0.95), nrow = 2, byrow = TRUE),
#'   emission_mean = c(0, 3), emission_sd = c(1, 1)
#' )
#' set.seed(1)
#' y <- c(rnorm(40, 0), rnorm(40, 3))
#' table(hmm_viterbi(model, y))
hmm_viterbi <- function(model, y) {
  if (!S7::S7_inherits(model, hmm)) {
    cli::cli_abort("`model` must be an {.cls hmm}.")
  }
  y <- .check_hmm_observations(y)
  n <- length(y)
  k <- model@n_states

  log_emit <- .hmm_log_emission(model, y)
  log_trans <- log(model@transition + .Machine$double.xmin)
  log_init <- log(model@init_prob + .Machine$double.xmin)

  delta <- matrix(-Inf, nrow = n, ncol = k)
  backptr <- matrix(0L, nrow = n, ncol = k)
  delta[1L, ] <- log_init + log_emit[1L, ]

  for (t in 2:n) {
    for (j in seq_len(k)) {
      scores <- delta[t - 1L, ] + log_trans[, j]
      backptr[t, j] <- which.max(scores)
      delta[t, j] <- max(scores) + log_emit[t, j]
    }
  }

  # Backtrack the optimal path -----------------------------------------------

  path <- integer(n)
  path[n] <- which.max(delta[n, ])
  for (t in rev(seq_len(n - 1L))) {
    path[t] <- backptr[t + 1L, path[t + 1L]]
  }
  path
}

#' A hidden Markov model fit
#'
#' An S7 object holding the output of [hmm_filter()]: the filtered and smoothed
#' state posteriors, the smoothed most-likely state at each step, and the model
#' log-likelihood.
#'
#' @param filtered An `n` by `k` matrix of filtered state posteriors.
#' @param smoothed An `n` by `k` matrix of smoothed state posteriors.
#' @param state_path An integer vector of the per-step smoothed modal state.
#' @param log_lik Numeric scalar — the model log-likelihood.
#' @param model The [hmm] that was filtered.
#'
#' @returns An S7 object of class `hmm_fit`.
#' @family hidden-markov
#' @export
hmm_fit <- S7::new_class(
  name = "hmm_fit",
  package = "kalmix",
  properties = list(
    filtered = S7::class_double,
    smoothed = S7::class_double,
    state_path = S7::class_integer,
    log_lik = S7::class_double,
    model = hmm
  )
)

#' @export
S7::method(print, hmm) <- function(x, ...) {
  cat(sprintf("<hmm>: %d-state Gaussian hidden Markov model\n", x@n_states))
  cat("  emission means :", paste(round(x@emission_mean, 3), collapse = ", "),
      "\n")
  cat("  emission sds   :", paste(round(x@emission_sd, 3), collapse = ", "),
      "\n")
  invisible(x)
}

#' @export
S7::method(print, hmm_fit) <- function(x, ...) {
  n <- nrow(x@filtered)
  cat(sprintf(
    "<hmm_fit>: %d step%s, %d states, log-likelihood %.3f\n",
    n, if (n == 1L) "" else "s", ncol(x@filtered), x@log_lik
  ))
  invisible(x)
}
