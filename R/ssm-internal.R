# ssm-internal.R -- matrix normalisation and validation for the ssm class.
#
# These helpers turn the user-friendly inputs accepted by ssm() (scalars,
# vectors, single matrices, or lists of matrices) into the canonical internal
# form -- a list of one or more square-or-rectangular matrices -- and then
# dimension-check the assembled model. Keeping them here keeps the ssm
# constructor and validator reading as their contract.

#' Coerce a covariance-or-system argument to a square matrix
#'
#' Promotes a scalar or length-`m` vector to a diagonal `m` by `m` matrix and
#' otherwise checks that a supplied matrix is `m` by `m`. Used for the initial
#' covariance, where a single square matrix is always expected.
#'
#' @param x A scalar, length-`m` vector, or `m` by `m` matrix.
#' @param m Integer scalar — the required dimension.
#' @param arg Character scalar — the argument name, for error messages.
#'
#' @returns An `m` by `m` numeric matrix.
#' @noRd
#' @keywords internal
.as_square_matrix <- function(x, m, arg) {
  if (is.matrix(x)) {
    if (nrow(x) != m || ncol(x) != m) {
      cli::cli_abort("`{arg}` must be a {m} by {m} matrix.")
    }
    return(x)
  }
  if (length(x) == 1L) {
    return(diag(as.numeric(x), nrow = m))
  }
  if (length(x) == m) {
    return(diag(as.numeric(x), nrow = m))
  }
  cli::cli_abort(
    "`{arg}` must be a scalar, a length-{m} vector, or a {m} by {m} matrix."
  )
}

#' Normalise a system argument to a list of matrices
#'
#' Accepts a single matrix (static), a list of matrices (time-varying), or a
#' bare scalar or numeric vector. A scalar becomes a one-by-one matrix and a
#' length-`m` vector a diagonal `m` by `m` matrix, so a diagonal covariance can
#' be written as a plain vector. Returns a list of matrices without yet
#' cross-checking dimensions against the rest of the model — that is the
#' validator's job once every component is in canonical form.
#'
#' @param x A matrix, a list of matrices, a scalar, or a numeric vector.
#' @param arg Character scalar — the argument name, for error messages.
#'
#' @returns A list of numeric matrices.
#' @noRd
#' @keywords internal
.as_matrix_list <- function(x, arg) {
  promote <- function(el) {
    if (is.matrix(el)) {
      return(el)
    }
    if (is.numeric(el) && length(el) == 1L) {
      return(matrix(as.numeric(el), nrow = 1L, ncol = 1L))
    }
    if (is.numeric(el) && length(el) > 1L) {
      return(diag(as.numeric(el), nrow = length(el)))
    }
    cli::cli_abort(
      "`{arg}` entries must be matrices, a scalar, or a numeric vector."
    )
  }
  if (is.list(x)) {
    if (length(x) == 0L) {
      cli::cli_abort("`{arg}` must not be an empty list.")
    }
    return(lapply(x, promote))
  }
  list(promote(x))
}

#' Normalise an intercept argument to a list of numeric vectors
#'
#' Accepts a scalar (recycled to length `len`), a length-`len` vector, or a
#' list of such (time-varying). Mirrors [.as_matrix_list()] for the additive
#' intercept terms, which are vectors rather than matrices.
#'
#' @param x A scalar, a numeric vector, or a list of numeric vectors.
#' @param len Integer scalar — the required vector length.
#' @param arg Character scalar — the argument name, for error messages.
#'
#' @returns A list of numeric vectors of length `len`.
#' @noRd
#' @keywords internal
.as_vector_list <- function(x, len, arg) {
  promote <- function(el) {
    if (!is.numeric(el)) {
      cli::cli_abort("`{arg}` entries must be numeric.")
    }
    el <- as.numeric(el)
    if (length(el) == 1L) {
      el <- rep(el, len)
    }
    if (length(el) != len) {
      cli::cli_abort("`{arg}` entries must be scalars or length-{len} vectors.")
    }
    el
  }
  if (is.list(x)) {
    if (length(x) == 0L) {
      cli::cli_abort("`{arg}` must not be an empty list.")
    }
    return(lapply(x, promote))
  }
  list(promote(x))
}

#' Validate the assembled properties of an ssm object
#'
#' Cross-checks that the state and observation dimensions implied by the system
#' matrices agree with each other and with the initial state, that every matrix
#' in a time-varying list shares the right shape, and that the covariance
#' matrices are square. Returns an error string (S7 validator contract) or
#' `NULL` when the model is well-formed.
#'
#' @param self An `ssm` object under construction.
#'
#' @returns `NULL` if valid, otherwise a character scalar describing the fault.
#' @noRd
#' @keywords internal
.validate_ssm_properties <- function(self) {
  m <- length(self@init_state)
  if (m < 1L) {
    return("`init_state` must have at least one element")
  }
  if (!all(is.finite(self@init_state))) {
    return("`init_state` must be finite")
  }
  if (nrow(self@init_cov) != m || ncol(self@init_cov) != m) {
    return("`init_cov` must be a square matrix matching `init_state`")
  }

  # Each transition matrix is m by m -----------------------------------------

  for (a in self@transition) {
    if (nrow(a) != m || ncol(a) != m) {
      return(sprintf("each `transition` matrix must be %d by %d", m, m))
    }
  }
  for (q in self@state_cov) {
    if (nrow(q) != m || ncol(q) != m) {
      return(sprintf("each `state_cov` matrix must be %d by %d", m, m))
    }
  }

  # Observation matrices are d by m and observation covariances d by d -------

  d <- nrow(self@observation[[1L]])
  for (b in self@observation) {
    if (nrow(b) != d || ncol(b) != m) {
      return(sprintf("each `observation` matrix must be %d by %d", d, m))
    }
  }
  for (r in self@obs_cov) {
    if (nrow(r) != d || ncol(r) != d) {
      return(sprintf("each `obs_cov` matrix must be %d by %d", d, d))
    }
  }

  # Intercepts are length-m (state) and length-d (observation) vectors --------

  for (ci in self@state_intercept) {
    if (length(ci) != m || !all(is.finite(ci))) {
      return(sprintf("each `state_intercept` must be a finite length-%d vector", m))
    }
  }
  for (di in self@obs_intercept) {
    if (length(di) != d || !all(is.finite(di))) {
      return(sprintf("each `obs_intercept` must be a finite length-%d vector", d))
    }
  }

  # Observation family: gaussian (df = Inf) or student_t (finite df > 2, so the
  # observation variance -- and the F-reference of the scale diagnostic -- exist.
  if (length(self@obs_family) != 1L ||
      !self@obs_family %in% c("gaussian", "student_t")) {
    return("`obs_family` must be \"gaussian\" or \"student_t\"")
  }
  if (identical(self@obs_family, "gaussian")) {
    if (!(length(self@obs_df) == 1L && is.infinite(self@obs_df))) {
      return("a gaussian `obs_family` requires `obs_df` = Inf")
    }
  } else {
    if (length(self@obs_df) != 1L || !is.finite(self@obs_df) || self@obs_df <= 2) {
      return("a student_t `obs_family` requires a finite `obs_df` greater than 2")
    }
  }
  NULL
}

#' Index a system-matrix list at a given time step
#'
#' Returns the single matrix for a static component (length-one list) or the
#' step-`t` matrix for a time-varying component. The filtering and smoothing
#' recursions call this so they need not branch on static vs time-varying.
#'
#' @param lst A list of one or more matrices, as stored on an `ssm`.
#' @param t Integer scalar — the one-based time step.
#'
#' @returns A numeric matrix.
#' @noRd
#' @keywords internal
.at_step <- function(lst, t) {
  if (length(lst) == 1L) lst[[1L]] else lst[[t]]
}

#' Validate an observation matrix passed to a filter against a model
#'
#' Checks that the observations supplied to [kalman_filter()] (and the other
#' recursions) form a numeric matrix with the model's observation dimension in
#' its columns, and that a time-varying model has exactly as many system
#' matrices as there are observation rows. Missing observations are admitted
#' row-wise: a fully `NA` row is a prediction-only step for the filters, while
#' a partially observed row is refused (a sub-vector update is not supported),
#' so a half-missing row cannot silently misbehave.
#'
#' @param y A numeric vector or matrix of observations.
#' @param model An [ssm].
#'
#' @returns The observations as an `n` by `d` numeric matrix.
#' @noRd
#' @keywords internal
.check_observations <- function(y, model) {
  if (is.numeric(y) && is.null(dim(y))) {
    y <- matrix(y, ncol = 1L)
  }
  if (!is.numeric(y) || !is.matrix(y)) {
    cli::cli_abort("`y` must be a numeric vector or matrix of observations.")
  }
  d <- model@obs_dim
  if (ncol(y) != d) {
    cli::cli_abort(
      "`y` must have {d} column{?s} to match the model's observation dimension."
    )
  }
  if (anyNA(y) && d > 1L) {
    n_missing <- rowSums(is.na(y))
    partial <- n_missing > 0L & n_missing < d
    if (any(partial)) {
      cli::cli_abort(c(
        "`y` rows must be fully observed or fully missing.",
        "i" = "Partially missing row{?s}: {which(partial)}.",
        "i" = "A fully `NA` row is skipped (prediction only); sub-vector
               updates for partially observed rows are not supported."
      ))
    }
  }
  n <- nrow(y)
  for (nm in c("transition", "observation", "state_cov", "obs_cov",
               "state_intercept", "obs_intercept")) {
    lst <- S7::prop(model, nm)
    if (length(lst) != 1L && length(lst) != n) {
      cli::cli_abort(c(
        "Time-varying `{nm}` has {length(lst)} entries but there are {n} observations.",
        "i" = "A time-varying component needs one entry per observation."
      ))
    }
  }
  y
}
