# ssm-class.R -- the linear-Gaussian state-space model specification.
#
# A single S7 class, `ssm`, carries the matrices of a discrete-time
# linear-Gaussian state-space model and is the common substrate the Kalman
# filter, the Rauch-Tung-Striebel smoother and the regime-switching mixture
# filter all run on. System matrices may be supplied static (a single matrix,
# reused at every step) or time-varying (a list with one matrix per step). The
# constructor normalises scalars and vectors to matrices, wraps a static matrix
# in a length-one list, and dimension-checks every form so the recursions can
# index a per-step matrix without re-checking shapes.
#
# Notation follows Durbin and Koopman (2012): a latent state evolves as
# x_t = A_t x_{t-1} + w_t with w_t ~ N(0, Q_t), and is observed through
# y_t = B_t x_t + v_t with v_t ~ N(0, R_t).

# The ssm class ---------------------------------------------------------------

#' A linear-Gaussian state-space model
#'
#' An S7 object specifying a discrete-time linear-Gaussian state-space model,
#' built by the constructor of the same name. The latent state of dimension
#' `m` evolves through a transition matrix, an optional intercept (control
#' input) and Gaussian process noise, and is observed through an observation
#' matrix, an optional observation intercept and Gaussian observation noise:
#' \deqn{x_t = A_t\, x_{t-1} + c_t + w_t, \qquad w_t \sim N(0, Q_t),}
#' \deqn{y_t = B_t\, x_t + d_t + v_t, \qquad v_t \sim N(0, R_t).}
#' The intercepts default to zero, recovering the intercept-free form; a
#' non-zero \eqn{c_t} is the standard way to give a stationary state a
#' non-zero long-run mean (for an AR(1) state with coefficient \eqn{\phi} and
#' mean \eqn{\theta}, set \eqn{c = \theta(1 - \phi)}).
#'
#' Each of `transition`, `observation`, `state_cov` and `obs_cov` is accepted
#' as a single matrix (static, reused at every step) or as a list of matrices,
#' one per time step (time-varying). Scalars and length-`m` vectors are
#' promoted to matrices, so a univariate local-level model can be written with
#' bare numbers. Internally every system matrix is stored as a list of one or
#' more matrices; the filtering and smoothing recursions index it per step.
#'
#' The dimensions must be mutually consistent: with state dimension `m` and
#' observation dimension `d`, the transition matrices are `m` by `m`, the
#' observation matrices `d` by `m`, the process covariances `m` by `m`, and the
#' observation covariances `d` by `d`. The initial covariance is `m` by `m`.
#'
#' @param transition The state transition matrix \eqn{A} (`m` by `m`), or a
#'   list of such matrices for a time-varying transition.
#' @param observation The observation matrix \eqn{B} (`d` by `m`, where `d` is
#'   the observation dimension), or a list of such matrices.
#' @param state_cov The process-noise covariance \eqn{Q} (`m` by `m`, symmetric
#'   positive semi-definite), or a list of such matrices.
#' @param obs_cov The observation-noise covariance \eqn{R} (`d` by `d`,
#'   symmetric positive definite), or a list of such matrices.
#' @param state_intercept The state intercept \eqn{c} (a scalar or length-`m`
#'   vector, or a list of such for a time-varying control input). Defaults to
#'   `0`.
#' @param obs_intercept The observation intercept \eqn{d} (a scalar or
#'   length-`d` vector, or a list of such). Defaults to `0`.
#' @param init_state The prior mean of the initial state \eqn{x_0}, a numeric
#'   vector of length `m`.
#' @param init_cov The prior covariance of the initial state, an `m` by `m`
#'   symmetric positive semi-definite matrix.
#' @param obs_family The observation-noise family, `"gaussian"` (the default) or
#'   `"student_t"`. A Student-t family gives a heavy-tailed observation model
#'   whose filter is robust to outliers; `obs_cov` is then the t *scale* matrix
#'   \eqn{R} rather than a covariance. See [kalman_filter()] for the recursion.
#' @param obs_df The Student-t degrees of freedom \eqn{\nu} (a finite scalar
#'   greater than two) when `obs_family = "student_t"`; ignored, and held at
#'   `Inf`, for the Gaussian family. [estimate_obs_df()] selects it by profile
#'   likelihood. As \eqn{\nu \to \infty} the Student-t model collapses to the
#'   Gaussian one.
#'
#' @returns An S7 object of class `ssm`.
#' @family state-space
#' @seealso [kalman_filter()], [rts_smoother()], [mixture_filter()].
#' @references
#' Durbin, J. and Koopman, S. J. (2012). *Time Series Analysis by State Space
#' Methods*, 2nd ed. Oxford University Press.
#' @export
#' @examples
#' ## A local-level (random-walk-plus-noise) model written with bare numbers.
#' model <- ssm(
#'   transition = 1,
#'   observation = 1,
#'   state_cov = 0.1,
#'   obs_cov = 1,
#'   init_state = 0,
#'   init_cov = 10
#' )
#' model
#'
#' ## A constant-velocity tracking model with a vector state.
#' velocity <- ssm(
#'   transition = matrix(c(1, 0, 1, 1), nrow = 2),
#'   observation = matrix(c(1, 0), nrow = 1),
#'   state_cov = diag(c(0.01, 0.01)),
#'   obs_cov = 1,
#'   init_state = c(0, 0),
#'   init_cov = diag(c(10, 10))
#' )
#' velocity
ssm <- S7::new_class(
  name = "ssm",
  package = "kalmix",
  properties = list(
    transition = S7::class_list,
    observation = S7::class_list,
    state_cov = S7::class_list,
    obs_cov = S7::class_list,
    state_intercept = S7::class_list,
    obs_intercept = S7::class_list,
    init_state = S7::class_double,
    init_cov = S7::class_double,
    obs_family = S7::new_property(S7::class_character, default = "gaussian"),
    obs_df = S7::new_property(S7::class_double, default = Inf),
    state_dim = S7::new_property(
      class = S7::class_integer,
      getter = function(self) length(self@init_state)
    ),
    obs_dim = S7::new_property(
      class = S7::class_integer,
      getter = function(self) nrow(self@observation[[1L]])
    )
  ),
  constructor = function(transition,
                         observation,
                         state_cov,
                         obs_cov,
                         init_state,
                         init_cov,
                         state_intercept = 0,
                         obs_intercept = 0,
                         obs_family = c("gaussian", "student_t"),
                         obs_df = Inf) {
    obs_family <- match.arg(obs_family)
    init_state <- as.numeric(init_state)
    m <- length(init_state)
    if (m < 1L) {
      cli::cli_abort("`init_state` must have at least one element.")
    }
    init_cov <- .as_square_matrix(init_cov, m, "init_cov")
    obs_df <- if (identical(obs_family, "gaussian")) Inf else as.numeric(obs_df)
    observation <- .as_matrix_list(observation, "observation")
    d <- nrow(observation[[1L]])

    S7::new_object(
      S7::S7_object(),
      transition = .as_matrix_list(transition, "transition"),
      observation = observation,
      state_cov = .as_matrix_list(state_cov, "state_cov"),
      obs_cov = .as_matrix_list(obs_cov, "obs_cov"),
      state_intercept = .as_vector_list(state_intercept, m, "state_intercept"),
      obs_intercept = .as_vector_list(obs_intercept, d, "obs_intercept"),
      init_state = init_state,
      init_cov = init_cov,
      obs_family = obs_family,
      obs_df = obs_df
    )
  },
  validator = function(self) {
    .validate_ssm_properties(self)
  }
)

#' @export
S7::method(print, ssm) <- function(x, ...) {
  shape <- function(lst) if (length(lst) == 1L) "static" else "time-varying"
  kind <- if (identical(x@obs_family, "student_t")) {
    sprintf("Student-t obs (df %.3g)", x@obs_df)
  } else {
    "linear-Gaussian"
  }
  cat(sprintf(
    "<ssm>: %s state-space model (state %d, obs %d)\n",
    kind, x@state_dim, x@obs_dim
  ))
  cat(sprintf("  transition  : %s\n", shape(x@transition)))
  cat(sprintf("  observation : %s\n", shape(x@observation)))
  cat(sprintf("  state_cov   : %s\n", shape(x@state_cov)))
  cat(sprintf("  obs_cov     : %s\n", shape(x@obs_cov)))
  invisible(x)
}
