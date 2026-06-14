# student-t.R -- selecting the Student-t degrees of freedom by profile likelihood.
#
# The robust observation model has one free shape parameter, the degrees of
# freedom nu. estimate_obs_df() profiles the model log-likelihood over a grid of
# nu (Inf included, so a series with no heavy tails selects the Gaussian model)
# and returns the maximiser. The caller builds the chosen ssm with
# obs_family = "student_t" and the returned df, so the adequacy guard then reads
# the verdict off the correctly specified family.

#' Choose the Student-t degrees of freedom by profile likelihood
#'
#' Profiles the one-step-ahead predictive log-likelihood of a state-space model
#' over a grid of Student-t degrees of freedom and returns the maximiser. The
#' grid includes `Inf` (the Gaussian model), so a series whose innovations are
#' not heavy-tailed selects `Inf` and the model collapses to the linear-Gaussian
#' filter -- the estimator never imposes heavy tails that the data do not show.
#'
#' The supplied `model` provides the dynamics (transition, process noise, prior)
#' and the observation *scale* `obs_cov`; only the observation family and its
#' degrees of freedom are varied. The log-likelihood is the Student-t predictive
#' density accumulated by [kalman_filter()] (an approximation that treats the
#' one-step-ahead distribution as Student-t; see [kalman_filter()]).
#'
#' @param model An [ssm] giving the dynamics and the observation scale. Its own
#'   `obs_family` is ignored; each grid value is evaluated on a clone.
#' @param y A numeric vector or `n` by `d` matrix of observations.
#' @param df_grid A numeric vector of candidate degrees of freedom to profile
#'   over. Must be positive; may include `Inf`. Defaults to
#'   `c(3, 4, 5, 7, 10, 15, 30, Inf)`.
#' @param lrt_margin The parsimony margin in log-likelihood units. The Gaussian
#'   model (`df = Inf`) is preferred unless a finite `df` improves the
#'   log-likelihood by more than this much, so a finite Gaussian sample -- whose
#'   profile peaks at a large but finite `df` by chance -- still selects the
#'   Gaussian. Defaults to half the 95% one-degree-of-freedom chi-squared
#'   quantile (a likelihood-ratio test at the 5% level), about `1.92`.
#'
#' @returns A list with `df` (the selected degrees of freedom), `family`
#'   (`"gaussian"` if `df` is `Inf`, else `"student_t"`), `log_lik` (the
#'   maximised log-likelihood), `grid` (the grid evaluated) and `profile` (the
#'   log-likelihood at each grid value).
#' @family causal
#' @seealso [ssm()] for building the chosen model; [aci()] for the causal
#'   read-out the chosen family grounds.
#' @references
#' Roth, M., Ozkan, E. and Gustafsson, F. (2013). A Student's t filter for heavy
#' tailed process and measurement noise. *ICASSP*, 5770--5774.
#' @export
#' @examples
#' ## A heavy-tailed series: the profile prefers a finite df over the Gaussian.
#' model <- ssm(
#'   transition = 1, observation = 1,
#'   state_cov = 0.04, obs_cov = 1,
#'   init_state = 0, init_cov = 10
#' )
#' set.seed(1)
#' level <- cumsum(rnorm(300, sd = 0.2))
#' y <- level + rt(300, df = 3)
#' est <- estimate_obs_df(model, y)
#' est$df
estimate_obs_df <- function(model, y,
                            df_grid = c(3, 4, 5, 7, 10, 15, 30, Inf),
                            lrt_margin = stats::qchisq(0.95, 1) / 2) {
  if (!S7::S7_inherits(model, ssm)) {
    cli::cli_abort("`model` must be an {.cls ssm}.")
  }
  df_grid <- as.numeric(df_grid)
  if (length(df_grid) < 1L || any(is.na(df_grid)) || any(df_grid <= 2)) {
    cli::cli_abort("`df_grid` must be positive degrees of freedom greater than 2.")
  }

  profile <- vapply(df_grid, function(nu) {
    family <- if (is.infinite(nu)) "gaussian" else "student_t"
    candidate <- .with_family(model, family, nu)
    tryCatch(kalman_filter(candidate, y)@log_lik, error = function(e) NA_real_)
  }, numeric(1L))

  if (all(!is.finite(profile))) {
    cli::cli_abort("No candidate degrees of freedom yielded a finite likelihood.")
  }
  best <- which.max(profile)
  # Parsimony: prefer the Gaussian model unless a finite df beats it by more than
  # the likelihood-ratio margin, so light-tailed data are not forced heavy-tailed.
  inf_idx <- which(is.infinite(df_grid))
  if (length(inf_idx) == 1L && is.finite(profile[inf_idx]) &&
      profile[best] - profile[inf_idx] < lrt_margin) {
    best <- inf_idx
  }
  df <- df_grid[best]

  list(
    df = df,
    family = if (is.infinite(df)) "gaussian" else "student_t",
    log_lik = profile[best],
    grid = df_grid,
    profile = profile
  )
}
