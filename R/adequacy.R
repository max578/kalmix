# adequacy.R -- model-adequacy diagnostics on the Kalman innovations, and the
# grounding read-out that lets Assimilative Causal Inference abstain.
#
# Assimilative Causal Inference reads a causal verdict off a state-space model
# under the maintained assumption that the model is the data-generating
# mechanism; it does not, on its own, test that assumption. A verdict computed
# from a model that cannot even describe the observed series is not an honest
# verdict, so kalmix runs the standard innovation diagnostics of state-space
# model checking (Harvey 1989; Durbin and Koopman 2012) over the one-step-ahead
# prediction errors and turns the result into an abstention: when the
# standardised innovations are not white, not correctly scaled, or not Gaussian,
# the model is inadequate and the causal read-out is labelled `[unverified]`.
#
# Model adequacy is a necessary, not a sufficient, condition for a grounded
# causal claim. Passing the diagnostics shows only that the model is a
# self-consistent description of the data; it does not verify the model against
# an external authority, which is what grounding means (provenance_vocabulary.md
# section 1: self-consistency does not ground a fact). The two axes are kept
# apart here: innovation_diagnostics() is the internal statistical check, and
# the grounding token is set to "grounded" only when an adequate model is also
# accompanied by declared, dated mechanism provenance. The default verdict is
# the honest "[unverified]".

# Canonical grounding tokens ---------------------------------------------------

# The two grounding states of provenance_vocabulary.md, hardcoded because kalmix
# is standalone (invariant 1) and cannot source the orchestra's
# integration/orchestra_manifest.R. The square brackets on the unverified token
# are load-bearing: they make an un-grounded label unmissable in a verdict line.

.grounding_grounded <- "grounded"
.grounding_unverified <- "[unverified]"

#' Worst-case combination of grounding tokens
#'
#' Combines grounding tokens so that a verdict built from several inputs is
#' grounded only if every input is grounded -- one un-grounded input taints the
#' whole verdict (provenance_vocabulary.md section 2b).
#'
#' @param ... Grounding token strings.
#'
#' @returns A single grounding token.
#' @noRd
#' @keywords internal
.combine_grounding <- function(...) {
  tokens <- c(...)
  if (length(tokens) == 0L || any(tokens != .grounding_grounded)) {
    .grounding_unverified
  } else {
    .grounding_grounded
  }
}

#' Is the supplied mechanism provenance grounded?
#'
#' A mechanism is grounded only when its provenance carries a non-`NA` date on
#' which the model was checked against the external authority that owns it
#' (provenance_vocabulary.md section 2a: a fact is grounded if and only if
#' `verified_on` is a non-`NA` `Date`).
#'
#' @param mechanism `NULL`, or a list with a `verified_on` `Date` entry.
#'
#' @returns A logical scalar.
#' @noRd
#' @keywords internal
.mechanism_grounded <- function(mechanism) {
  if (is.null(mechanism)) {
    return(FALSE)
  }
  vo <- mechanism[["verified_on"]]
  !is.null(vo) && length(vo) == 1L && inherits(vo, "Date") && !is.na(vo)
}

# Innovation diagnostics -------------------------------------------------------

#' Standardised one-step-ahead innovations
#'
#' Whitens the Kalman innovations by their own predicted covariance,
#' \eqn{w_t = S_t^{-1/2} e_t} through the Cholesky factor of \eqn{S_t}, so that
#' under a correctly specified linear-Gaussian model the rows of the result are
#' independent standard Gaussian vectors. This is the quantity every
#' state-space residual diagnostic is built on.
#'
#' @param filter A [kalman_fit].
#'
#' @returns An `n` by `d` matrix of standardised innovations.
#' @noRd
#' @keywords internal
.standardised_innovations <- function(filter) {
  e <- filter@innovation
  covs <- filter@innovation_cov
  n <- nrow(e)
  d <- ncol(e)
  w <- matrix(0, nrow = n, ncol = d)
  for (t in seq_len(n)) {
    s_chol <- .safe_chol(covs[[t]], "innovation covariance")
    w[t, ] <- backsolve(s_chol, e[t, ], transpose = TRUE)
  }
  w
}

#' Jarque-Bera normality p-value
#'
#' The Jarque and Bera (1980) test statistic
#' \eqn{n(\hat{S}^2/6 + (\hat{K}-3)^2/24)} is asymptotically \eqn{\chi^2_2}
#' under normality; its upper-tail probability is returned. Computed natively
#' from the sample skewness and kurtosis so the package takes no dependency
#' for it.
#'
#' @param v A numeric vector.
#'
#' @returns A numeric p-value, or `NA_real_` when there are too few finite
#'   values to form the statistic.
#' @noRd
#' @keywords internal
.jarque_bera_p <- function(v) {
  v <- v[is.finite(v)]
  n <- length(v)
  if (n < 8L) {
    return(NA_real_)
  }
  centred <- v - mean(v)
  s2 <- mean(centred^2)
  if (s2 <= 0) {
    return(NA_real_)
  }
  skew <- mean(centred^3) / s2^1.5
  kurt <- mean(centred^4) / s2^2
  jb <- n * (skew^2 / 6 + (kurt - 3)^2 / 24)
  stats::pchisq(jb, df = 2, lower.tail = FALSE)
}

#' Innovation diagnostics for a fitted state-space model
#'
#' An S7 object holding the model-adequacy diagnostics of
#' [innovation_diagnostics()]: the standardised one-step-ahead innovations and
#' the three residual tests read off them -- whiteness, scale calibration and
#' normality -- together with the combined adequacy verdict.
#'
#' @param standardised An `n` by `d` matrix of standardised innovations.
#' @param whiteness_p Numeric scalar -- the Ljung-Box whiteness p-value,
#'   Bonferroni-combined across the observation dimensions.
#' @param calibration_p Numeric scalar -- the two-sided p-value of the
#'   normalised-innovation-squared scale test.
#' @param normality_p Numeric scalar -- the Jarque-Bera normality p-value of the
#'   pooled standardised innovations.
#' @param p_value Numeric scalar -- the combined adequacy p-value, the
#'   Bonferroni minimum across the three test families.
#' @param adequate Logical scalar -- `TRUE` when no test rejects at `alpha`,
#'   `FALSE` when the model is rejected, `NA` when there is too little data to
#'   judge.
#' @param alpha Numeric scalar -- the per-battery significance level.
#' @param lags Integer scalar -- the Ljung-Box lag.
#' @param n_obs Integer scalar -- the series length.
#' @param obs_dim Integer scalar -- the observation dimension.
#' @param obs_family Character scalar -- the observation-noise family the
#'   diagnostics were taken under (`"gaussian"` or `"student_t"`), which selects
#'   the reference distribution of the scale and family-fit tests.
#' @param reason Character scalar -- a short label for the verdict
#'   (`"adequate"`, `"model_inadequate"` or `"insufficient_data"`).
#'
#' @returns An S7 object of class `innov_diag`.
#' @family causal
#' @export
innov_diag <- S7::new_class(
  name = "innov_diag",
  package = "kalmix",
  properties = list(
    standardised = S7::class_double,
    whiteness_p = S7::new_property(S7::class_double, default = NA_real_),
    calibration_p = S7::new_property(S7::class_double, default = NA_real_),
    normality_p = S7::new_property(S7::class_double, default = NA_real_),
    p_value = S7::new_property(S7::class_double, default = NA_real_),
    adequate = S7::new_property(S7::class_logical, default = NA),
    alpha = S7::new_property(S7::class_double, default = 0.05),
    lags = S7::new_property(S7::class_integer, default = NA_integer_),
    n_obs = S7::new_property(S7::class_integer, default = NA_integer_),
    obs_dim = S7::new_property(S7::class_integer, default = NA_integer_),
    obs_family = S7::new_property(S7::class_character, default = "gaussian"),
    reason = S7::new_property(S7::class_character, default = NA_character_)
  )
)

#' Model-adequacy diagnostics for a fitted state-space model
#'
#' Runs the standard innovation diagnostics of linear-Gaussian state-space model
#' checking over a completed Kalman pass and reduces them to a single adequacy
#' verdict. Under a correctly specified model the standardised one-step-ahead
#' innovations are independent standard Gaussian, so three departures from that
#' are tested: serial correlation (a Ljung-Box whiteness test on each
#' observation dimension), a mis-scaled noise level (a two-sided test of the
#' total normalised innovation squared against its \eqn{\chi^2} reference) and
#' non-Gaussianity (a Jarque-Bera test on the pooled innovations). The model is
#' declared inadequate when any family rejects at `alpha`, combined
#' conservatively across the families by the Bonferroni rule.
#'
#' The diagnostics are what [aci()] uses to decide whether to trust its causal
#' read-out: a model whose innovations reject it is not a credible
#' data-generating mechanism, so the read-out is labelled `[unverified]` and the
#' decision abstains. Passing the diagnostics establishes self-consistency only,
#' which is necessary but not sufficient for a grounded causal claim (see
#' [aci()] on declaring mechanism provenance).
#'
#' @param object A [kalman_fit].
#' @param lags Integer scalar -- the Ljung-Box lag. Defaults to `NULL`, which
#'   uses one fifth of the series length, capped between one and ten.
#' @param alpha Numeric scalar -- the per-battery significance level the
#'   adequacy verdict is taken at. Defaults to `0.05`.
#'
#' @returns An [innov_diag].
#' @family causal
#' @seealso [aci()], which calls this to ground its causal verdict.
#' @references
#' Harvey, A. C. (1989). *Forecasting, Structural Time Series Models and the
#' Kalman Filter*. Cambridge University Press.
#'
#' Durbin, J. and Koopman, S. J. (2012). *Time Series Analysis by State Space
#' Methods*. 2nd ed. Oxford University Press.
#'
#' Jarque, C. M. and Bera, A. K. (1980). Efficient tests for normality,
#' homoscedasticity and serial independence of regression residuals.
#' *Economics Letters*, 6(3), 255--259.
#' @export
#' @examples
#' ## A correctly specified local level: the innovations pass the diagnostics.
#' model <- ssm(
#'   transition = 1, observation = 1,
#'   state_cov = 0.04, obs_cov = 1,
#'   init_state = 0, init_cov = 10
#' )
#' set.seed(1)
#' level <- cumsum(rnorm(300, sd = 0.2))
#' y <- level + rnorm(300)
#' innovation_diagnostics(kalman_filter(model, y))
innovation_diagnostics <- function(object, lags = NULL, alpha = 0.05) {
  if (!S7::S7_inherits(object, kalman_fit)) {
    cli::cli_abort("`object` must be a {.cls kalman_fit}.")
  }
  if (length(alpha) != 1L || !is.finite(alpha) || alpha <= 0 || alpha >= 1) {
    cli::cli_abort("`alpha` must be a scalar in (0, 1).")
  }

  w <- .standardised_innovations(object)
  n <- nrow(w)
  d <- ncol(w)

  if (is.null(lags)) {
    lags <- max(1L, min(10L, as.integer(n / 5L)))
  }
  lags <- as.integer(lags)
  if (length(lags) != 1L || is.na(lags) || lags < 1L) {
    cli::cli_abort("`lags` must be a positive integer scalar.")
  }

  model <- object@model
  is_t <- identical(model@obs_family, "student_t") && is.finite(model@obs_df)
  family <- if (is_t) "student_t" else "gaussian"

  # Too short to judge: fail closed to an indeterminate verdict rather than
  # assert adequacy the data cannot support.
  if (n < max(20L, 3L * lags)) {
    return(innov_diag(
      standardised = w, lags = lags, alpha = alpha,
      n_obs = n, obs_dim = d, obs_family = family,
      adequate = NA, reason = "insufficient_data"
    ))
  }

  # Whiteness: a Ljung-Box test per observation dimension, Bonferroni-combined.
  # Serial correlation is family-agnostic, so this test is unchanged for both.
  white_each <- vapply(
    seq_len(d),
    function(j) stats::Box.test(w[, j], lag = lags, type = "Ljung-Box")$p.value,
    numeric(1L)
  )
  p_white <- min(min(white_each) * d, 1)

  if (is_t) {
    # Student-t family: the standardised innovations are t-, not Gaussian-,
    # distributed, so the scale and family-fit tests read off the t-PIT.
    tt <- .t_scale_family_tests(w, model@obs_df)
    p_cal <- tt$scale_p
    p_norm <- tt$family_p
  } else {
    # Scale: the total normalised innovation squared is chi-squared with n*d
    # degrees of freedom under correct specification; a two-sided tail catches an
    # over- or under-stated noise level.
    nis <- sum(w^2)
    df_nis <- n * d
    p_cal <- 2 * min(
      stats::pchisq(nis, df_nis),
      stats::pchisq(nis, df_nis, lower.tail = FALSE)
    )

    # Normality of the pooled standardised innovations.
    p_norm <- .jarque_bera_p(as.numeric(w))
  }

  finite_p <- c(p_white, p_cal, p_norm)
  finite_p <- finite_p[is.finite(finite_p)]
  p_value <- min(min(finite_p) * length(finite_p), 1)
  adequate <- p_value >= alpha

  innov_diag(
    standardised = w,
    whiteness_p = p_white,
    calibration_p = p_cal,
    normality_p = p_norm,
    p_value = p_value,
    adequate = adequate,
    alpha = alpha,
    lags = lags,
    n_obs = n,
    obs_dim = d,
    obs_family = family,
    reason = if (adequate) "adequate" else "model_inadequate"
  )
}

#' The grounding verdict of an Assimilative Causal Inference read-out
#'
#' Combines the model-adequacy diagnostics with any declared mechanism
#' provenance into a grounding token and a short reason. An inadequate model
#' abstains outright; an adequate model is grounded only when its provenance has
#' been verified against an external authority, and is otherwise the honest
#' `[unverified]`.
#'
#' @param adequacy An [innov_diag].
#' @param mechanism `NULL`, or a list with a `verified_on` `Date`.
#'
#' @returns A list with `grounding` and `reason` character scalars.
#' @noRd
#' @keywords internal
.aci_grounding <- function(adequacy, mechanism) {
  if (isTRUE(adequacy@adequate)) {
    if (.mechanism_grounded(mechanism)) {
      list(grounding = .grounding_grounded, reason = "grounded")
    } else {
      list(grounding = .grounding_unverified, reason = "mechanism_unverified")
    }
  } else if (is.na(adequacy@adequate)) {
    list(grounding = .grounding_unverified, reason = adequacy@reason)
  } else {
    list(grounding = .grounding_unverified, reason = "model_inadequate")
  }
}

#' @export
S7::method(print, innov_diag) <- function(x, ...) {
  verdict <- if (is.na(x@adequate)) {
    "indeterminate"
  } else if (x@adequate) {
    "adequate"
  } else {
    "INADEQUATE"
  }
  is_t <- identical(x@obs_family, "student_t")
  scale_lab <- if (is_t) "scale (NIS, t-PIT)     " else "scale (NIS chi-squared) "
  fam_lab <- if (is_t) "family-fit (t-PIT JB)  " else "normality (Jarque-Bera) "
  cat(sprintf(
    "<innov_diag>: state-space model adequacy -- %s%s\n", verdict,
    if (is_t) " (Student-t obs)" else ""
  ))
  cat(sprintf(
    "  whiteness (Ljung-Box)   : %s\n", .format_p(x@whiteness_p)
  ))
  cat(sprintf(
    "  %s: %s\n", scale_lab, .format_p(x@calibration_p)
  ))
  cat(sprintf(
    "  %s: %s\n", fam_lab, .format_p(x@normality_p)
  ))
  if (!is.na(x@p_value)) {
    cat(sprintf(
      "  combined                : p = %.3g (alpha = %.3g)\n",
      x@p_value, x@alpha
    ))
  } else {
    cat(sprintf("  combined                : %s\n", x@reason))
  }
  invisible(x)
}

#' Format a diagnostic p-value for printing
#'
#' @param p A numeric p-value or `NA`.
#'
#' @returns A character scalar.
#' @noRd
#' @keywords internal
.format_p <- function(p) {
  if (length(p) != 1L || is.na(p)) {
    "not assessed"
  } else {
    sprintf("p = %.3g", p)
  }
}
