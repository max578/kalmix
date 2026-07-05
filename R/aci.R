# aci.R -- the time-resolved Assimilative Causal Inference (ACI) read-out.
#
# Assimilative Causal Inference (Andreou, Chen and Bollt 2026) casts causality
# as a Bayesian data-assimilation inverse problem: a hidden state is a cause of
# an observed effect at time t exactly when incorporating the future of the
# effect sharpens the estimate of the state, that is when the smoothing
# distribution of the state departs from its filtering distribution. The causal
# information a measurement carries about the hidden state at each step is the
# relative entropy (Kullback-Leibler divergence) between the smoother and the
# filter state distributions (their eq. 7); time-aggregated, the divergence
# profile against an expanding future window yields a threshold-free causal
# information rate and a decision lead-time (their eq. 9).
#
# The divergence is a standard closed form for two Gaussians and is
# implemented natively here in a handful of lines; when kernR is installed its
# `relative_entropy()` is used as an independent oracle in the test suite (it
# is not a runtime dependency). aci() runs the Kalman filter and the
# Rauch-Tung-Striebel smoother of kalman.R and reads the causal information
# series off the two passes; causal_information_rate() (cir.R) forms the
# objective CIR and lead-time from an expanding-window divergence profile.

# The aci_fit class -----------------------------------------------------------

#' An Assimilative Causal Inference read-out
#'
#' An S7 object holding the output of [aci()]: the per-step causal information
#' series, its time aggregate (the mean causal information), the objective
#' causal information rate and the implied decision lead-time, alongside the
#' filter and smoother passes the read-out was computed from.
#'
#' The causal information at step `t` is the relative entropy of the smoothing
#' distribution from the filtering distribution of the latent state,
#' \deqn{I_t = D_{\mathrm{KL}}\!\left(
#'   N(x_{t\mid n}, P_{t\mid n}) \,\Vert\, N(x_{t\mid t}, P_{t\mid t}) \right),}
#' the information the future of the series adds to the present estimate of the
#' state. It is non-negative and zero only where the future is uninformative
#' about the state at `t`.
#'
#' @param causal_information A numeric vector of length `n` -- the per-step
#'   causal information series \eqn{I_t}, the smoother-from-filter relative
#'   entropy at each step.
#' @param mean_causal_information Numeric scalar -- the time-average of
#'   `causal_information`, a single summary of how much the future informs the
#'   state over the series.
#' @param causal_information_rate Numeric scalar -- the objective,
#'   threshold-free causal information rate (a lead-time in the time units of
#'   the series), or `NA_real_` if it was not requested.
#' @param lead_time Numeric scalar -- the decision lead-time, an alias of the
#'   causal information rate in time units, or `NA_real_` if not requested.
#' @param dt Numeric scalar -- the sampling interval the lead-time is expressed
#'   in. Defaults to `1`, so the lead-time is in steps.
#' @param filter The forward pass the read-out was computed from: a
#'   [kalman_fit], or a [regime_fit] for a regime-switching model.
#' @param smoother The smoothing pass: an [rts_fit], or a [regime_smooth]
#'   for a regime-switching model.
#' @param grounding Character scalar -- the grounding token of the causal
#'   verdict, `"grounded"` or `"[unverified]"`; see the token definitions in
#'   [innovation_diagnostics()].
#' @param grounding_reason Character scalar -- a short label for why the verdict
#'   carries that token (`"grounded"`, `"mechanism_unverified"`,
#'   `"model_inadequate"` or `"insufficient_data"`).
#' @param adequacy The [innov_diag] model-adequacy diagnostics the grounding was
#'   read from.
#'
#' @returns An S7 object of class `aci_fit`.
#' @family causal
#' @export
aci_fit <- S7::new_class(
  name = "aci_fit",
  package = "kalmix",
  properties = list(
    causal_information = S7::class_double,
    mean_causal_information = S7::class_double,
    causal_information_rate = S7::new_property(
      class = S7::class_double,
      default = NA_real_
    ),
    lead_time = S7::new_property(
      class = S7::class_double,
      default = NA_real_
    ),
    dt = S7::new_property(class = S7::class_double, default = 1),
    filter = S7::class_any,
    smoother = S7::class_any,
    grounding = S7::new_property(
      class = S7::class_character,
      default = "[unverified]"
    ),
    grounding_reason = S7::new_property(
      class = S7::class_character,
      default = NA_character_
    ),
    adequacy = S7::class_any
  ),
  validator = function(self) {
    if (anyNA(self@causal_information) ||
        any(self@causal_information < -1e-8)) {
      return("`causal_information` must be non-negative and free of NA")
    }
    if (length(self@mean_causal_information) != 1L) {
      return("`mean_causal_information` must be a scalar")
    }
    if (length(self@dt) != 1L || !is.finite(self@dt) || self@dt <= 0) {
      return("`dt` must be a finite positive scalar")
    }
    NULL
  }
)

# The causal metric -----------------------------------------------------------

#' Relative entropy between two Gaussians (the ACI causal metric)
#'
#' Computes the Kullback-Leibler divergence
#' \eqn{D_{\mathrm{KL}}(p \Vert q)} between two `k`-variate Gaussians
#' \eqn{p = N(\mu_p, \Sigma_p)} and \eqn{q = N(\mu_q, \Sigma_q)} in closed form,
#' \deqn{D_{\mathrm{KL}}(p \Vert q) = \tfrac12 \left[
#'   \operatorname{tr}(\Sigma_q^{-1}\Sigma_p)
#'   + (\mu_q - \mu_p)^\top \Sigma_q^{-1} (\mu_q - \mu_p) - k
#'   + \log\frac{\lvert\Sigma_q\rvert}{\lvert\Sigma_p\rvert} \right].}
#' This is the operational statistic of Assimilative Causal Inference (Andreou,
#' Chen and Bollt 2026): with `p` the smoother posterior of the latent state and
#' `q` the filter posterior, a non-zero value identifies the future of the
#' observed series as informative about -- a cause of -- the state at that step.
#'
#' The divergence is non-negative, zero exactly when the two Gaussians coincide,
#' and asymmetric in its arguments. This is canonical Gaussian arithmetic, not a
#' reimplementation of a divergence API: the package ships it natively so the
#' read-out is self-contained, and cross-checks the value against
#' `kernR::relative_entropy()` in its test suite when kernR is installed.
#'
#' @param mu_p Numeric vector -- the mean of `p`, the integrating density (the
#'   smoother posterior in eq. 7, the complete smoother in eq. 8).
#' @param sigma_p A `k` by `k` covariance matrix of `p`.
#' @param mu_q Numeric vector -- the mean of `q`, the reference density (the
#'   filter posterior in eq. 7, the lagged estimate in eq. 8).
#' @param sigma_q A `k` by `k` covariance matrix of `q`, positive definite (it
#'   is inverted through its Cholesky factor).
#'
#' @returns Numeric scalar -- the relative entropy, floored at zero to absorb
#'   rounding noise.
#' @references
#' Kullback, S. and Leibler, R. A. (1951). On information and sufficiency.
#' *Annals of Mathematical Statistics*, 22(1), 79--86.
#'
#' Andreou, M., Chen, N. and Bollt, E. (2026). Assimilative causal inference.
#' *Nature Communications*, 17, 1854.
#' @noRd
#' @keywords internal
.gaussian_relative_entropy <- function(mu_p, sigma_p, mu_q, sigma_q) {
  mu_p <- as.numeric(mu_p)
  mu_q <- as.numeric(mu_q)
  k <- length(mu_p)
  sigma_p <- as.matrix(sigma_p)
  sigma_q <- as.matrix(sigma_q)

  q_chol <- .safe_chol(sigma_q, "reference covariance")
  q_inv <- chol2inv(q_chol)
  delta <- mu_q - mu_p

  tr_term <- sum(q_inv * sigma_p)
  quad <- sum(delta * as.numeric(q_inv %*% delta))
  log_det_q <- 2 * sum(log(diag(q_chol)))
  log_det_p <- as.numeric(determinant(sigma_p, logarithm = TRUE)$modulus)

  val <- 0.5 * (tr_term + quad - k + (log_det_q - log_det_p))
  max(val, 0)
}

# The read-out verb -----------------------------------------------------------

#' Assimilative Causal Inference read-out for a state-space model
#'
#' Runs the closed-form Kalman filter and Rauch-Tung-Striebel smoother of
#' [kalman_filter()] and [rts_smoother()] over an [ssm], then reads the
#' time-resolved causal information series off the two passes following Andreou,
#' Chen and Bollt (2026). The causal information at each step is the relative
#' entropy of the smoothing distribution of the latent state from its filtering
#' distribution,
#' \deqn{I_t = D_{\mathrm{KL}}\!\left(
#'   N(x_{t\mid n}, P_{t\mid n}) \,\Vert\, N(x_{t\mid t}, P_{t\mid t}) \right),}
#' the information the future of the series adds to the present estimate of the
#' state (their eq. 7). The series is summarised by its time-average, the mean
#' causal information.
#'
#' By default the objective causal information rate and the decision lead-time
#' are also computed (`lead_time = TRUE`), through [causal_information_rate()]:
#' the divergence of the complete smoother from an expanding-future-window
#' smoother is integrated into a single threshold-free lead-time in the time
#' units of the series (their eqs. 8--9). Set `lead_time = FALSE` to return the
#' causal information series alone, which is cheaper for long series.
#'
#' Under a Student-t observation family (`obs_family = "student_t"`) the
#' filter pass is the outlier-robust variational update, whose per-step
#' posteriors are Gaussian approximations. The relative-entropy chain treats
#' those approximations as exact Gaussians, so the causal information and
#' lead-time inherit the variational approximation; the Gaussian family
#' involves no approximation.
#'
#' A regime-switching model is supplied as a list of [ssm] regimes (with the
#' optional `transition` and `init_prob` of [mixture_filter()]): the read-out
#' then compares the Kim smoother's collapsed posteriors with the GPB1
#' filter's. Both are the model class's own Gaussian collapses of the true
#' mixtures, so the causal information inherits that approximation; the
#' verdict is grounded through the mixture-PIT diagnostics, and the
#' lead-time uses the expanding-window construction (the online recursion is
#' linear-Gaussian only).
#'
#' The recovered state is interpretable as a cause of the observed series only
#' under the maintained assumption that the [ssm] is the data-generating
#' mechanism. That assumption is not taken on trust: the read-out is grounded
#' through [innovation_diagnostics()], and a model whose standardised
#' innovations are not white, not correctly scaled, or not Gaussian is declared
#' inadequate, so the verdict abstains with a `[unverified]` grounding token.
#' Passing the diagnostics establishes self-consistency only; the verdict is
#' labelled `"grounded"` solely when an adequate model is accompanied by
#' declared, dated `mechanism` provenance, and is `"[unverified]"` otherwise.
#'
#' @param model An [ssm], or a list of two or more [ssm] regimes for a
#'   regime-switching read-out.
#' @param y A numeric vector or `n` by `d` matrix of observations.
#' @param lead_time Logical scalar -- whether to compute the objective causal
#'   information rate and lead-time as well as the causal information series.
#'   Defaults to `TRUE`.
#' @param transition,init_prob The regime-transition matrix and initial
#'   regime distribution when `model` is a list of regimes (see
#'   [mixture_filter()]); ignored for an [ssm]. Default to `NULL`.
#' @param dt Numeric scalar -- the sampling interval, so the lead-time is in
#'   meaningful time units (for example years for an annual ENSO index).
#'   Defaults to `1`, giving the lead-time in steps.
#' @param eval_points An optional integer vector of step indices at which to
#'   anchor the lead-time computation; passed to [causal_information_rate()].
#'   Defaults to `NULL`, which spreads a handful of anchors across the interior
#'   of the series.
#' @param max_lag Integer scalar -- the largest future-window lag, in steps,
#'   used for the lead-time. Defaults to one fifth of the series length.
#' @param mechanism Optional declared provenance for the model as the
#'   data-generating mechanism: `NULL` (the default, an unverified
#'   mechanism), or a list carrying a `verified_on` `Date` on which the model
#'   was checked against the external authority that owns it. Only a verified
#'   mechanism can raise an adequate read-out to a `"grounded"` verdict.
#' @param alpha Numeric scalar -- the significance level the model-adequacy
#'   diagnostics abstain at. Defaults to `0.05`.
#'
#' @returns An [aci_fit].
#' @family causal
#' @seealso [causal_information_rate()] for the lead-time alone;
#'   [kalman_filter()] and [rts_smoother()] for the underlying passes.
#' @references
#' Andreou, M., Chen, N. and Bollt, E. (2026). Assimilative causal inference.
#' *Nature Communications*, 17, 1854.
#'
#' Kalman, R. E. (1960). A new approach to linear filtering and prediction
#' problems. *Journal of Basic Engineering*, 82(1), 35--45.
#'
#' Rauch, H. E., Tung, F. and Striebel, C. T. (1965). Maximum likelihood
#' estimates of linear dynamic systems. *AIAA Journal*, 3(8), 1445--1450.
#' @export
#' @examples
#' ## A latent oscillation seen through a slow, noisy observation channel: the
#' ## future of the series is informative about the hidden state, so the causal
#' ## information is positive.
#' model <- ssm(
#'   transition = matrix(c(0.6, -0.5, 0.5, 0.6), nrow = 2),
#'   observation = matrix(c(1, 0), nrow = 1),
#'   state_cov = diag(c(0.2, 0.2)),
#'   obs_cov = 1,
#'   init_state = c(0, 0),
#'   init_cov = diag(c(5, 5))
#' )
#' set.seed(1)
#' x <- numeric(200)
#' s <- c(0, 0)
#' for (t in seq_len(200)) {
#'   s <- as.numeric(model@transition[[1]] %*% s) + rnorm(2, sd = 0.45)
#'   x[t] <- s[1] + rnorm(1)
#' }
#' fit <- aci(model, x, lead_time = FALSE)
#' fit@mean_causal_information
aci <- function(model,
                y,
                lead_time = TRUE,
                dt = 1,
                eval_points = NULL,
                max_lag = NULL,
                mechanism = NULL,
                alpha = 0.05,
                transition = NULL,
                init_prob = NULL) {
  is_regime <- is.list(model) && !S7::S7_inherits(model, ssm)
  if (!is_regime && !S7::S7_inherits(model, ssm)) {
    cli::cli_abort(
      "`model` must be an {.cls ssm} or a list of {.cls ssm} regimes."
    )
  }
  if (length(dt) != 1L || !is.finite(dt) || dt <= 0) {
    cli::cli_abort("`dt` must be a finite positive scalar.")
  }

  if (is_regime) {
    filter <- mixture_filter(model, y, transition = transition,
                             init_prob = init_prob)
    smoother <- mixture_smoother(model, y, transition = transition,
                                 init_prob = init_prob)
    info <- vapply(
      seq_len(nrow(filter@state_mean)),
      function(t) {
        .gaussian_relative_entropy(
          smoother@smoothed_mean[t, ], smoother@smoothed_cov[[t]],
          filter@state_mean[t, ], filter@state_cov[[t]]
        )
      },
      numeric(1L)
    )
  } else {
    filter <- kalman_filter(model, y)
    smoother <- rts_smoother(filter)
    info <- .causal_information_series(filter, smoother)
  }

  rate <- NA_real_
  lead <- NA_real_
  if (isTRUE(lead_time)) {
    rate <- causal_information_rate(
      model, y,
      dt = dt, eval_points = eval_points, max_lag = max_lag,
      transition = transition, init_prob = init_prob
    )
    lead <- rate
  }

  adequacy <- innovation_diagnostics(filter, alpha = alpha)
  grounding <- .aci_grounding(adequacy, mechanism)

  aci_fit(
    causal_information = info,
    mean_causal_information = mean(info),
    causal_information_rate = rate,
    lead_time = lead,
    dt = dt,
    filter = filter,
    smoother = smoother,
    grounding = grounding$grounding,
    grounding_reason = grounding$reason,
    adequacy = adequacy
  )
}

# The per-step series ---------------------------------------------------------

#' The per-step causal information series
#'
#' Reads the smoother-from-filter relative entropy off a matched [kalman_fit]
#' and [rts_fit] at every step, the time-resolved causal information of Andreou,
#' Chen and Bollt (2026, eq. 7). Factored out so [aci()] and the lead-time
#' anchors share one definition.
#'
#' @param filter A [kalman_fit].
#' @param smoother The [rts_fit] from the same model and series.
#'
#' @returns A numeric vector of length `n`.
#' @noRd
#' @keywords internal
.causal_information_series <- function(filter, smoother) {
  n <- nrow(filter@filtered_mean)
  vapply(
    seq_len(n),
    function(t) {
      .gaussian_relative_entropy(
        smoother@smoothed_mean[t, ], smoother@smoothed_cov[[t]],
        filter@filtered_mean[t, ], filter@filtered_cov[[t]]
      )
    },
    numeric(1L)
  )
}

#' @export
S7::method(print, aci_fit) <- function(x, ...) {
  n <- length(x@causal_information)
  cat(sprintf(
    "<aci_fit>: assimilative causal inference (%d step%s)\n",
    n, if (n == 1L) "" else "s"
  ))
  cat(sprintf(
    "  mean causal information : %.4f nats\n",
    x@mean_causal_information
  ))
  if (is.finite(x@lead_time)) {
    unit <- if (isTRUE(all.equal(x@dt, 1))) "steps" else "time units"
    cat(sprintf(
      "  decision lead-time      : %.3f %s\n",
      x@lead_time, unit
    ))
  } else {
    cat("  decision lead-time      : not computed\n")
  }
  reason <- if (is.na(x@grounding_reason)) {
    ""
  } else {
    sprintf(" (%s)", x@grounding_reason)
  }
  cat(sprintf(
    "  grounding               : %s%s\n", x@grounding, reason
  ))
  invisible(x)
}
