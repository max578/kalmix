# ssm-fit.R -- maximum-likelihood estimation of ssm noise parameters.
#
# ssm_fit() estimates the process and observation covariances (and optionally
# the transition matrix) of a linear-Gaussian or Student-t ssm by maximising
# the prediction-error log-likelihood of kalman_filter(). Covariances are
# parameterised through their log-Cholesky factors, so the optimiser works on
# an unconstrained scale and every iterate is positive definite; the fit is
# repeated from jittered starts so a poor template does not strand the search
# in a local optimum. Standard errors come from the numerical Hessian at the
# optimum, mapped back to the covariance entries by a numerical delta method.

# The ssm_mle class -------------------------------------------------------

#' A maximum-likelihood ssm fit
#'
#' An S7 object holding the output of [ssm_fit()]: the fitted [ssm] (ready
#' for every verb that accepts one), the maximised log-likelihood, the
#' estimates and standard errors of the estimated entries, and the optimiser
#' diagnostics.
#'
#' @param model The fitted [ssm].
#' @param log_lik Numeric scalar — the maximised prediction-error
#'   log-likelihood.
#' @param estimates A named list of the estimated matrices (`state_cov`,
#'   `obs_cov` and, when requested, `transition`).
#' @param std_errors A named list matching `estimates`, each entry a matrix of
#'   delta-method standard errors for the corresponding estimates.
#' @param convergence Integer scalar — the [stats::optim()] convergence code
#'   of the winning start (`0` means converged).
#' @param n_start Integer scalar — the number of optimisation starts run.
#' @param start_log_lik A numeric vector with the best log-likelihood reached
#'   from each start.
#'
#' @returns An S7 object of class `ssm_mle`.
#' @family state-space
#' @export
ssm_mle <- S7::new_class(
  name = "ssm_mle",
  package = "kalmix",
  properties = list(
    model = ssm,
    log_lik = S7::class_double,
    estimates = S7::class_list,
    std_errors = S7::class_list,
    convergence = S7::class_integer,
    n_start = S7::class_integer,
    start_log_lik = S7::class_double
  )
)

# The fitting verb --------------------------------------------------------

#' Fit ssm noise parameters by maximum likelihood
#'
#' Estimates the process covariance, the observation covariance and
#' (optionally) the transition matrix of an [ssm] by maximising the
#' prediction-error log-likelihood of [kalman_filter()] with [stats::optim()].
#' The supplied model is the template: its non-estimated components, its prior
#' and its observation family are kept as given, and its current values seed
#' the first optimisation start.
#'
#' Covariances are optimised through their log-Cholesky factors (the
#' logarithm of the factor's diagonal, the off-diagonal entries free), so the
#' search is unconstrained and every iterate is a valid covariance. The fit
#' is repeated from `n_start` starts (the template plus a deterministic
#' ladder of log-scale shifts, so no seed is involved) and the best optimum
#' is kept, guarding against local optima at a modest cost.
#' Standard errors are read from the numerical Hessian at the optimum and
#' mapped to the covariance entries by a numerical delta method; they are
#' asymptotic and can be optimistic on short series, where the likelihood
#' separates the two noise sources weakly (a wide interval there is a finding
#' about the data, not a defect of the fit).
#'
#' Missing observations (fully `NA` rows) are handled by the filter as
#' prediction-only steps, so a gappy series can be fitted directly. Under a
#' Student-t observation family the same machinery maximises the robust
#' filter's likelihood at the template's fixed degrees of freedom; combine
#' with [estimate_obs_df()] to choose the degrees of freedom.
#'
#' @param model An [ssm] — the template supplying structure, prior, family
#'   and starting values. Estimated components must be time-invariant.
#' @param y A numeric vector or `n` by `d` matrix of observations (`NA` rows
#'   allowed).
#' @param estimate A character vector naming the components to estimate, from
#'   `"state_cov"` and `"obs_cov"`. Defaults to both.
#' @param transition Logical scalar — whether to estimate the transition
#'   matrix as well (unconstrained; a fit whose spectral radius reaches one
#'   is reported with a warning). Defaults to `FALSE`.
#' @param n_start Integer scalar — the number of optimisation starts.
#'   Defaults to `3`.
#' @param control A list passed to [stats::optim()]'s `control` argument.
#'
#' @returns An [ssm_mle].
#' @family state-space
#' @seealso [kalman_filter()] for the likelihood; [estimate_obs_df()] for the
#'   Student-t degrees of freedom; [ssm()] for the template.
#' @references
#' Harvey, A. C. (1989). *Forecasting, Structural Time Series Models and the
#' Kalman Filter*. Cambridge University Press. (Prediction-error
#' decomposition maximum likelihood.)
#'
#' Durbin, J. and Koopman, S. J. (2012). *Time Series Analysis by State Space
#' Methods*. 2nd ed. Oxford University Press.
#' @export
#' @examples
#' ## Recover the noise split of a local-level series. One start keeps the
#' ## example fast; the default multi-start is the safer everyday setting.
#' set.seed(1)
#' level <- cumsum(rnorm(150, sd = 0.5))
#' y <- level + rnorm(150)
#' template <- ssm(
#'   transition = 1, observation = 1, state_cov = 1, obs_cov = 1,
#'   init_state = y[1], init_cov = 10
#' )
#' fit <- ssm_fit(template, y, n_start = 1L)
#' fit@estimates$state_cov
#' fit@estimates$obs_cov
ssm_fit <- function(model,
                    y,
                    estimate = c("state_cov", "obs_cov"),
                    transition = FALSE,
                    n_start = 3L,
                    control = list()) {
  if (!S7::S7_inherits(model, ssm)) {
    cli::cli_abort("`model` must be an {.cls ssm}.")
  }
  estimate <- match.arg(estimate, c("state_cov", "obs_cov"),
                        several.ok = TRUE)
  if (!is.logical(transition) || length(transition) != 1L ||
      is.na(transition)) {
    cli::cli_abort("`transition` must be `TRUE` or `FALSE`.")
  }
  n_start <- as.integer(n_start)
  if (length(n_start) != 1L || is.na(n_start) || n_start < 1L) {
    cli::cli_abort("`n_start` must be a positive integer scalar.")
  }
  y <- .check_observations(y, model)

  parts <- estimate
  if (transition) {
    parts <- c(parts, "transition")
  }
  for (nm in parts) {
    if (length(S7::prop(model, nm)) != 1L) {
      cli::cli_abort(
        "Estimated component `{nm}` must be time-invariant in the template."
      )
    }
  }

  pack <- .ssm_fit_pack(model, parts)
  objective <- function(theta) {
    fitted <- .ssm_fit_apply(model, parts, theta, pack)
    ll <- tryCatch(
      kalman_filter(fitted, y)@log_lik,
      error = function(e) -Inf
    )
    if (!is.finite(ll)) {
      return(1e10)
    }
    -ll
  }

  # Multi-start: the template's own values first, then a deterministic ladder
  # of shifted copies (variances scaled up and down on the log scale), so the
  # fit is reproducible without a seed.
  starts <- vector("list", n_start)
  starts[[1L]] <- pack$theta
  if (n_start > 1L) {
    for (i in 2L:n_start) {
      shift <- if (i %% 2L == 0L) 1.5 * (i %/% 2L) else -1.5 * (i %/% 2L)
      starts[[i]] <- pack$theta + shift
    }
  }

  best <- NULL
  start_ll <- rep(NA_real_, n_start)
  for (i in seq_len(n_start)) {
    opt <- tryCatch(
      stats::optim(starts[[i]], objective, method = "BFGS",
                   control = control),
      error = function(e) NULL
    )
    if (is.null(opt)) {
      next
    }
    start_ll[i] <- -opt$value
    if (is.null(best) || opt$value < best$value) {
      best <- opt
    }
  }
  if (is.null(best)) {
    cli::cli_abort("Every optimisation start failed; check the template.")
  }

  fitted <- .ssm_fit_apply(model, parts, best$par, pack)
  if (transition) {
    a_hat <- fitted@transition[[1L]]
    if (max(Mod(eigen(a_hat, only.values = TRUE)$values)) >= 1) {
      cli::cli_warn(
        "The fitted transition matrix is not stable (spectral radius >= 1)."
      )
    }
  }

  # Delta-method standard errors: numerical Hessian on the unconstrained
  # scale, numerical Jacobian of each estimated entry with respect to theta.
  hess <- tryCatch(
    stats::optimHess(best$par, objective, control = control),
    error = function(e) NULL
  )
  vcov_theta <- if (!is.null(hess)) {
    tryCatch(solve(hess), error = function(e) NULL)
  } else {
    NULL
  }
  std_errors <- .ssm_fit_std_errors(model, parts, best$par, pack, vcov_theta)

  estimates <- lapply(parts, function(nm) S7::prop(fitted, nm)[[1L]])
  names(estimates) <- parts

  ssm_mle(
    model = fitted,
    log_lik = -best$value,
    estimates = estimates,
    std_errors = std_errors,
    convergence = as.integer(best$convergence),
    n_start = n_start,
    start_log_lik = start_ll
  )
}

#' @export
S7::method(print, ssm_mle) <- function(x, ...) {
  cat(sprintf(
    "<ssm_mle>: log-likelihood %.3f, %d start%s, convergence %d\n",
    x@log_lik, x@n_start, if (x@n_start == 1L) "" else "s", x@convergence
  ))
  for (nm in names(x@estimates)) {
    cat(sprintf("  %s:\n", nm))
    print(round(x@estimates[[nm]], 4))
  }
  invisible(x)
}

# Parameterisation helpers ------------------------------------------------

#' Pack the estimated components of a template into one parameter vector
#'
#' Covariances enter through their log-Cholesky factor (log diagonal, free
#' lower triangle of the transposed upper factor); a transition matrix enters
#' unconstrained. Returns the packed vector plus the index map used to unpack.
#'
#' @param model The template [ssm].
#' @param parts Character vector of component names.
#'
#' @returns A list with `theta` and per-part index/shape metadata.
#' @noRd
#' @keywords internal
.ssm_fit_pack <- function(model, parts) {
  theta <- numeric(0L)
  meta <- vector("list", length(parts))
  names(meta) <- parts
  for (nm in parts) {
    mat <- S7::prop(model, nm)[[1L]]
    k <- nrow(mat)
    if (nm == "transition") {
      vals <- as.numeric(mat)
    } else {
      l <- t(chol(mat + diag(1e-10, k)))
      vals <- c(log(diag(l)), l[lower.tri(l)])
    }
    meta[[nm]] <- list(
      offset = length(theta), length = length(vals), dim = k
    )
    theta <- c(theta, vals)
  }
  list(theta = theta, meta = meta)
}

#' Rebuild an ssm from a packed parameter vector
#'
#' @param model The template [ssm].
#' @param parts Character vector of component names.
#' @param theta The packed parameter vector.
#' @param pack The index map from [.ssm_fit_pack()].
#'
#' @returns An [ssm] with the estimated components replaced.
#' @noRd
#' @keywords internal
.ssm_fit_apply <- function(model, parts, theta, pack) {
  args <- list(
    transition = model@transition[[1L]],
    observation = model@observation[[1L]],
    state_cov = model@state_cov[[1L]],
    obs_cov = model@obs_cov[[1L]],
    state_intercept = model@state_intercept[[1L]],
    obs_intercept = model@obs_intercept[[1L]],
    init_state = model@init_state,
    init_cov = model@init_cov,
    obs_family = model@obs_family,
    obs_df = model@obs_df
  )
  # A time-varying non-estimated component is carried through as the list it
  # already is; only time-invariant components can be estimated.
  for (nm in c("transition", "observation", "state_cov", "obs_cov",
               "state_intercept", "obs_intercept")) {
    lst <- S7::prop(model, nm)
    if (length(lst) != 1L) {
      args[[nm]] <- lst
    }
  }
  for (nm in parts) {
    mi <- pack$meta[[nm]]
    vals <- theta[mi$offset + seq_len(mi$length)]
    k <- mi$dim
    if (nm == "transition") {
      args[[nm]] <- matrix(vals, k, k)
    } else {
      l <- matrix(0, k, k)
      diag(l) <- exp(vals[seq_len(k)])
      l[lower.tri(l)] <- vals[-seq_len(k)]
      args[[nm]] <- .symmetrise(l %*% t(l))
    }
  }
  do.call(ssm, args)
}

#' Delta-method standard errors for the estimated entries
#'
#' Maps the covariance of the unconstrained parameters to each estimated
#' matrix entry through a central-difference Jacobian. Returns a list of
#' SE matrices matching the estimates, or `NA` matrices when the Hessian was
#' unavailable or singular.
#'
#' @param model The template [ssm].
#' @param parts Character vector of component names.
#' @param theta The optimum on the unconstrained scale.
#' @param pack The index map from [.ssm_fit_pack()].
#' @param vcov_theta The inverse Hessian at the optimum, or `NULL`.
#'
#' @returns A named list of SE matrices.
#' @noRd
#' @keywords internal
.ssm_fit_std_errors <- function(model, parts, theta, pack, vcov_theta) {
  out <- vector("list", length(parts))
  names(out) <- parts
  for (nm in parts) {
    k <- pack$meta[[nm]]$dim
    out[[nm]] <- matrix(NA_real_, k, k)
  }
  if (is.null(vcov_theta)) {
    return(out)
  }

  entry_values <- function(th) {
    unlist(lapply(parts, function(nm) {
      mi <- pack$meta[[nm]]
      vals <- th[mi$offset + seq_len(mi$length)]
      k <- mi$dim
      if (nm == "transition") {
        vals
      } else {
        l <- matrix(0, k, k)
        diag(l) <- exp(vals[seq_len(k)])
        l[lower.tri(l)] <- vals[-seq_len(k)]
        as.numeric(l %*% t(l))
      }
    }))
  }

  p <- length(theta)
  base <- entry_values(theta)
  jac <- matrix(0, nrow = length(base), ncol = p)
  h <- pmax(1e-5, 1e-5 * abs(theta))
  for (j in seq_len(p)) {
    up <- theta
    dn <- theta
    up[j] <- up[j] + h[j]
    dn[j] <- dn[j] - h[j]
    jac[, j] <- (entry_values(up) - entry_values(dn)) / (2 * h[j])
  }
  vv <- diag(jac %*% vcov_theta %*% t(jac))
  vv[vv < 0] <- NA_real_
  se_flat <- sqrt(vv)

  pos <- 0L
  for (nm in parts) {
    k <- pack$meta[[nm]]$dim
    out[[nm]] <- matrix(se_flat[pos + seq_len(k * k)], k, k)
    pos <- pos + k * k
  }
  out
}
