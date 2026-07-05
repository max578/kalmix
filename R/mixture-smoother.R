# mixture-smoother.R -- the Kim smoother for regime-switching models.
#
# mixture_smoother() runs the GPB1 forward pass of mixture_filter() and the
# backward recursion of Kim (1994): the smoothed regime probabilities follow
# from the filtered and chain-prior probabilities, and the smoothed states
# from per-regime-pair Rauch-Tung-Striebel updates collapsed under the
# smoothed pair weights. Like the forward pass, each backward step collapses
# a mixture to one Gaussian per regime, so the smoother is the standard
# approximation for this model class, not an exact mixture smoother (the
# exact one is exponential in the series length).

# The regime_smooth class ---------------------------------------------------

#' A regime-switching smoother pass
#'
#' An S7 object holding the output of [mixture_smoother()]: the smoothed
#' regime probabilities and the collapsed smoothed state moments, each
#' conditioning on the whole observed series.
#'
#' @param smoothed_prob An `n` by `k` matrix of smoothed regime
#'   probabilities, row `t` summing to one.
#' @param smoothed_mean An `n` by `m` matrix of collapsed smoothed state
#'   means.
#' @param smoothed_cov A list of `n` collapsed smoothed state covariances.
#' @param log_lik Numeric scalar — the mixture model log-likelihood, carried
#'   through from the forward pass.
#' @param n_regimes Integer scalar — the number of regimes `k`.
#'
#' @returns An S7 object of class `regime_smooth`.
#' @family state-space
#' @export
regime_smooth <- S7::new_class(
  name = "regime_smooth",
  package = "kalmix",
  properties = list(
    smoothed_prob = S7::class_double,
    smoothed_mean = S7::class_double,
    smoothed_cov = S7::class_list,
    log_lik = S7::class_double,
    n_regimes = S7::class_integer
  ),
  validator = function(self) {
    if (self@n_regimes < 1L) {
      return("`n_regimes` must be a positive integer")
    }
    NULL
  }
)

# The smoothing verb --------------------------------------------------------

#' Regime-switching (Kim) smoother
#'
#' Smooths a regime-switching series: runs the Gaussian pseudo-Bayesian
#' (GPB1) forward pass of [mixture_filter()] and the backward recursion of
#' Kim (1994), returning the regime probabilities and collapsed state
#' moments conditioned on the whole series. Smoothed regime probabilities
#' use the future as well as the past, so a regime that only becomes
#' apparent in hindsight is attributed to the steps where it held.
#'
#' Each backward step collapses a Gaussian mixture to one Gaussian per
#' regime under the smoothed pair weights, mirroring the forward collapse,
#' so the smoother is the standard approximation for Markov-switching
#' state-space models rather than the (exponentially costly) exact mixture
#' smoother. Missing observations are handled by the forward pass as
#' prediction-only steps and flow through the backward pass unchanged.
#'
#' @param models A list of two or more [ssm] objects, one per regime,
#'   sharing a common state and observation dimension.
#' @param y A numeric vector or `n` by `d` matrix of observations (fully
#'   `NA` rows allowed).
#' @param transition A `k` by `k` regime-transition probability matrix, or
#'   `NULL` for the sticky default of [mixture_filter()]. Defaults to
#'   `NULL`.
#' @param init_prob A length-`k` vector of initial regime probabilities, or
#'   `NULL` for uniform. Defaults to `NULL`.
#'
#' @returns A [regime_smooth].
#' @family state-space
#' @seealso [mixture_filter()] for the forward pass; [rts_smoother()] for
#'   the single-regime smoother.
#' @references
#' Kim, C.-J. (1994). Dynamic linear models with Markov-switching.
#' *Journal of Econometrics*, 60(1--2), 1--22.
#' @export
#' @examples
#' calm <- ssm(
#'   transition = 1, observation = 1, state_cov = 0.01, obs_cov = 0.25,
#'   init_state = 0, init_cov = 1
#' )
#' wild <- ssm(
#'   transition = 1, observation = 1, state_cov = 0.01, obs_cov = 4,
#'   init_state = 0, init_cov = 1
#' )
#' set.seed(1)
#' y <- c(rnorm(60, sd = 0.5), rnorm(60, sd = 2))
#' sm <- mixture_smoother(list(calm, wild), y)
#' tail(round(sm@smoothed_prob, 2))
mixture_smoother <- function(models,
                             y,
                             transition = NULL,
                             init_prob = NULL) {
  .check_regime_models(models)
  k <- length(models)
  y <- .check_observations(y, models[[1L]])
  transition <- .check_regime_transition(transition, k)
  init_prob <- .check_regime_init(init_prob, k)

  fwd <- .mixture_forward(models, y, transition, init_prob)
  .kim_backward(models, fwd, transition)
}

#' The Kim (1994) backward pass over a stored GPB1 forward pass
#'
#' Smoothed regime probabilities follow the classical recursion
#' \eqn{p_{t|n}(j) = p_{t|t}(j) \sum_i T_{ji}\, p_{t+1|n}(i) /
#' p_{t+1|t}(i)}; smoothed states run a per-pair (j, i)
#' Rauch-Tung-Striebel update from the per-regime filtered moments and are
#' collapsed under the smoothed pair weights, first over the future regime
#' `i`, then over `j` for the reported collapsed moments.
#'
#' @param models The validated list of regime [ssm] objects.
#' @param fwd The list returned by `.mixture_forward()`.
#' @param transition The `k` by `k` regime-transition matrix.
#'
#' @returns A [regime_smooth].
#' @noRd
#' @keywords internal
.kim_backward <- function(models, fwd, transition) {
  k <- fwd$n_regimes
  n <- nrow(fwd$regime_prob)
  m <- ncol(fwd$state_mean)

  smoothed_prob <- fwd$regime_prob
  smoothed_mean <- fwd$state_mean
  smoothed_cov <- fwd$state_cov

  # Per-regime smoothed moments, initialised at the last filtered step.
  x_s <- fwd$regime_means[[n]]
  p_s <- fwd$regime_covs[[n]]

  for (t in rev(seq_len(n - 1L))) {
    # Chain-prior probabilities of the next step given data to t.
    p_next_pred <- as.numeric(fwd$regime_prob[t, ] %*% transition)

    # Smoothed pair weights w[j, i] = P(S_t = j, S_t+1 = i | y_1:n).
    w <- matrix(0, k, k)
    for (j in seq_len(k)) {
      for (i in seq_len(k)) {
        w[j, i] <- fwd$regime_prob[t, j] * transition[j, i] *
          smoothed_prob[t + 1L, i] /
          max(p_next_pred[i], .Machine$double.xmin)
      }
    }
    smoothed_prob[t, ] <- rowSums(w)

    # Per-pair RTS updates from the per-regime filtered moments, collapsed
    # over the future regime i under the pair weights.
    x_new <- vector("list", k)
    p_new <- vector("list", k)
    for (j in seq_len(k)) {
      wj <- w[j, ]
      wj <- if (sum(wj) > 0) wj / sum(wj) else rep(1 / k, k)
      xj <- fwd$regime_means[[t]][[j]]
      pj <- fwd$regime_covs[[t]][[j]]

      xs_i <- vector("list", k)
      ps_i <- vector("list", k)
      for (i in seq_len(k)) {
        a_i <- .at_step(models[[i]]@transition, t + 1L)
        q_i <- .at_step(models[[i]]@state_cov, t + 1L)
        c_i <- .at_step(models[[i]]@state_intercept, t + 1L)
        x_pred <- as.numeric(a_i %*% xj) + c_i
        p_pred <- a_i %*% pj %*% t(a_i) + q_i
        gain <- pj %*% t(a_i) %*% chol2inv(
          .safe_chol(p_pred, "predicted covariance")
        )
        xs_i[[i]] <- xj + as.numeric(gain %*% (x_s[[i]] - x_pred))
        ps_i[[i]] <- .symmetrise(
          pj + gain %*% (p_s[[i]] - p_pred) %*% t(gain)
        )
      }
      x_new[[j]] <- Reduce(`+`, Map(function(p, mu) p * mu, wj, xs_i))
      p_new[[j]] <- .symmetrise(Reduce(`+`, Map(function(p, mu, sig) {
        gap <- mu - x_new[[j]]
        p * (sig + outer(gap, gap))
      }, wj, xs_i, ps_i)))
    }
    x_s <- x_new
    p_s <- p_new

    # Collapse over the current regime for the reported moments.
    pr <- smoothed_prob[t, ]
    x_mix <- Reduce(`+`, Map(function(p, mu) p * mu, pr, x_s))
    p_mix <- Reduce(`+`, Map(function(p, mu, sig) {
      gap <- mu - x_mix
      p * (sig + outer(gap, gap))
    }, pr, x_s, p_s))
    smoothed_mean[t, ] <- x_mix
    smoothed_cov[[t]] <- .symmetrise(p_mix)
  }

  regime_smooth(
    smoothed_prob = smoothed_prob,
    smoothed_mean = smoothed_mean,
    smoothed_cov = smoothed_cov,
    log_lik = fwd$log_lik,
    n_regimes = as.integer(k)
  )
}

#' @export
S7::method(print, regime_smooth) <- function(x, ...) {
  n <- nrow(x@smoothed_prob)
  cat(sprintf(
    "<regime_smooth>: %d step%s, %d regimes, log-likelihood %.3f\n",
    n, if (n == 1L) "" else "s", x@n_regimes, x@log_lik
  ))
  invisible(x)
}
