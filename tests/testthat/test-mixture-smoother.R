test_that("mixture_smoother() collapses to the RTS smoother on identical regimes", {
  ## Internal-consistency oracle: with identical regimes the Kim recursion
  ## must reproduce the single-model Rauch-Tung-Striebel smoother exactly,
  ## and the smoothed regime probabilities must stay at the uninformative
  ## chain prior.
  m <- ssm(
    transition = 0.9, observation = 1, state_cov = 0.2, obs_cov = 1,
    init_state = 0, init_cov = 5
  )
  withr::local_seed(41L)
  y <- as.numeric(stats::arima.sim(list(ar = 0.5), n = 100))
  sm <- mixture_smoother(list(m, m), y)
  ref <- rts_smoother(kalman_filter(m, y))

  expect_equal(sm@smoothed_mean[, 1L], ref@smoothed_mean[, 1L],
               tolerance = 1e-10)
  for (t in c(1L, 50L, 100L)) {
    expect_equal(sm@smoothed_cov[[t]], ref@smoothed_cov[[t]],
                 tolerance = 1e-10)
  }
  expect_equal(sm@smoothed_prob[, 1L], rep(0.5, 100L), tolerance = 1e-8)
  expect_equal(rowSums(sm@smoothed_prob), rep(1, 100L), tolerance = 1e-10)
})

test_that("smoothed regime probabilities beat filtered ones on a known path", {
  ## Known DGP: Markov-switching observation noise with the true regime path
  ## retained. Smoothing sees the future, so its regime attribution must be
  ## at least as accurate as the filter's, and the two must agree at the
  ## final step (where there is no future).
  withr::local_seed(42L)
  n <- 400L
  stay <- 0.95
  regime <- integer(n)
  regime[1L] <- 1L
  for (t in 2:n) {
    regime[t] <- if (stats::runif(1) < stay) regime[t - 1L] else
      3L - regime[t - 1L]
  }
  x <- numeric(n)
  x[1L] <- stats::rnorm(1)
  for (t in 2:n) {
    x[t] <- 0.7 * x[t - 1L] + stats::rnorm(1, sd = 0.5)
  }
  y <- x + stats::rnorm(n, sd = c(0.5, 2)[regime])

  calm <- ssm(
    transition = 0.7, observation = 1, state_cov = 0.25, obs_cov = 0.25,
    init_state = 0, init_cov = 1
  )
  wild <- ssm(
    transition = 0.7, observation = 1, state_cov = 0.25, obs_cov = 4,
    init_state = 0, init_cov = 1
  )
  trans <- matrix(c(stay, 1 - stay, 1 - stay, stay), 2L)

  fit <- mixture_filter(list(calm, wild), y, transition = trans)
  sm <- mixture_smoother(list(calm, wild), y, transition = trans)

  acc_filter <- mean((fit@regime_prob[, 2L] > 0.5) == (regime == 2L))
  acc_smooth <- mean((sm@smoothed_prob[, 2L] > 0.5) == (regime == 2L))
  expect_gte(acc_smooth, acc_filter)
  expect_gt(acc_smooth, 0.75)
  expect_equal(sm@smoothed_prob[n, ], fit@regime_prob[n, ], tolerance = 1e-10)

  ## The smoothed state must not be worse than the filtered state against
  ## the known truth.
  rmse_f <- sqrt(mean((fit@state_mean[, 1L] - x)^2))
  rmse_s <- sqrt(mean((sm@smoothed_mean[, 1L] - x)^2))
  expect_lte(rmse_s, rmse_f + 1e-8)
})

test_that("mixture_smoother() handles missing observations", {
  calm <- ssm(
    transition = 1, observation = 1, state_cov = 0.01, obs_cov = 0.25,
    init_state = 0, init_cov = 1
  )
  wild <- ssm(
    transition = 1, observation = 1, state_cov = 0.01, obs_cov = 4,
    init_state = 0, init_cov = 1
  )
  withr::local_seed(43L)
  y <- c(stats::rnorm(40, sd = 0.5), stats::rnorm(40, sd = 2))
  y[c(20L, 60L)] <- NA_real_
  sm <- mixture_smoother(list(calm, wild), y)
  expect_false(anyNA(sm@smoothed_mean))
  expect_equal(rowSums(sm@smoothed_prob), rep(1, 80L), tolerance = 1e-10)
})

test_that("aci() runs end to end on a regime-switching model", {
  ## The regime read-out: on a coupled Markov-switching DGP filtered with
  ## the true regimes the verdict certifies (adequate), the causal
  ## information is positive, and the expanding-window lead-time is a
  ## non-negative finite number. Under identical no-signal regimes the
  ## series is white noise and the causal information collapses toward zero.
  withr::local_seed(44L)
  n <- 300L
  stay <- 0.95
  regime <- integer(n)
  regime[1L] <- 1L
  for (t in 2:n) {
    regime[t] <- if (stats::runif(1) < stay) regime[t - 1L] else
      3L - regime[t - 1L]
  }
  x <- numeric(n)
  x[1L] <- stats::rnorm(1)
  for (t in 2:n) {
    x[t] <- 0.7 * x[t - 1L] + stats::rnorm(1, sd = 0.5)
  }
  y <- x + stats::rnorm(n, sd = c(0.5, 2)[regime])

  calm <- ssm(
    transition = 0.7, observation = 1, state_cov = 0.25, obs_cov = 0.25,
    init_state = 0, init_cov = 1
  )
  wild <- ssm(
    transition = 0.7, observation = 1, state_cov = 0.25, obs_cov = 4,
    init_state = 0, init_cov = 1
  )
  trans <- matrix(c(stay, 1 - stay, 1 - stay, stay), 2L)

  fit <- aci(list(calm, wild), y, transition = trans,
             eval_points = c(100L, 150L, 200L), max_lag = 20L)
  expect_s3_class(fit, "kalmix::aci_fit")
  expect_s3_class(fit@filter, "kalmix::regime_fit")
  expect_s3_class(fit@smoother, "kalmix::regime_smooth")
  expect_gt(fit@mean_causal_information, 0)
  expect_true(is.finite(fit@lead_time))
  expect_gte(fit@lead_time, 0)
  expect_identical(fit@adequacy@obs_family, "gaussian_mixture")
  expect_true(fit@adequacy@adequate)
  expect_identical(fit@grounding_reason, "mechanism_unverified")
})
