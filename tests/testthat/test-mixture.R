test_that("mixture_filter() recovers a known volatility switch", {
  ## Known DGP: a calm half (sd 0.5) then a turbulent half (sd 3). The filter
  ## must assign low probability to the wild regime in the calm half and high
  ## probability in the wild half. The regime labels are the known truth.
  calm <- ssm(
    transition = 1, observation = 1, state_cov = 0.01, obs_cov = 0.25,
    init_state = 0, init_cov = 1
  )
  wild <- ssm(
    transition = 1, observation = 1, state_cov = 0.01, obs_cov = 9,
    init_state = 0, init_cov = 1
  )
  withr::local_seed(1L)
  y <- c(stats::rnorm(80, sd = 0.5), stats::rnorm(80, sd = 3))
  fit <- mixture_filter(list(calm, wild), y)

  expect_s3_class(fit, "kalmix::regime_fit")
  expect_equal(dim(fit@regime_prob), c(160L, 2L))
  expect_equal(rowSums(fit@regime_prob), rep(1, 160L), tolerance = 1e-8)

  expect_lt(mean(fit@regime_prob[1:80, 2L]), 0.25)
  expect_gt(mean(fit@regime_prob[81:160, 2L]), 0.6)
})

test_that("mixture_filter() with identical regimes returns the prior", {
  ## When every regime is the same model, the data carry no information about
  ## which regime is active, so the posterior must stay at the (uniform) prior.
  m <- ssm(
    transition = 1, observation = 1, state_cov = 0.1, obs_cov = 1,
    init_state = 0, init_cov = 5
  )
  withr::local_seed(2L)
  y <- stats::rnorm(60)
  fit <- mixture_filter(list(m, m), y)
  expect_equal(fit@regime_prob[, 1L], rep(0.5, 60L), tolerance = 1e-6)
})

test_that("mixture_filter() collapsed state matches a single Kalman filter", {
  ## With identical regimes the collapsed mixture state must equal the plain
  ## single-regime Kalman filtered state -- an internal-consistency oracle.
  m <- ssm(
    transition = 0.9, observation = 1, state_cov = 0.2, obs_cov = 1,
    init_state = 0, init_cov = 5
  )
  withr::local_seed(3L)
  y <- as.numeric(stats::arima.sim(list(ar = 0.5), n = 80))
  mix <- mixture_filter(list(m, m), y)
  one <- kalman_filter(m, y)
  expect_equal(mix@state_mean[, 1L], one@filtered_mean[, 1L], tolerance = 1e-6)
})

test_that("collapsed regime innovations match the Kalman innovations when regimes agree", {
  ## Internal-consistency oracle for the new innovation exposure: with
  ## identical regimes the mixture predictive collapses to the single model's
  ## predictive, so the collapsed innovations and their covariances must
  ## equal the plain Kalman filter's.
  m <- ssm(
    transition = 0.9, observation = 1, state_cov = 0.2, obs_cov = 1,
    init_state = 0, init_cov = 5
  )
  withr::local_seed(14L)
  y <- as.numeric(stats::arima.sim(list(ar = 0.5), n = 80))
  mix <- mixture_filter(list(m, m), y)
  one <- kalman_filter(m, y)
  expect_equal(mix@innovation[, 1L], one@innovation[, 1L], tolerance = 1e-8)
  for (t in c(1L, 40L, 80L)) {
    expect_equal(mix@innovation_cov[[t]], one@innovation_cov[[t]],
                 tolerance = 1e-8)
  }
})

test_that("a regime fit is certified by innovation_diagnostics()", {
  ## Known DGP: one latent AR(1) state observed under Markov-switching
  ## observation noise, filtered with the true regime models. The mixture-PIT
  ## battery must certify the correct specification, and a frozen model that
  ## cannot track the series must fail it. This is the certification path the
  ## regime read-outs needed: standardising by the collapsed Gaussian would
  ## wrongly fail the correct mixture (its innovations are heavy-tailed
  ## relative to one Gaussian by construction).
  withr::local_seed(15L)
  n <- 400L
  stay <- 0.95
  regime <- integer(n)
  regime[1L] <- 1L
  for (t in 2:n) {
    regime[t] <- if (stats::runif(1) < stay) regime[t - 1L] else
      3L - regime[t - 1L]
  }
  x <- numeric(n)
  x[1L] <- stats::rnorm(1, sd = 1)
  for (t in 2:n) {
    x[t] <- 0.7 * x[t - 1L] + stats::rnorm(1, sd = 0.5)
  }
  obs_sd <- c(0.5, 2)[regime]
  y <- x + stats::rnorm(n, sd = obs_sd)

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
  diag_ok <- innovation_diagnostics(fit)
  expect_identical(diag_ok@obs_family, "gaussian_mixture")
  expect_true(diag_ok@adequate)

  ## The regime posterior must also track the true noise regime.
  expect_gt(mean((fit@regime_prob[, 2L] > 0.5) == (regime == 2L)), 0.7)

  frozen <- ssm(
    transition = 1, observation = 1, state_cov = 1e-8, obs_cov = 0.01,
    init_state = 0, init_cov = 0.01
  )
  bad <- mixture_filter(list(frozen, frozen), y)
  diag_bad <- innovation_diagnostics(bad)
  expect_false(diag_bad@adequate)
})

test_that("mixture_filter() log-likelihood matches a single Kalman filter", {
  ## With identical regimes the per-step mixture evidence is the single
  ## model's predictive density (the regime weights sum to one over identical
  ## components), so the log-likelihoods must agree exactly, not just the
  ## collapsed means.
  m <- ssm(
    transition = 0.9, observation = 1, state_cov = 0.2, obs_cov = 1,
    init_state = 0, init_cov = 5
  )
  withr::local_seed(3L)
  y <- as.numeric(stats::arima.sim(list(ar = 0.5), n = 80))
  mix <- mixture_filter(list(m, m), y)
  one <- kalman_filter(m, y)
  expect_equal(mix@log_lik, one@log_lik, tolerance = 1e-8)
})

test_that("mixture_filter() validates models, transition and init_prob", {
  m <- ssm(
    transition = 1, observation = 1, state_cov = 0.1, obs_cov = 1,
    init_state = 0, init_cov = 5
  )
  expect_error(mixture_filter(list(m), 1:10), "at least two")
  expect_error(mixture_filter(list(m, "x"), 1:10), "ssm")
  expect_error(
    mixture_filter(list(m, m), 1:10, transition = matrix(0.5, 3, 3)),
    "2 by 2"
  )
  expect_error(
    mixture_filter(list(m, m), 1:10, init_prob = c(0.3, 0.3)),
    "sum to one"
  )
})

