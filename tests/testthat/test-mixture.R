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

