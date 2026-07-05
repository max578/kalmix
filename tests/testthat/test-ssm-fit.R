test_that("ssm_fit() recovers known noise parameters in simulation", {
  ## Known DGP: local level with q = 0.25, r = 1. The fitted variances must
  ## land near the truth on a long series, and the delta-method intervals
  ## must cover it.
  withr::local_seed(21L)
  n <- 800L
  q_true <- 0.25
  r_true <- 1
  level <- cumsum(stats::rnorm(n, sd = sqrt(q_true)))
  y <- level + stats::rnorm(n, sd = sqrt(r_true))

  template <- ssm(
    transition = 1, observation = 1, state_cov = 1, obs_cov = 1,
    init_state = y[1L], init_cov = 10
  )
  fit <- ssm_fit(template, y)

  expect_s3_class(fit, "kalmix::ssm_mle")
  expect_identical(fit@convergence, 0L)
  q_hat <- fit@estimates$state_cov[1L, 1L]
  r_hat <- fit@estimates$obs_cov[1L, 1L]
  q_se <- fit@std_errors$state_cov[1L, 1L]
  r_se <- fit@std_errors$obs_cov[1L, 1L]
  expect_gt(q_hat, 0)
  expect_gt(r_hat, 0)
  expect_lt(abs(q_hat - q_true), 4 * q_se)
  expect_lt(abs(r_hat - r_true), 4 * r_se)

  ## The fitted model is a working ssm: the filter runs and its innovations
  ## pass the adequacy battery on the correctly specified series.
  expect_true(innovation_diagnostics(kalman_filter(fit@model, y))@adequate)
})

test_that("ssm_fit() agrees with StructTS on the local level", {
  ## Independent competitor oracle: stats::StructTS fits the same local-level
  ## model by maximum likelihood with its own diffuse-style initialisation,
  ## so on a long series (where the initialisation washes out) the fitted
  ## variances must agree to within a few per cent.
  withr::local_seed(22L)
  n <- 600L
  level <- cumsum(stats::rnorm(n, sd = 0.6))
  y <- level + stats::rnorm(n, sd = 1.2)

  ref <- stats::StructTS(y, type = "level")
  template <- ssm(
    transition = 1, observation = 1, state_cov = 1, obs_cov = 1,
    init_state = y[1L], init_cov = 1e7
  )
  fit <- ssm_fit(template, y)

  expect_equal(
    fit@estimates$state_cov[1L, 1L],
    unname(ref$coef["level"]),
    tolerance = 0.05
  )
  expect_equal(
    fit@estimates$obs_cov[1L, 1L],
    unname(ref$coef["epsilon"]),
    tolerance = 0.05
  )
})

test_that("ssm_fit() estimates are equivariant under scaling", {
  ## Metamorphic property: scaling the series by c scales both fitted
  ## variances by c^2 and leaves nothing else to chance.
  withr::local_seed(23L)
  n <- 400L
  y <- cumsum(stats::rnorm(n, sd = 0.5)) + stats::rnorm(n)
  template <- ssm(
    transition = 1, observation = 1, state_cov = 1, obs_cov = 1,
    init_state = y[1L], init_cov = 10
  )
  fit1 <- ssm_fit(template, y)

  cc <- 3
  template2 <- ssm(
    transition = 1, observation = 1, state_cov = 1, obs_cov = 1,
    init_state = cc * y[1L], init_cov = cc^2 * 10
  )
  fit2 <- ssm_fit(template2, cc * y)

  expect_equal(
    fit2@estimates$state_cov[1L, 1L],
    cc^2 * fit1@estimates$state_cov[1L, 1L],
    tolerance = 1e-3
  )
  expect_equal(
    fit2@estimates$obs_cov[1L, 1L],
    cc^2 * fit1@estimates$obs_cov[1L, 1L],
    tolerance = 1e-3
  )
})

test_that("ssm_fit() recovers a stable transition when asked", {
  withr::local_seed(24L)
  n <- 800L
  a_true <- 0.8
  x <- as.numeric(stats::arima.sim(list(ar = a_true), n = n, sd = sqrt(0.4)))
  y <- x + stats::rnorm(n, sd = 0.5)
  template <- ssm(
    transition = 0.5, observation = 1, state_cov = 1, obs_cov = 1,
    init_state = 0, init_cov = 5
  )
  fit <- ssm_fit(template, y, transition = TRUE)
  a_hat <- fit@estimates$transition[1L, 1L]
  expect_lt(abs(a_hat - a_true), 0.1)
})

test_that("ssm_fit() fits a gappy series through the missing-data filter", {
  withr::local_seed(25L)
  n <- 500L
  level <- cumsum(stats::rnorm(n, sd = 0.5))
  y <- level + stats::rnorm(n)
  y[seq(25L, 475L, by = 25L)] <- NA_real_
  template <- ssm(
    transition = 1, observation = 1, state_cov = 1, obs_cov = 1,
    init_state = 0, init_cov = 10
  )
  fit <- ssm_fit(template, y)
  expect_identical(fit@convergence, 0L)
  expect_gt(fit@estimates$state_cov[1L, 1L], 0.05)
  expect_lt(abs(fit@estimates$obs_cov[1L, 1L] - 1), 0.5)
})

test_that("ssm_fit() validates its inputs", {
  y <- cumsum(stats::rnorm(50))
  template <- ssm(
    transition = 1, observation = 1, state_cov = 1, obs_cov = 1,
    init_state = 0, init_cov = 10
  )
  expect_error(ssm_fit("nope", y), "ssm")
  expect_error(ssm_fit(template, y, estimate = "banana"), "arg")
  expect_error(ssm_fit(template, y, n_start = 0L), "positive integer")
  tv <- ssm(
    transition = rep(list(matrix(1)), 50L), observation = 1,
    state_cov = 1, obs_cov = 1, init_state = 0, init_cov = 10
  )
  expect_error(
    ssm_fit(tv, y, transition = TRUE),
    "time-invariant"
  )
})
