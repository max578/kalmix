test_that("intercept models equal their manually de-meaned counterparts", {
  ## Independent oracle by exact algebra: with state intercept c and
  ## observation intercept d, writing x_t = z_t + m_t against the
  ## deterministic mean path m_t = A m_{t-1} + c (m_0 = the prior mean)
  ## reduces the model to the intercept-free one on y'_t = y_t - B m_t - d,
  ## with the filtered means shifted by m_t and every covariance unchanged.
  a <- 0.85
  b <- 1.3
  q <- 0.3
  r <- 0.7
  mu0 <- 2
  p0 <- 4
  ci <- 0.6
  di <- -1.1
  n <- 60L

  withr::local_seed(31L)
  y <- cumsum(stats::rnorm(n))

  with_int <- ssm(
    transition = a, observation = b, state_cov = q, obs_cov = r,
    init_state = mu0, init_cov = p0,
    state_intercept = ci, obs_intercept = di
  )
  fit1 <- kalman_filter(with_int, y)

  m_path <- numeric(n)
  m_prev <- mu0
  for (t in seq_len(n)) {
    m_path[t] <- a * m_prev + ci
    m_prev <- m_path[t]
  }
  plain <- ssm(
    transition = a, observation = b, state_cov = q, obs_cov = r,
    init_state = 0, init_cov = p0
  )
  fit2 <- kalman_filter(plain, y - b * m_path - di)

  expect_equal(fit1@filtered_mean[, 1L], fit2@filtered_mean[, 1L] + m_path,
               tolerance = 1e-10)
  expect_equal(fit1@innovation, fit2@innovation, tolerance = 1e-10)
  expect_equal(fit1@log_lik, fit2@log_lik, tolerance = 1e-10)
  for (t in c(1L, 30L, 60L)) {
    expect_equal(fit1@filtered_cov[[t]], fit2@filtered_cov[[t]],
                 tolerance = 1e-12)
  }

  ## The smoother inherits the same shift.
  sm1 <- rts_smoother(fit1)
  sm2 <- rts_smoother(fit2)
  expect_equal(sm1@smoothed_mean[, 1L], sm2@smoothed_mean[, 1L] + m_path,
               tolerance = 1e-10)

  ## A time-varying observation intercept equals pre-subtracting it.
  d_path <- sin(seq_len(n) / 5)
  tv <- ssm(
    transition = a, observation = b, state_cov = q, obs_cov = r,
    init_state = 0, init_cov = p0,
    obs_intercept = as.list(d_path)
  )
  fit_tv <- kalman_filter(tv, y)
  fit_ref <- kalman_filter(plain, y - d_path)
  expect_equal(fit_tv@filtered_mean, fit_ref@filtered_mean, tolerance = 1e-10)
})

test_that("the intercept-form regime filter equals the deviation-form spread filter", {
  ## Cross-check of the two shipped mechanisms: kalmix_filter() de-means the
  ## spread and shifts the state back, while the native intercept form gives
  ## each regime the mean-reversion pull c = theta (1 - phi) directly. The
  ## two are the same linear system in different coordinates, so the filtered
  ## states and regime probabilities must agree.
  x <- 5 + .sim_ou(n = 1500, seed = 32L)
  model <- ou_fit(x)
  fit_dev <- kalmix_filter(model, x)

  phi <- exp(-model@kappa * model@dt)
  innov_var <- model@sigma^2 * (1 - phi^2) / (2 * model@kappa)
  base_obs <- max(innov_var, .Machine$double.eps)
  regimes <- lapply(c(0.5, 1), function(mult) {
    ssm(
      transition = phi, observation = 1,
      state_cov = innov_var, obs_cov = base_obs * mult,
      state_intercept = model@theta * (1 - phi),
      init_state = model@theta,
      init_cov = model@sigma^2 / (2 * model@kappa)
    )
  })
  fit_int <- mixture_filter(regimes, x)

  expect_equal(fit_int@state_mean[, 1L], fit_dev@state_mean[, 1L],
               tolerance = 1e-8)
  expect_equal(fit_int@regime_prob, fit_dev@regime_prob, tolerance = 1e-8)
})

test_that("ssm_fit() carries intercepts through the rebuild", {
  withr::local_seed(33L)
  n <- 400L
  theta <- 4
  phi <- 0.8
  x <- numeric(n)
  x[1L] <- theta
  for (t in 2:n) {
    x[t] <- theta * (1 - phi) + phi * x[t - 1L] + stats::rnorm(1, sd = 0.5)
  }
  y <- x + stats::rnorm(n, sd = 0.6)

  template <- ssm(
    transition = phi, observation = 1, state_cov = 1, obs_cov = 1,
    state_intercept = theta * (1 - phi),
    init_state = theta, init_cov = 5
  )
  fit <- ssm_fit(template, y)
  expect_identical(fit@convergence, 0L)
  expect_equal(fit@model@state_intercept[[1L]], theta * (1 - phi))
  expect_lt(abs(fit@estimates$state_cov[1L, 1L] - 0.25), 0.15)
  expect_lt(abs(fit@estimates$obs_cov[1L, 1L] - 0.36), 0.2)
})

test_that("intercept arguments are validated", {
  expect_error(
    ssm(
      transition = 1, observation = 1, state_cov = 1, obs_cov = 1,
      init_state = c(0, 0), init_cov = diag(2),
      state_intercept = c(1, 2, 3)
    ),
    "length"
  )
  expect_error(
    ssm(
      transition = 1, observation = 1, state_cov = 1, obs_cov = 1,
      init_state = 0, init_cov = 1,
      obs_intercept = "x"
    ),
    "numeric"
  )
})
