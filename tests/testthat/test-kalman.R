test_that("kalman_filter() returns the right shapes", {
  m <- ssm(
    transition = 1, observation = 1, state_cov = 0.01, obs_cov = 1,
    init_state = 0, init_cov = 10
  )
  sim <- .sim_local_level(n = 120, seed = 1L)
  fit <- kalman_filter(m, sim$y)
  expect_s3_class(fit, "kalmix::kalman_fit")
  expect_equal(dim(fit@filtered_mean), c(120L, 1L))
  expect_length(fit@filtered_cov, 120L)
  expect_true(is.finite(fit@log_lik))
  expect_output(print(fit), "kalman_fit")
})

test_that("kalman_filter() recovers a known latent level", {
  ## Known DGP: the synthetic level is the truth the filter must track. The
  ## oracle is the simulated state, not the code.
  sim <- .sim_local_level(n = 300, state_sd = 0.3, obs_sd = 1, seed = 2L)
  m <- ssm(
    transition = 1, observation = 1,
    state_cov = sim$state_sd^2, obs_cov = sim$obs_sd^2,
    init_state = sim$y[1L], init_cov = 10
  )
  fit <- kalman_filter(m, sim$y)
  rmse_filter <- sqrt(mean((fit@filtered_mean[, 1L] - sim$level)^2))
  rmse_naive <- sqrt(mean((sim$y - sim$level)^2))
  ## The filter must beat the raw observation as an estimate of the level.
  expect_lt(rmse_filter, rmse_naive)
  expect_gt(cor(fit@filtered_mean[, 1L], sim$level), 0.8)
})

test_that("kalman_filter() recovers a known constant velocity", {
  ## Known DGP: position grows at a constant velocity of 0.5 per step; the
  ## two-state model must recover that slope.
  withr::local_seed(3L)
  velocity <- 0.5
  truth <- velocity * seq_len(150)
  y <- truth + stats::rnorm(150, sd = 1)
  m <- ssm(
    transition = matrix(c(1, 0, 1, 1), nrow = 2L),
    observation = matrix(c(1, 0), nrow = 1L),
    state_cov = diag(c(1e-5, 1e-5)),
    obs_cov = 1,
    init_state = c(0, 1),
    init_cov = diag(c(10, 10))
  )
  fit <- kalman_filter(m, y)
  expect_equal(fit@filtered_mean[150L, 2L], velocity, tolerance = 0.1)
})

test_that("kalman_filter() log-likelihood matches a direct Gaussian density", {
  ## For a static-mean model the innovations are independent Gaussians, so the
  ## accumulated prediction-error log-likelihood must equal the explicit sum of
  ## dnorm() log-densities -- an oracle independent of the recursion.
  withr::local_seed(4L)
  y <- stats::rnorm(50, mean = 0, sd = 2)
  m <- ssm(
    transition = 0, observation = 1, state_cov = 0, obs_cov = 4,
    init_state = 0, init_cov = 0
  )
  fit <- kalman_filter(m, y)
  direct <- sum(stats::dnorm(y, mean = 0, sd = 2, log = TRUE))
  expect_equal(fit@log_lik, direct, tolerance = 1e-8)
})

test_that("rts_smoother() is no worse than the filter and accepts both inputs", {
  sim <- .sim_local_level(n = 250, state_sd = 0.3, obs_sd = 1, seed = 5L)
  m <- ssm(
    transition = 1, observation = 1,
    state_cov = sim$state_sd^2, obs_cov = sim$obs_sd^2,
    init_state = sim$y[1L], init_cov = 10
  )
  fit <- kalman_filter(m, sim$y)
  sm <- rts_smoother(fit)
  expect_s3_class(sm, "kalmix::rts_fit")

  rmse_filter <- sqrt(mean((fit@filtered_mean[, 1L] - sim$level)^2))
  rmse_smooth <- sqrt(mean((sm@smoothed_mean[, 1L] - sim$level)^2))
  ## Smoothing uses the whole series, so it cannot be worse than filtering.
  expect_lte(rmse_smooth, rmse_filter + 1e-8)

  ## Passing the model plus observations must give the same answer.
  sm2 <- rts_smoother(m, sim$y)
  expect_equal(sm@smoothed_mean, sm2@smoothed_mean)
})

test_that("rts_smoother() covariances never exceed the filter's (trace-wise)", {
  ## Conditioning on the whole series cannot add uncertainty: at every step
  ## the smoothed covariance is dominated by the filtered covariance, so its
  ## trace must be no larger. Checked as a property over random stable one-
  ## and two-dimensional models rather than a single fixture (the covariance
  ## recursions do not depend on the observed values, only on the model).
  for (seed in 1:5) {
    withr::local_seed(seed)
    m_dim <- 1L + seed %% 2L
    a <- matrix(stats::rnorm(m_dim^2), m_dim)
    a <- a / (max(Mod(eigen(a, only.values = TRUE)$values)) + 0.1)
    q_root <- matrix(stats::rnorm(m_dim^2, sd = 0.5), m_dim)
    model <- ssm(
      transition = a,
      observation = matrix(stats::rnorm(m_dim), nrow = 1L),
      state_cov = crossprod(q_root) + diag(0.05, m_dim),
      obs_cov = 0.5 + stats::runif(1),
      init_state = numeric(m_dim),
      init_cov = diag(1, m_dim)
    )
    fit <- kalman_filter(model, stats::rnorm(60))
    sm <- rts_smoother(fit)
    trace_filter <- vapply(
      fit@filtered_cov, function(p) sum(diag(p)), numeric(1L)
    )
    trace_smooth <- vapply(
      sm@smoothed_cov, function(p) sum(diag(p)), numeric(1L)
    )
    expect_true(all(trace_smooth <= trace_filter + 1e-10))
  }
})

test_that("kalman_filter() and rts_smoother() validate their inputs", {
  m <- ssm(
    transition = 1, observation = 1, state_cov = 0.1, obs_cov = 1,
    init_state = 0, init_cov = 10
  )
  expect_error(kalman_filter("not a model", 1:10), "ssm")
  expect_error(kalman_filter(m, matrix(1, 5, 2)), "column")
  expect_error(kalman_filter(m, c(1, NA, 3)), "missing")
  expect_error(rts_smoother(m), "required")
  expect_error(rts_smoother("nope"), "kalman_fit")
})
