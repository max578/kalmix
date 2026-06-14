## The online (fixed-point smoother) engine and the expanding (re-smoothing)
## engine compute the same objective causal information rate of Andreou, Chen
## and Bollt (2026, eq. 9) by two independent constructions, so their agreement
## is a self-checking oracle for the recursion. The online engine additionally
## handles time-varying models, which the truncating expanding construction
## cannot.

test_that("online and expanding CIR engines agree to numerical precision", {
  model <- ssm(
    transition = matrix(c(0.6, -0.5, 0.5, 0.6), nrow = 2L),
    observation = matrix(c(1, 0), nrow = 1L),
    state_cov = diag(c(0.2, 0.2)), obs_cov = 1,
    init_state = c(0, 0), init_cov = diag(c(5, 5))
  )
  withr::local_seed(1L)
  x <- numeric(300L)
  s <- c(0, 0)
  for (t in seq_len(300L)) {
    s <- as.numeric(model@transition[[1L]] %*% s) + stats::rnorm(2L, sd = 0.45)
    x[t] <- s[1L] + stats::rnorm(1L)
  }
  on <- causal_information_rate(model, x, engine = "online")
  ex <- causal_information_rate(model, x, engine = "expanding")
  expect_equal(on, ex, tolerance = 1e-8)
  expect_gt(on, 0)
})

test_that("the engines agree on a univariate local level too", {
  model <- ssm(
    transition = 1, observation = 1, state_cov = 0.04, obs_cov = 1,
    init_state = 0, init_cov = 10
  )
  withr::local_seed(3L)
  y <- cumsum(stats::rnorm(500L, sd = 0.2)) + stats::rnorm(500L)
  expect_equal(
    causal_information_rate(model, y, engine = "online"),
    causal_information_rate(model, y, engine = "expanding"),
    tolerance = 1e-8
  )
})

test_that("the online engine handles a time-varying model", {
  n <- 120L
  a_list <- lapply(
    seq_len(n),
    function(t) matrix(0.9 + 0.05 * sin(t / 10), 1L, 1L)
  )
  model <- ssm(
    transition = a_list, observation = 1, state_cov = 0.1, obs_cov = 1,
    init_state = 0, init_cov = 10
  )
  withr::local_seed(4L)
  y <- as.numeric(stats::arima.sim(list(ar = 0.9), n = n)) + stats::rnorm(n)
  cir <- causal_information_rate(model, y)
  expect_true(is.finite(cir) && cir >= 0)
  expect_error(
    causal_information_rate(model, y, engine = "expanding"), "time-invariant"
  )
})

test_that("aci() uses the online engine, matching the expanding cross-check", {
  model <- ssm(
    transition = matrix(c(0.6, -0.5, 0.5, 0.6), nrow = 2L),
    observation = matrix(c(1, 0), nrow = 1L),
    state_cov = diag(c(0.2, 0.2)), obs_cov = 1,
    init_state = c(0, 0), init_cov = diag(c(5, 5))
  )
  withr::local_seed(2L)
  x <- numeric(250L)
  s <- c(0, 0)
  for (t in seq_len(250L)) {
    s <- as.numeric(model@transition[[1L]] %*% s) + stats::rnorm(2L, sd = 0.45)
    x[t] <- s[1L] + stats::rnorm(1L)
  }
  fit <- aci(model, x)
  expect_equal(
    fit@lead_time,
    causal_information_rate(model, x, engine = "expanding"),
    tolerance = 1e-8
  )
})
