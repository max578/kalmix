# Tests for the Student-t observation model: the robust filter, the profile-
# likelihood degrees-of-freedom estimator, and the family-aware adequacy guard.

local_level <- function(obs_family = "gaussian", obs_df = Inf) {
  ssm(
    transition = 1, observation = 1,
    state_cov = 0.04, obs_cov = 1,
    init_state = 0, init_cov = 10,
    obs_family = obs_family, obs_df = obs_df
  )
}

# --- Class surface ------------------------------------------------------------

test_that("ssm defaults to a gaussian observation family", {
  m <- local_level()
  expect_identical(m@obs_family, "gaussian")
  expect_identical(m@obs_df, Inf)
})

test_that("ssm validates the observation family and degrees of freedom", {
  expect_error(
    ssm(transition = 1, observation = 1, state_cov = 0.04, obs_cov = 1,
        init_state = 0, init_cov = 10, obs_family = "student_t", obs_df = 2),
    "greater than 2"
  )
  expect_error(
    ssm(transition = 1, observation = 1, state_cov = 0.04, obs_cov = 1,
        init_state = 0, init_cov = 10, obs_family = "student_t", obs_df = Inf),
    "greater than 2"
  )
  # A gaussian family ignores a supplied df and holds it at Inf.
  m <- ssm(transition = 1, observation = 1, state_cov = 0.04, obs_cov = 1,
           init_state = 0, init_cov = 10, obs_family = "gaussian", obs_df = 5)
  expect_identical(m@obs_df, Inf)
})

# --- Continuity with the Gaussian filter --------------------------------------

test_that("the Student-t filter approaches the Gaussian filter as df grows", {
  set.seed(1)
  y <- cumsum(rnorm(300, sd = 0.2)) + rnorm(300)
  g <- kalman_filter(local_level(), y)
  t_big <- kalman_filter(local_level("student_t", 1e6), y)
  expect_equal(t_big@log_lik, g@log_lik, tolerance = 1e-3)
  expect_equal(t_big@filtered_mean, g@filtered_mean, tolerance = 1e-3)
})

test_that("the gaussian path is unchanged by the new branch", {
  set.seed(7)
  y <- cumsum(rnorm(200, sd = 0.3)) + rnorm(200)
  fit <- kalman_filter(local_level(), y)
  # Reference value of the prediction-error-decomposition log-likelihood.
  expect_equal(fit@log_lik, sum(vapply(seq_len(200), function(t) {
    -0.5 * (log(2 * pi) + log(fit@innovation_cov[[t]][1, 1]) +
              fit@innovation[t, 1]^2 / fit@innovation_cov[[t]][1, 1])
  }, numeric(1))))
})

# --- Independent oracle for the t density -------------------------------------

test_that("the multivariate-t log-density matches an independent implementation", {
  skip_if_not_installed("mvtnorm")
  set.seed(3)
  scale <- matrix(c(2, 0.5, 0.5, 1), 2)
  e <- c(1.3, -0.7)
  for (nu in c(3, 7, 20)) {
    expect_equal(
      kalmix:::.mvt_logpdf(e, scale, nu),
      mvtnorm::dmvt(e, delta = c(0, 0), sigma = scale, df = nu, log = TRUE),
      tolerance = 1e-6
    )
  }
})

# --- Robustness to outliers ---------------------------------------------------

test_that("the Student-t filter is more robust to observation outliers", {
  set.seed(11)
  n <- 300
  level <- cumsum(rnorm(n, sd = 0.1))
  y <- level + rnorm(n, sd = 0.5)
  # Inject heavy outliers the Gaussian filter cannot discount.
  hits <- c(60, 120, 180, 240)
  y[hits] <- y[hits] + c(12, -10, 14, -11)
  m_obs <- 0.25
  g <- kalman_filter(ssm(1, 1, 0.01, m_obs, 0, 10), y)
  t <- kalman_filter(ssm(1, 1, 0.01, m_obs, 0, 10,
                         obs_family = "student_t", obs_df = 4), y)
  rmse <- function(fit) sqrt(mean((fit@filtered_mean[, 1] - level)^2))
  expect_lt(rmse(t), rmse(g))
})

# --- Degrees-of-freedom estimation --------------------------------------------

test_that("estimate_obs_df selects a finite df on heavy-tailed data", {
  set.seed(21)
  y <- cumsum(rnorm(400, sd = 0.2)) + rt(400, df = 3)
  est <- estimate_obs_df(local_level(), y)
  expect_identical(est$family, "student_t")
  expect_true(is.finite(est$df) && est$df <= 10)
  expect_length(est$profile, length(est$grid))
})

test_that("estimate_obs_df selects the Gaussian model on light-tailed data", {
  set.seed(22)
  y <- cumsum(rnorm(400, sd = 0.2)) + rnorm(400)
  est <- estimate_obs_df(local_level(), y)
  expect_identical(est$family, "gaussian")
  expect_identical(est$df, Inf)
})

# --- The headline: the guard certifies under the right family -----------------

test_that("the adequacy guard abstains under gaussian but certifies under t", {
  set.seed(31)
  level <- cumsum(rnorm(500, sd = 0.15))
  y <- level + rt(500, df = 4)
  ad_g <- innovation_diagnostics(kalman_filter(local_level(), y))
  expect_false(isTRUE(ad_g@adequate))
  expect_identical(ad_g@obs_family, "gaussian")

  est <- estimate_obs_df(local_level(), y)
  ad_t <- innovation_diagnostics(
    kalman_filter(local_level("student_t", est$df), y)
  )
  expect_true(isTRUE(ad_t@adequate))
  expect_identical(ad_t@obs_family, "student_t")
})

test_that("aci grounds a heavy-tailed series only under the t family", {
  set.seed(41)
  level <- cumsum(rnorm(500, sd = 0.15))
  y <- level + rt(500, df = 4)
  fit_g <- aci(local_level(), y, lead_time = FALSE)
  expect_identical(fit_g@grounding, "[unverified]")

  est <- estimate_obs_df(local_level(), y)
  fit_t <- aci(local_level("student_t", est$df), y, lead_time = FALSE,
               mechanism = list(verified_on = as.Date("2026-06-14")))
  expect_true(isTRUE(fit_t@adequacy@adequate))
  expect_identical(fit_t@grounding, "grounded")
})
