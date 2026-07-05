test_that("gappy filter and smoother match direct joint-Gaussian conditioning", {
  ## Independent oracle: for a linear-Gaussian model the states and observed
  ## entries are jointly Gaussian, so the filtered and smoothed moments at
  ## every step follow from one multivariate-normal conditioning -- a
  ## derivation path with no recursion in it. Build the joint moments directly
  ## from the model algebra and condition.
  a <- 0.9
  b <- 1.2
  q <- 0.3
  r <- 0.8
  mu0 <- 1.5
  p0 <- 2
  n <- 8L
  miss <- c(3L, 6L)

  withr::local_seed(9L)
  y <- numeric(n)
  x_true <- mu0 + stats::rnorm(1L) * sqrt(p0)
  for (t in seq_len(n)) {
    x_true <- a * x_true + stats::rnorm(1L, sd = sqrt(q))
    y[t] <- b * x_true + stats::rnorm(1L, sd = sqrt(r))
  }
  y[miss] <- NA_real_
  obs <- setdiff(seq_len(n), miss)

  ## Joint moments of (x_1, ..., x_n): mean a^t mu0; Cov(x_s, x_t) for s <= t
  ## is a^(t-s) V_s with V_t = a^2 V_{t-1} + q from V_0 = p0.
  v <- numeric(n)
  v_prev <- p0
  for (t in seq_len(n)) {
    v[t] <- a^2 * v_prev + q
    v_prev <- v[t]
  }
  cov_xx <- matrix(0, n, n)
  for (s in seq_len(n)) {
    for (t in seq_len(n)) {
      cov_xx[s, t] <- a^(abs(t - s)) * v[min(s, t)]
    }
  }
  mean_x <- a^(seq_len(n)) * mu0
  cov_xy <- b * cov_xx[, obs, drop = FALSE]
  cov_yy <- b^2 * cov_xx[obs, obs, drop = FALSE] + diag(r, length(obs))
  mean_y <- b * mean_x[obs]

  model <- ssm(
    transition = a, observation = b, state_cov = q, obs_cov = r,
    init_state = mu0, init_cov = p0
  )
  fit <- kalman_filter(model, y)
  sm <- rts_smoother(fit)

  condition_on <- function(t, idx) {
    keep <- match(idx, obs)
    s22 <- cov_yy[keep, keep, drop = FALSE]
    s12 <- cov_xy[t, keep, drop = FALSE]
    gain <- s12 %*% solve(s22)
    list(
      mean = mean_x[t] + as.numeric(gain %*% (y[idx] - mean_y[keep])),
      var = cov_xx[t, t] - as.numeric(gain %*% t(s12))
    )
  }

  for (t in seq_len(n)) {
    past <- obs[obs <= t]
    oracle_f <- condition_on(t, past)
    expect_equal(fit@filtered_mean[t, 1L], oracle_f$mean, tolerance = 1e-6)
    expect_equal(fit@filtered_cov[[t]][1L, 1L], oracle_f$var, tolerance = 1e-6)
    oracle_s <- condition_on(t, obs)
    expect_equal(sm@smoothed_mean[t, 1L], oracle_s$mean, tolerance = 1e-6)
    expect_equal(sm@smoothed_cov[[t]][1L, 1L], oracle_s$var, tolerance = 1e-6)
  }

  ## The log-likelihood must equal the joint Gaussian density of the observed
  ## entries -- the prediction-error decomposition against the direct form.
  skip_if_not_installed("mvtnorm")
  direct <- mvtnorm::dmvnorm(y[obs], mean = mean_y, sigma = cov_yy, log = TRUE)
  expect_equal(fit@log_lik, direct, tolerance = 1e-8)
})

test_that("a missing step is exactly a prediction-only step", {
  model <- ssm(
    transition = 1, observation = 1, state_cov = 0.1, obs_cov = 1,
    init_state = 0, init_cov = 5
  )
  withr::local_seed(10L)
  y <- cumsum(stats::rnorm(40, sd = 0.4)) + stats::rnorm(40)
  y[15L] <- NA_real_
  fit <- kalman_filter(model, y)
  expect_equal(fit@filtered_mean[15L, ], fit@predicted_mean[15L, ])
  expect_equal(fit@filtered_cov[[15L]], fit@predicted_cov[[15L]])
  expect_true(is.na(fit@innovation[15L, ]))
  expect_false(anyNA(fit@filtered_mean))

  ## Student-t family: the robust filter skips the same way.
  tmodel <- ssm(
    transition = 1, observation = 1, state_cov = 0.1, obs_cov = 1,
    init_state = 0, init_cov = 5, obs_family = "student_t", obs_df = 5
  )
  tfit <- kalman_filter(tmodel, y)
  expect_equal(tfit@filtered_mean[15L, ], tfit@predicted_mean[15L, ])
  expect_true(is.na(tfit@innovation[15L, ]))
})

test_that("partially missing multivariate rows are refused", {
  model <- ssm(
    transition = diag(2), observation = diag(2),
    state_cov = diag(0.1, 2), obs_cov = diag(1, 2),
    init_state = c(0, 0), init_cov = diag(5, 2)
  )
  y <- matrix(stats::rnorm(20), ncol = 2L)
  y[4L, 2L] <- NA_real_
  expect_error(kalman_filter(model, y), "fully observed or fully missing")
})

test_that("mixture_filter() propagates the chain prior through a gap", {
  calm <- ssm(
    transition = 1, observation = 1, state_cov = 0.01, obs_cov = 0.25,
    init_state = 0, init_cov = 1
  )
  wild <- ssm(
    transition = 1, observation = 1, state_cov = 0.01, obs_cov = 4,
    init_state = 0, init_cov = 1
  )
  withr::local_seed(11L)
  y <- c(stats::rnorm(30, sd = 0.5), stats::rnorm(30, sd = 2))
  y[20L] <- NA_real_
  trans <- matrix(c(0.95, 0.05, 0.05, 0.95), 2L)
  fit <- mixture_filter(list(calm, wild), y, transition = trans)

  ## At the gap the data carry no evidence, so the regime posterior is the
  ## previous posterior pushed through the chain alone.
  expect_equal(
    fit@regime_prob[20L, ],
    as.numeric(fit@regime_prob[19L, ] %*% trans),
    tolerance = 1e-10
  )
  expect_true(anyNA(fit@innovation[20L, ]))
  expect_false(anyNA(fit@state_mean))
  expect_equal(rowSums(fit@regime_prob), rep(1, 60L), tolerance = 1e-8)
})

test_that("hmm_filter() and hmm_viterbi() treat NA as a flat emission", {
  sim <- .sim_hmm(n = 150, means = c(0, 4), sd = 1, seed = 12L)
  y <- sim$y
  y[c(40L, 41L, 90L)] <- NA_real_
  model <- hmm(
    init_prob = c(0.5, 0.5), transition = sim$transition,
    emission_mean = sim$means, emission_sd = c(sim$sd, sim$sd)
  )
  fit <- hmm_filter(model, y)
  expect_false(anyNA(fit@filtered))
  expect_true(is.finite(fit@log_lik))

  ## At a gap the filtered state distribution is the previous one pushed
  ## through the transition alone.
  expect_equal(
    fit@filtered[40L, ],
    as.numeric(fit@filtered[39L, ] %*% sim$transition),
    tolerance = 1e-10
  )

  path <- hmm_viterbi(model, y)
  expect_equal(length(path), length(y))
  expect_gt(mean(path == sim$state), 0.8)

  expect_error(hmm_filter(model, rep(NA_real_, 10)), "at least one observed")
})

test_that("innovation_diagnostics() skips missing steps and reports the rest", {
  model <- ssm(
    transition = 1, observation = 1, state_cov = 0.04, obs_cov = 1,
    init_state = 0, init_cov = 10
  )
  withr::local_seed(13L)
  level <- cumsum(stats::rnorm(300, sd = 0.2))
  y <- level + stats::rnorm(300)
  y[c(50L, 150L, 250L)] <- NA_real_
  diag_fit <- innovation_diagnostics(kalman_filter(model, y))
  expect_true(diag_fit@adequate)
  expect_identical(diag_fit@n_obs, 297L)
  expect_identical(sum(is.na(diag_fit@standardised)), 3L)
})
