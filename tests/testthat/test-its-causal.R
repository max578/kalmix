test_that("its_causal() recovers a known post-intervention lift", {
  ## Known DGP: a smooth pre-period, then a post-period with a genuine +3 lift
  ## added. The estimated average effect must recover that lift, and the
  ## interval must exclude zero.
  withr::local_seed(1L)
  pre <- cumsum(stats::rnorm(60, sd = 0.2))
  post <- pre[60L] + cumsum(stats::rnorm(40, sd = 0.2)) + 3
  fit <- its_causal(c(pre, post), intervention = 61L)
  expect_s3_class(fit, "kalmix::its_fit")
  expect_equal(fit@avg_effect, 3, tolerance = 1)
  expect_gt(fit@avg_effect_lower, 0)
  expect_output(print(fit), "average effect")
})

test_that("its_causal() reports a null effect honestly", {
  ## Known DGP: a pure random walk with NO intervention effect. A correctly
  ## calibrated interval must straddle zero rather than assert an effect.
  withr::local_seed(2L)
  x <- cumsum(stats::rnorm(120, sd = 0.3))
  fit <- its_causal(x, intervention = 61L)
  expect_lte(fit@avg_effect_lower, 0)
  expect_gte(fit@avg_effect_upper, 0)
})

test_that("its_causal() interval coverage is near nominal under the null", {
  ## Repeated no-effect random walks: the 95% interval should straddle zero in
  ## roughly 95% of replicates. This is the calibration oracle for the verb --
  ## the maximum-likelihood noise split is what makes it hold. With 300
  ## replicates the binomial standard error around 0.95 is about 0.013, so a
  ## two-sided acceptance band of [0.92, 0.985] sits at roughly 2.5 standard
  ## errors: an anticonservative interval drags the coverage below the floor,
  ## an inflated one pushes it above the ceiling.
  withr::local_seed(5L)
  covered <- vapply(seq_len(300L), function(i) {
    x <- cumsum(stats::rnorm(120, sd = 0.3))
    f <- its_causal(x, intervention = 61L)
    f@avg_effect_lower <= 0 && f@avg_effect_upper >= 0
  }, logical(1L))
  expect_gte(mean(covered), 0.92)
  expect_lte(mean(covered), 0.985)
})

test_that("its_causal() average-effect variance matches a simulation oracle", {
  ## Independent Monte-Carlo oracle for the running-sum forecast variance. The
  ## verb propagates an augmented state (x_t, c_t) whose process noise enters
  ## both blocks as the same draw w_t, so the augmented noise covariance is
  ## [[Q, Q], [Q, Q]] -- a zero off-diagonal block drops the state/running-sum
  ## noise correlation and understates the interval. The oracle simulates the
  ## post-period state paths directly from the fitted pre-period model (no
  ## augmented recursion involved) and compares the variance of the summed
  ## counterfactual with the value implied by the reported interval.
  for (trend in c(FALSE, TRUE)) {
    withr::local_seed(7L)
    slope <- if (trend) 0.3 * seq_len(120L) else 0
    x <- slope + cumsum(stats::rnorm(120, sd = 0.3))
    fit <- its_causal(x, intervention = 81L, trend = trend)
    n_post <- 40L
    z <- stats::qnorm(0.975)
    v_fit <- ((fit@avg_effect_upper - fit@avg_effect) / z * n_post)^2

    ## Rebuild the fitted pre-period model and its final filtered covariance,
    ## exactly as the verb seeds its forecast.
    pre <- x[seq_len(80L)]
    model <- kalmix:::.its_pre_model(pre, trend)
    kf <- kalman_filter(model, pre)
    p0 <- kf@filtered_cov[[nrow(kf@filtered_mean)]]
    a <- model@transition[[1L]]
    b <- model@observation[[1L]]
    q <- model@state_cov[[1L]]
    r <- as.numeric(model@obs_cov[[1L]])
    m <- nrow(a)

    ## Simulate the state deviation paths in bulk: one matrix multiply per
    ## post-period step across all replicates.
    n_sim <- 200000L
    p_chol <- t(chol(p0 + diag(1e-12, m)))
    q_chol <- t(chol(q + diag(1e-12, m)))
    dev <- p_chol %*% matrix(stats::rnorm(m * n_sim), nrow = m)
    dev_sum <- matrix(0, nrow = m, ncol = n_sim)
    for (s in seq_len(n_post)) {
      dev <- a %*% dev + q_chol %*% matrix(stats::rnorm(m * n_sim), nrow = m)
      dev_sum <- dev_sum + dev
    }
    v_oracle <- stats::var(as.numeric(b %*% dev_sum)) + n_post * r

    ## 200,000 replicates put the Monte-Carlo standard error near 0.3% of the
    ## variance, so a 1% relative tolerance is wide against the noise yet
    ## tight against the understatement a dropped noise block produces.
    expect_equal(v_fit, v_oracle, tolerance = 0.01)
  }
})

test_that("its_causal() pointwise and cumulative effects are consistent", {
  withr::local_seed(3L)
  pre <- cumsum(stats::rnorm(50, sd = 0.2))
  post <- pre[50L] + cumsum(stats::rnorm(30, sd = 0.2)) + 2
  fit <- its_causal(c(pre, post), intervention = 51L)
  expect_equal(fit@observed - fit@counterfactual, fit@pointwise_effect)
  expect_equal(cumsum(fit@pointwise_effect), fit@cumulative_effect)
  expect_equal(mean(fit@pointwise_effect), fit@avg_effect)
})

test_that("its_causal() trend model carries a pre-period slope forward", {
  ## A pre-period with a clear upward slope and an unchanged post-period should,
  ## under the trend model, forecast a rising counterfactual and so detect a
  ## negative effect (the series fell short of its trajectory).
  withr::local_seed(6L)
  pre <- 0.5 * seq_len(40) + stats::rnorm(40, sd = 1)
  post <- rep(pre[40L], 30L) + stats::rnorm(30, sd = 1)
  fit <- its_causal(c(pre, post), intervention = 41L, trend = TRUE)
  expect_lt(fit@avg_effect, 0)
})

test_that("its_causal() validates its inputs", {
  x <- cumsum(stats::rnorm(60, sd = 0.3))
  expect_error(its_causal(x, intervention = 5L), "ten")
  expect_error(its_causal(x, intervention = 100L), "exceed")
  expect_error(its_causal(x, intervention = 30L, level = 1.5), "level")
  expect_error(its_causal(c(1, NA, 3), intervention = 2L), "missing")
})
