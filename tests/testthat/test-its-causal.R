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
  ## the maximum-likelihood noise split is what makes it hold.
  withr::local_seed(5L)
  covered <- vapply(seq_len(120L), function(i) {
    x <- cumsum(stats::rnorm(120, sd = 0.3))
    f <- its_causal(x, intervention = 61L)
    f@avg_effect_lower <= 0 && f@avg_effect_upper >= 0
  }, logical(1L))
  expect_gt(mean(covered), 0.85)
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
