test_that("kalmix_filter() has a stable signature", {
  expect_setequal(
    names(formals(kalmix_filter)),
    c("model", "y", "n_regimes", "vol_multipliers", "engine")
  )
})

test_that("kalmix_filter() runs the native regime engine standalone", {
  ## The native engine must work with no sibling package installed. Filtering a
  ## fitted spread returns a regime_fit over its volatility regimes.
  model <- ou_fit(.sim_ou(n = 400, seed = 1L))
  fit <- kalmix_filter(model, model@theta + .sim_ou(n = 400, seed = 1L))
  expect_s3_class(fit, "kalmix::regime_fit")
  expect_equal(fit@n_regimes, 2L)
  expect_equal(rowSums(fit@regime_prob), rep(1, nrow(fit@regime_prob)),
               tolerance = 1e-8)
})

test_that("kalmix_filter() tracks a known volatility regime in the spread", {
  ## Known DGP: a mean-reverting spread that is calm then turbulent. The filter
  ## should raise the turbulent-regime probability in the second half.
  withr::local_seed(2L)
  calm <- as.numeric(stats::arima.sim(list(ar = 0.7), n = 150, sd = 0.4))
  wild <- as.numeric(stats::arima.sim(list(ar = 0.7), n = 150, sd = 2))
  spread <- c(calm, wild)
  model <- ou_fit(spread)
  fit <- kalmix_filter(model, spread, n_regimes = 2L)
  expect_gt(
    mean(fit@regime_prob[151:300, 2L]),
    mean(fit@regime_prob[1:150, 2L])
  )
})

test_that("kalmix_filter() tracks a spread with a non-zero long-run mean", {
  ## Known DGP: a mean-reverting spread centred at 5 rather than zero. The
  ## regime models operate on deviations from the fitted long-run mean, so the
  ## filtered state must track the spread around theta instead of being pulled
  ## toward zero by the autoregressive contraction.
  x <- 5 + .sim_ou(n = 2000, seed = 7L)
  model <- ou_fit(x)
  fit <- kalmix_filter(model, x)
  expect_equal(mean(fit@state_mean), model@theta, tolerance = 0.02)
  expect_gt(stats::cor(as.numeric(fit@state_mean), x), 0.9)

  ## The single-regime path restores the same observed scale.
  fit_one <- kalmix_filter(model, x, n_regimes = 1L)
  expect_equal(mean(fit_one@filtered_mean), model@theta, tolerance = 0.02)

  ## Regime tracking still works on the shifted spread: calm-then-turbulent
  ## noise must raise the turbulent-regime probability in the second half.
  withr::local_seed(8L)
  calm <- as.numeric(stats::arima.sim(list(ar = 0.7), n = 150, sd = 0.4))
  wild <- as.numeric(stats::arima.sim(list(ar = 0.7), n = 150, sd = 2))
  spread <- 5 + c(calm, wild)
  model_shift <- ou_fit(spread)
  fit_shift <- kalmix_filter(model_shift, spread, n_regimes = 2L)
  expect_gt(
    mean(fit_shift@regime_prob[151:300, 2L]),
    mean(fit_shift@regime_prob[1:150, 2L])
  )
})

test_that("kalmix_filter() with one regime is a plain Kalman filter", {
  model <- ou_fit(.sim_ou(n = 300, seed = 3L))
  fit <- kalmix_filter(model, .sim_ou(n = 300, seed = 3L), n_regimes = 1L)
  expect_s3_class(fit, "kalmix::kalman_fit")
})

test_that("kalmix_filter() validates its arguments", {
  model <- ou_fit(.sim_ou(n = 400, seed = 4L))
  expect_error(kalmix_filter("not a model", 1:10), "spread_model")
  expect_error(kalmix_filter(model, 1), "at least two")
  expect_error(kalmix_filter(model, 1:10, n_regimes = 0L), "positive integer")
  expect_error(
    kalmix_filter(model, 1:10, n_regimes = 2L, vol_multipliers = c(1, -1)),
    "positive"
  )
})

test_that("the proxymix mixture-Kalman primitive is not relied upon", {
  ## kalmix is standalone-functional: the native filter is the default, so the
  ## absent proxymix primitive must never block a filter call.
  expect_false(kalmix:::.proxymix_has_filter())
  model <- ou_fit(.sim_ou(n = 200, seed = 5L))
  expect_no_error(kalmix_filter(model, .sim_ou(n = 200, seed = 5L)))
})
