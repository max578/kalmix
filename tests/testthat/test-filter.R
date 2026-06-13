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
