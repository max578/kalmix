test_that("spread_series() recovers the hedge ratio", {
  pair <- .sim_pair(n = 600, beta = 2, seed = 1L)
  s <- spread_series(pair$a, pair$b)
  expect_length(s, 600L)
  expect_equal(attr(s, "hedge_ratio"), pair$beta, tolerance = 0.05)
})

test_that("spread_series() validates its inputs", {
  expect_error(spread_series(1:5, 1:6), "same length")
  expect_error(spread_series("a", "b"), "numeric")
  expect_error(spread_series(c(1, NA, 3), c(1, 2, 3)), "missing")
})

test_that("ou_fit() recovers Ornstein-Uhlenbeck parameters", {
  x <- .sim_ou(n = 2000, ar = 0.8, seed = 2L)
  model <- ou_fit(x)
  expect_s3_class(model, "kalmix::spread_model")
  ## AR coef 0.8 -> kappa = -log(0.8) ~ 0.223 (finite-sample band)
  expect_equal(model@kappa, -log(0.8), tolerance = 0.15)
  expect_gt(model@sigma, 0)
  expect_gt(half_life(model), 0)
  expect_equal(model@n_obs, 2000L)
})

test_that("ou_fit() refuses a series inconsistent with Ornstein-Uhlenbeck", {
  ## An anti-persistent series has an AR(1) coefficient below zero, outside
  ## the (0, 1) range an OU process implies.
  anti <- withr::with_seed(7L, as.numeric(stats::arima.sim(list(ar = -0.6), n = 500)))
  expect_error(ou_fit(anti), "Ornstein-Uhlenbeck")
})

test_that("ou_fit()'s refusal carries the orchestra refusal contract classes", {
  ## KM-D1 (orchestra fitness audit, 2026-08-25): ou_fit() was the one verb of
  ## thirty that declined with a bare stop()-shaped condition -- an error with
  ## no class a cross-member caller could recognise as a structured refusal.
  ## The refusal contract (integration/refusal_contract.R) reads any class
  ## ending in `_refusal`, or the explicit `orchestra_refusal` marker, as a
  ## member declining to produce a result at all.
  anti <- withr::with_seed(7L, as.numeric(stats::arima.sim(list(ar = -0.6), n = 500)))
  cnd <- tryCatch(ou_fit(anti), error = function(e) e)
  expect_s3_class(cnd, "error")
  expect_true(inherits(cnd, "kalmix_refusal"))
  expect_true(inherits(cnd, "orchestra_refusal"))
  expect_error(ou_fit(anti), class = "kalmix_refusal")
})

test_that("ou_fit() carries the hedge ratio through from spread_series()", {
  pair <- .sim_pair(n = 800, beta = 1.5, seed = 4L)
  model <- ou_fit(spread_series(pair$a, pair$b))
  expect_equal(model@hedge_ratio, 1.5, tolerance = 0.1)
})

test_that("spread_model prints", {
  model <- ou_fit(.sim_ou(n = 500, seed = 5L))
  expect_output(print(model), "Ornstein-Uhlenbeck")
  expect_output(print(model), "half-life")
})
