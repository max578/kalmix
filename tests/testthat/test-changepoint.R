test_that("detect_changepoint() recovers known shifts in the mean", {
  ## Known DGP: mean shifts at indices 50 and 100. The detector must place its
  ## change-points near those known truths.
  withr::local_seed(1L)
  x <- c(stats::rnorm(50, 0), stats::rnorm(50, 4), stats::rnorm(50, 1))
  fit <- detect_changepoint(x)
  expect_s3_class(fit, "kalmix::changepoint_fit")
  expect_equal(length(fit@changepoints), 2L)
  expect_true(all(abs(fit@changepoints - c(50L, 100L)) <= 3L))
  expect_equal(fit@n_segments, 3L)
  expect_output(print(fit), "change-point")
})

test_that("detect_changepoint() reports no change on a homogeneous series", {
  ## A single homogeneous Gaussian segment must not be split under the default
  ## penalty -- the over-segmentation guard that the penalty exists to enforce.
  withr::local_seed(2L)
  x <- stats::rnorm(200)
  fit <- detect_changepoint(x)
  expect_length(fit@changepoints, 0L)
  expect_equal(fit@n_segments, 1L)
})

test_that("detect_changepoint() segment means match the truth", {
  ## The fitted per-segment means must recover the known segment means.
  withr::local_seed(3L)
  x <- c(stats::rnorm(60, -2), stats::rnorm(60, 5))
  fit <- detect_changepoint(x)
  expect_length(fit@changepoints, 1L)
  expect_equal(fit@segment_means[1L], -2, tolerance = 0.4)
  expect_equal(fit@segment_means[2L], 5, tolerance = 0.4)
})

test_that("detect_changepoint() penalty controls sensitivity", {
  ## A weaker penalty admits more change-points than a stronger one on the same
  ## series; this monotonicity is the contract of the penalty argument.
  withr::local_seed(4L)
  x <- c(stats::rnorm(40, 0), stats::rnorm(40, 1.5), stats::rnorm(40, 0))
  n_mbic <- length(detect_changepoint(x, penalty = "mbic")@changepoints)
  n_sic <- length(detect_changepoint(x, penalty = "sic")@changepoints)
  expect_gte(n_sic, n_mbic)
})

test_that("detect_changepoint() validates its inputs", {
  expect_error(detect_changepoint(matrix(1, 5, 2)), "vector")
  expect_error(detect_changepoint(c(1, NA, 3, 4, 5, 6, 7, 8, 9, 10)), "missing")
  expect_error(detect_changepoint(stats::rnorm(50), penalty = "nope"), "mbic")
  expect_error(detect_changepoint(stats::rnorm(50), penalty = -1), "negative")
  expect_error(detect_changepoint(1:4, min_segment = 5L), "observations")
})
