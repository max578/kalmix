test_that("kalmix_filter() has a stable signature", {
  expect_setequal(
    names(formals(kalmix_filter)),
    c("model", "y", "n_regimes", "...")
  )
})

test_that("kalmix_filter() validates its arguments before reaching proxymix", {
  model <- ou_fit(.sim_ou(n = 400, seed = 1L))
  expect_error(kalmix_filter("not a model", 1:10), "spread_model")
  expect_error(kalmix_filter(model, 1), "at least two")
  expect_error(kalmix_filter(model, 1:10, n_regimes = 0L), "positive integer")
})

test_that("kalmix_filter() is dormant until proxymix::gmm_filter() ships", {
  model <- ou_fit(.sim_ou(n = 400, seed = 2L))
  ## proxymix::gmm_filter() does not exist yet, so the adapter must error
  ## with a clear commissioning notice rather than dispatch.
  expect_false(kalmix:::.proxymix_has_filter())
  expect_error(kalmix_filter(model, model@theta + stats::rnorm(50)),
               "proxymix|not available|gmm_filter")
})
