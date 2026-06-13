test_that("ssm() builds from scalars and reports its dimensions", {
  m <- ssm(
    transition = 1, observation = 1, state_cov = 0.1, obs_cov = 1,
    init_state = 0, init_cov = 10
  )
  expect_s3_class(m, "kalmix::ssm")
  expect_equal(m@state_dim, 1L)
  expect_equal(m@obs_dim, 1L)
  expect_length(m@transition, 1L)
  expect_output(print(m), "linear-Gaussian")
})

test_that("ssm() builds a vector-state model and promotes diagonal covs", {
  m <- ssm(
    transition = matrix(c(1, 0, 1, 1), nrow = 2L),
    observation = matrix(c(1, 0), nrow = 1L),
    state_cov = c(0.01, 0.02),
    obs_cov = 1,
    init_state = c(0, 0),
    init_cov = c(10, 10)
  )
  expect_equal(m@state_dim, 2L)
  expect_equal(m@obs_dim, 1L)
  expect_equal(diag(m@state_cov[[1L]]), c(0.01, 0.02))
})

test_that("ssm() accepts time-varying system matrices", {
  a_list <- replicate(5L, matrix(1), simplify = FALSE)
  m <- ssm(
    transition = a_list, observation = 1, state_cov = 0.1, obs_cov = 1,
    init_state = 0, init_cov = 10
  )
  expect_length(m@transition, 5L)
})

test_that("ssm() rejects dimension mismatches", {
  expect_error(
    ssm(
      transition = matrix(1, 2, 2), observation = 1, state_cov = 0.1,
      obs_cov = 1, init_state = 0, init_cov = 10
    ),
    "transition"
  )
  expect_error(
    ssm(
      transition = 1, observation = matrix(1, 1, 2), state_cov = 0.1,
      obs_cov = 1, init_state = 0, init_cov = 10
    ),
    "observation"
  )
})
