test_that("hmm() builds and validates a Gaussian HMM", {
  m <- hmm(
    init_prob = c(0.5, 0.5),
    transition = matrix(c(0.9, 0.1, 0.1, 0.9), nrow = 2L, byrow = TRUE),
    emission_mean = c(0, 3), emission_sd = c(1, 1)
  )
  expect_s3_class(m, "kalmix::hmm")
  expect_equal(m@n_states, 2L)
  expect_output(print(m), "hidden Markov")

  expect_error(
    hmm(c(0.4, 0.4), matrix(c(0.9, 0.1, 0.1, 0.9), 2, byrow = TRUE),
        c(0, 3), c(1, 1)),
    "sum to one"
  )
  expect_error(
    hmm(c(0.5, 0.5), matrix(c(0.9, 0.1, 0.1, 0.9), 2, byrow = TRUE),
        c(0, 3), c(1, -1)),
    "positive"
  )
})

test_that("hmm_filter() and hmm_viterbi() recover a known state path", {
  ## Known DGP: a simulated two-state mean-shift sequence with the true state
  ## path returned. Recovery is measured against that truth, not the code.
  sim <- .sim_hmm(n = 200, means = c(0, 4), sd = 1, stay = 0.95, seed = 1L)
  model <- hmm(
    init_prob = c(0.5, 0.5),
    transition = sim$transition,
    emission_mean = sim$means, emission_sd = c(sim$sd, sim$sd)
  )
  fit <- hmm_filter(model, sim$y)
  expect_s3_class(fit, "kalmix::hmm_fit")
  expect_equal(rowSums(fit@filtered), rep(1, 200L), tolerance = 1e-8)
  expect_equal(rowSums(fit@smoothed), rep(1, 200L), tolerance = 1e-8)

  ## Well-separated states (means 0 and 4, sd 1) should be recovered almost
  ## perfectly by both the smoothed mode and the Viterbi path.
  expect_gt(mean(fit@state_path == sim$state), 0.9)
  vp <- hmm_viterbi(model, sim$y)
  expect_gt(mean(vp == sim$state), 0.9)
})

test_that("hmm_filter() log-likelihood matches the mixture marginal at p=0", {
  ## With a frozen chain (no transitions, equal start) and identical emissions,
  ## the HMM log-likelihood reduces to the i.i.d. Gaussian-mixture marginal,
  ## an oracle independent of the forward recursion.
  withr::local_seed(2L)
  y <- stats::rnorm(40, mean = 1, sd = 1)
  model <- hmm(
    init_prob = c(0.5, 0.5),
    transition = matrix(c(0.5, 0.5, 0.5, 0.5), nrow = 2L),
    emission_mean = c(1, 1), emission_sd = c(1, 1)
  )
  fit <- hmm_filter(model, y)
  direct <- sum(stats::dnorm(y, mean = 1, sd = 1, log = TRUE))
  expect_equal(fit@log_lik, direct, tolerance = 1e-8)
})

test_that("hmm_viterbi() is jointly optimal on a hand worked example", {
  ## A tiny example with an obvious optimal path: two clearly separated points
  ## must decode to states 1 then 2.
  model <- hmm(
    init_prob = c(0.5, 0.5),
    transition = matrix(c(0.5, 0.5, 0.5, 0.5), nrow = 2L),
    emission_mean = c(0, 10), emission_sd = c(1, 1)
  )
  expect_equal(hmm_viterbi(model, c(0, 10)), c(1L, 2L))
  expect_equal(hmm_viterbi(model, c(10, 0)), c(2L, 1L))
})

test_that("hmm verbs validate observations", {
  model <- hmm(
    init_prob = c(0.5, 0.5),
    transition = matrix(c(0.9, 0.1, 0.1, 0.9), nrow = 2L, byrow = TRUE),
    emission_mean = c(0, 3), emission_sd = c(1, 1)
  )
  expect_error(hmm_filter("x", 1:10), "hmm")
  expect_error(hmm_filter(model, 1), "at least two")
  ## An NA is a missing observation (flat emission), not an input error.
  expect_no_error(hmm_viterbi(model, c(1, NA, 3)))
})
