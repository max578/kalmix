## Model-adequacy diagnostics and the grounding / abstention read-out.
##
## The innovation diagnostics are checked against models whose adequacy we know
## by construction (a correctly specified local level versus a frozen-state
## misfit), and the grounding logic is checked as the [unverified] -> abstain
## falsifier: an adequate model is still only self-consistent, a model the data
## reject abstains outright, and only declared, dated mechanism provenance lifts
## an adequate verdict to "grounded".

.local_level <- function(state_cov = 0.04) {
  ssm(
    transition = 1, observation = 1,
    state_cov = state_cov, obs_cov = 1,
    init_state = 0, init_cov = 10
  )
}

.sim_local_level <- function(n, seed) {
  withr::with_seed(seed, {
    level <- cumsum(stats::rnorm(n, sd = 0.2))
    level + stats::rnorm(n)
  })
}

test_that("innovation_diagnostics() passes a correctly specified model", {
  ## The model IS the data-generating process, so at alpha = 0.05 it is rejected
  ## only by chance; the large majority of independent seeds must pass.
  model <- .local_level()
  adequate <- vapply(seq_len(8L), function(s) {
    isTRUE(innovation_diagnostics(
      kalman_filter(model, .sim_local_level(400L, s))
    )@adequate)
  }, logical(1L))
  expect_gte(sum(adequate), 6L)
})

test_that("innovation_diagnostics() rejects a misspecified model", {
  ## A near-frozen state cannot track a random walk: the standardised
  ## innovations are strongly autocorrelated and mis-scaled, so every seed must
  ## be declared inadequate.
  frozen <- .local_level(state_cov = 1e-6)
  inadequate <- vapply(seq_len(8L), function(s) {
    isFALSE(innovation_diagnostics(
      kalman_filter(frozen, .sim_local_level(400L, s))
    )@adequate)
  }, logical(1L))
  expect_equal(sum(inadequate), 8L)
})

test_that("innovation_diagnostics() is indeterminate on too little data", {
  adq <- innovation_diagnostics(
    kalman_filter(.local_level(), .sim_local_level(10L, 1L))
  )
  expect_true(is.na(adq@adequate))
  expect_identical(adq@reason, "insufficient_data")
  expect_true(is.na(adq@p_value))
})

test_that("a correct model's standardised innovations are white and scaled", {
  ## Independent sanity on the internal whitening: under the true model the
  ## standardised innovations have unit variance and negligible autocorrelation.
  w <- kalmix:::.standardised_innovations(
    kalman_filter(.local_level(), .sim_local_level(2000L, 7L))
  )
  expect_equal(stats::var(w[, 1L]), 1, tolerance = 0.1)
  expect_lt(abs(stats::acf(w[, 1L], plot = FALSE)$acf[2L]), 0.1)
})

test_that("aci() abstains under model inadequacy (the falsifier)", {
  ## Known DGP: the coupled OSSE ENSO series. The matched model is adequate and
  ## carries the honest mechanism-unverified label; a frozen-state local level
  ## on the same series is inadequate and ACI abstains outright.
  sim <- .sim_enso_osse(n = 600L, kappa = 0.6, seed = 3L)

  matched <- aci(.enso_osse_model(sim, kappa = 0.6), sim$x, dt = sim$dt)
  expect_true(isTRUE(matched@adequacy@adequate))
  expect_identical(matched@grounding, "[unverified]")
  expect_identical(matched@grounding_reason, "mechanism_unverified")

  frozen <- ssm(
    transition = 1, observation = 1, state_cov = 1e-6, obs_cov = 1,
    init_state = 0, init_cov = 10
  )
  abstained <- aci(frozen, sim$x, dt = sim$dt)
  expect_false(isTRUE(abstained@adequacy@adequate))
  expect_identical(abstained@grounding, "[unverified]")
  expect_identical(abstained@grounding_reason, "model_inadequate")
})

test_that("aci() grounds a verdict only with verified mechanism provenance", {
  sim <- .sim_enso_osse(n = 600L, kappa = 0.6, seed = 5L)
  model <- .enso_osse_model(sim, kappa = 0.6)

  ## No provenance: an adequate model is still only self-consistent.
  expect_identical(aci(model, sim$x, dt = sim$dt)@grounding, "[unverified]")

  ## Declared, dated provenance against an external authority: grounded.
  grounded <- aci(
    model, sim$x, dt = sim$dt,
    mechanism = list(source = "OSSE truth", verified_on = as.Date("2026-06-14"))
  )
  expect_identical(grounded@grounding, "grounded")
  expect_identical(grounded@grounding_reason, "grounded")

  ## A mechanism without a verified_on date does not ground it.
  ungrounded <- aci(
    model, sim$x, dt = sim$dt,
    mechanism = list(source = "claimed but unchecked")
  )
  expect_identical(ungrounded@grounding, "[unverified]")
})

test_that("grounding helpers combine worst-case", {
  expect_identical(
    kalmix:::.combine_grounding("grounded", "grounded"), "grounded"
  )
  expect_identical(
    kalmix:::.combine_grounding("grounded", "[unverified]"), "[unverified]"
  )
  expect_identical(kalmix:::.combine_grounding(), "[unverified]")

  expect_false(kalmix:::.mechanism_grounded(NULL))
  expect_false(kalmix:::.mechanism_grounded(list(source = "x")))
  expect_false(kalmix:::.mechanism_grounded(list(verified_on = NA)))
  expect_true(
    kalmix:::.mechanism_grounded(list(verified_on = as.Date("2026-06-14")))
  )
})

test_that("the adequacy verdict prints", {
  adq <- innovation_diagnostics(
    kalman_filter(.local_level(), .sim_local_level(300L, 2L))
  )
  expect_output(print(adq), "model adequacy")
  expect_output(print(adq), "Ljung-Box")
})
