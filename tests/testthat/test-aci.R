test_that(".gaussian_relative_entropy() matches analytic one-dimensional KL", {
  ## Independent oracle: the KL between two univariate Gaussians has the
  ## elementary closed forms below, computed here from first principles rather
  ## than from the function under test. The tolerance absorbs the small
  ## positive-definiteness jitter the Cholesky helper adds to the reference
  ## covariance.
  kl <- kalmix:::.gaussian_relative_entropy

  ## Pure mean shift of one standard deviation: KL(N(1,1) || N(0,1)) = 0.5.
  expect_equal(kl(1, matrix(1), 0, matrix(1)), 0.5, tolerance = 1e-6)

  ## Pure variance change: KL(N(0,2) || N(0,1)) = 0.5 (s2 - 1 - log s2) with
  ## s2 = 2.
  expect_equal(
    kl(0, matrix(2), 0, matrix(1)),
    0.5 * (2 - 1 - log(2)),
    tolerance = 1e-6
  )

  ## Identical distributions give exactly zero.
  expect_equal(kl(c(2, -3), diag(2), c(2, -3), diag(2)), 0)

  ## Non-negativity and asymmetry of the divergence.
  a <- kl(c(0, 0), diag(2), c(1, 0), diag(c(2, 1)))
  b <- kl(c(1, 0), diag(c(2, 1)), c(0, 0), diag(2))
  expect_gt(a, 0)
  expect_gt(b, 0)
  expect_false(isTRUE(all.equal(a, b)))
})

test_that(".gaussian_relative_entropy() agrees with kernR::relative_entropy()", {
  ## Independent-oracle cross-check (charter Independent Oracle Principle): the
  ## native closed form is verified against kernR's relative_entropy(), an
  ## implementation kalmix did not author. Skipped when kernR is absent, so the
  ## package stays standalone-functional and checks green without it.
  skip_if_not_installed("kernR")
  kl <- kalmix:::.gaussian_relative_entropy

  withr::local_seed(20260613L)
  for (i in seq_len(20L)) {
    k <- sample(1:3, 1L)
    mu_p <- stats::rnorm(k)
    mu_q <- stats::rnorm(k)
    ## Random symmetric positive-definite covariances.
    rand_spd <- function() {
      m <- matrix(stats::rnorm(k * k), k, k)
      crossprod(m) + diag(k)
    }
    sp <- rand_spd()
    sq <- rand_spd()
    native <- kl(mu_p, sp, mu_q, sq)
    oracle <- kernR::relative_entropy(mu_p, sp, mu_q, sq)
    ## Two independent implementations of the same closed form: agreement to
    ## floating-point accumulation, not bit-identity.
    expect_equal(native, oracle, tolerance = 1e-6)
  }
})

test_that("aci() returns a well-formed aci_fit", {
  sim <- .sim_enso_osse(n = 300L, kappa = 0.6, seed = 1L)
  model <- .enso_osse_model(sim, kappa = 0.6)
  fit <- aci(model, sim$x, dt = sim$dt)

  expect_s3_class(fit, "kalmix::aci_fit")
  expect_length(fit@causal_information, 300L)
  expect_true(all(fit@causal_information >= 0))
  expect_equal(fit@mean_causal_information, mean(fit@causal_information))
  expect_s3_class(fit@filter, "kalmix::kalman_fit")
  expect_s3_class(fit@smoother, "kalmix::rts_fit")
  expect_true(is.finite(fit@lead_time))
  expect_output(print(fit), "assimilative causal inference")
  expect_output(print(fit), "lead-time")
})

test_that("aci() lead_time = FALSE skips the rate and prints honestly", {
  sim <- .sim_enso_osse(n = 250L, kappa = 0.6, seed = 2L)
  fit <- aci(.enso_osse_model(sim, 0.6), sim$x, lead_time = FALSE)
  expect_true(is.na(fit@causal_information_rate))
  expect_true(is.na(fit@lead_time))
  expect_output(print(fit), "not computed")
})

test_that("OSSE: causal information is zero under the null, positive when coupled", {
  ## Known DGP. Under the null (kappa = 0) the observation carries no
  ## information about the hidden cause, so the smoother cannot improve on the
  ## filter and the causal information must be (numerically) zero. With a real
  ## coupling the future of the series sharpens the hidden state, so the causal
  ## information must be substantially positive. The expected behaviour comes
  ## from the simulated truth, not from the code.
  null_ci <- vapply(seq_len(8L), function(s) {
    sim <- .sim_enso_osse(n = 500L, kappa = 0, seed = s)
    aci(.enso_osse_model(sim, kappa = 0), sim$x,
        lead_time = FALSE)@mean_causal_information
  }, numeric(1L))

  coupled_ci <- vapply(seq_len(8L), function(s) {
    sim <- .sim_enso_osse(n = 500L, kappa = 0.6, seed = s)
    aci(.enso_osse_model(sim, kappa = 0.6), sim$x,
        lead_time = FALSE)@mean_causal_information
  }, numeric(1L))

  expect_lt(max(null_ci), 1e-6)
  expect_gt(min(coupled_ci), 0.3)
  expect_gt(mean(coupled_ci), mean(null_ci) + 0.3)
})

test_that("OSSE: the smoother recovers the hidden ENSO state from the effect", {
  ## Known DGP. The inverse-problem essence of ACI: from the observed effect
  ## series alone the smoother must reconstruct the hidden SST-anomaly state it
  ## was simulated from. The true state is used only to score, never to fit.
  sim <- .sim_enso_osse(n = 600L, kappa = 0.6, seed = 77L)
  fit <- aci(.enso_osse_model(sim, kappa = 0.6), sim$x, lead_time = FALSE)
  recovery <- stats::cor(fit@smoother@smoothed_mean[, 1L], sim$z[, 1L])
  expect_gt(recovery, 0.8)
})

test_that("OSSE: the lead-time is finite, positive and rises with noise", {
  ## Known DGP. The decision lead-time must be a finite, positive horizon for a
  ## genuinely coupled system, and -- the qualitative oracle -- it must grow as
  ## the observation channel gets noisier: a noisier effect needs more of the
  ## future assimilated before the hidden state is pinned, so the cause keeps
  ## informing the estimate over a longer window.
  lead_at <- function(obs_sd, seed) {
    sim <- .sim_enso_osse(n = 600L, kappa = 0.6, obs_sd = obs_sd, seed = seed)
    aci(.enso_osse_model(sim, kappa = 0.6), sim$x, dt = sim$dt)@lead_time
  }
  lead_low <- lead_at(0.3, 5L)
  lead_high <- lead_at(1.2, 5L)

  ## A sensible ENSO predictability horizon: positive, finite, under a few
  ## years even in the noisy case.
  expect_gt(lead_low, 0)
  expect_lt(lead_high, 5)
  expect_gt(lead_high, lead_low)
})

test_that("OSSE: the null lead-time collapses to zero", {
  ## Known DGP. With no coupling there is no recoverable future information, so
  ## the objective causal information rate must report a zero lead-time rather
  ## than a spurious horizon from a flat, near-zero divergence profile.
  sim <- .sim_enso_osse(n = 500L, kappa = 0, seed = 11L)
  fit <- aci(.enso_osse_model(sim, kappa = 0), sim$x, dt = sim$dt)
  expect_equal(fit@lead_time, 0)
})

test_that("causal_information_rate() validates its inputs", {
  sim <- .sim_enso_osse(n = 200L, kappa = 0.6, seed = 3L)
  model <- .enso_osse_model(sim, 0.6)
  expect_error(causal_information_rate("nope", sim$x), "ssm")
  expect_error(causal_information_rate(model, sim$x, dt = -1), "positive")
  expect_error(
    causal_information_rate(model, sim$x, eval_points = c(1L, 9999L)),
    "seq_len"
  )
})

test_that("causal_information_rate() handles time-varying models online only", {
  ## The online engine assimilates forward in place, so a time-varying model is
  ## fine; the expanding engine reruns the recursions on truncated series, which
  ## a time-varying model cannot honour, so it must say so rather than fail
  ## obscurely. The causal information series itself is unaffected either way.
  n <- 60L
  a_list <- rep(list(matrix(1)), n)
  model <- ssm(
    transition = a_list, observation = 1, state_cov = 0.1, obs_cov = 1,
    init_state = 0, init_cov = 10
  )
  withr::local_seed(4L)
  y <- cumsum(stats::rnorm(n, sd = 0.3)) + stats::rnorm(n)
  cir <- causal_information_rate(model, y)
  expect_true(is.finite(cir) && cir >= 0)
  expect_error(
    causal_information_rate(model, y, engine = "expanding"), "time-invariant"
  )
  ## The series path still works.
  fit <- aci(model, y, lead_time = FALSE)
  expect_length(fit@causal_information, n)
})

test_that("real ENSO index: ACI runs on the NOAA ONI as an external oracle", {
  ## Independent real-data oracle (charter Independent Oracle Principle). The
  ## OSSE above proves self-consistency against a known DGP; this pass confirms
  ## the read-out behaves on a real ENSO index nobody in this package authored.
  ## A local-linear-oscillation state-space model is fitted to the recorded
  ## NOAA Oceanic Nino Index and the ACI read-out is computed; the assertions
  ## are sanity bounds the physics implies (a positive, multi-month but
  ## sub-decadal predictability horizon, and a smoother that tracks the index),
  ## not point values. Gated so offline and CRAN runs skip it.
  skip_on_cran()
  skip_if_not_installed("curl")
  skip_if_offline()

  fixture <- testthat::test_path("fixtures", "oni_noaa.csv")
  skip_if_not(file.exists(fixture), "ONI fixture not present")
  oni <- utils::read.csv(fixture, stringsAsFactors = FALSE)$oni
  expect_gt(length(oni), 100L)

  ## A two-state local oscillator: a damped rotation captures the quasi-periodic
  ## ENSO swing, observed through the index with measurement noise. The phases
  ## are unobserved; only the index is seen.
  dt_year <- 1 / 12
  angle <- 2 * pi * dt_year / 3.8           # ~3.8-year ENSO period
  rho <- 0.96                                # weak damping
  rot <- rho * matrix(
    c(cos(angle), -sin(angle), sin(angle), cos(angle)),
    nrow = 2L, byrow = TRUE
  )
  model <- ssm(
    transition = rot,
    observation = matrix(c(1, 0), nrow = 1L),
    state_cov = diag(c(0.05, 0.05)),
    obs_cov = 0.2,
    init_state = c(oni[1L], 0),
    init_cov = diag(c(1, 1))
  )

  fit <- aci(model, oni, dt = dt_year)

  ## The read-out must be well-formed and physically sensible.
  expect_true(all(fit@causal_information >= 0))
  expect_gt(fit@mean_causal_information, 0)
  expect_true(is.finite(fit@lead_time))
  expect_gt(fit@lead_time, 0)
  expect_lt(fit@lead_time, 5)

  ## The smoother should track the observed index (the observation maps the
  ## first state directly), an honest sanity check on the real-data pass.
  tracked <- stats::cor(fit@smoother@smoothed_mean[, 1L], oni)
  expect_gt(tracked, 0.7)
})
