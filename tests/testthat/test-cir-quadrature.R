## The quadrature grid and the censoring contract of the causal influence
## range. The oracles here are arithmetic rather than external: the range is an
## integral of the divergence profile, so evaluating that profile at every lag
## the window contains is the exact quadrature and a thinned grid is an
## approximation to it, and an integral over the lags the record supports
## cannot exceed the largest of those lags.

.cir_quadrature_case <- function(n = 600L, seed = 5L) {
  sim <- .sim_enso_osse(n = n, kappa = 0.6, seed = seed)
  list(sim = sim, model = .enso_osse_model(sim, kappa = 0.6))
}

test_that("the default lag grid resolves the window, a 12-point grid does not", {
  case <- .cir_quadrature_case()
  max_lag <- 120L

  ## Oracle: the fully resolved quadrature, every integer lag in the window.
  dense <- suppressWarnings(causal_information_rate(
    case$model, case$sim$x, max_lag = max_lag, n_lag = max_lag + 1L
  ))
  coarse <- suppressWarnings(causal_information_rate(
    case$model, case$sim$x, max_lag = max_lag, n_lag = 12L
  ))
  default <- suppressWarnings(causal_information_rate(
    case$model, case$sim$x, max_lag = max_lag
  ))

  ## A twelve-point trapezoid over a convex decaying profile over-estimates the
  ## integral, and does so by far more than the tolerance a lead-time reported
  ## to two significant figures can carry.
  expect_gt(abs(coarse - dense) / dense, 0.05)
  ## The default must land on the resolved value, not on the coarse one.
  expect_lt(abs(default - dense) / dense, 0.02)
  expect_equal(attr(default, "n_lag"), max_lag + 1L)
  expect_equal(attr(dense, "n_lag"), max_lag + 1L)
})

test_that("a grid too coarse to have converged warns", {
  case <- .cir_quadrature_case()
  expect_warning(
    causal_information_rate(
      case$model, case$sim$x, max_lag = 120L, n_lag = 12L
    ),
    "converged"
  )
  ## A grid that already carries every lag has nothing finer to converge to,
  ## so the check is silent there.
  expect_no_warning(
    causal_information_rate(
      case$model, case$sim$x, max_lag = 120L, n_lag = 121L
    )
  )
})

test_that("aci() exposes the quadrature grid so the lead-time can be refined", {
  case <- .cir_quadrature_case()
  max_lag <- 120L
  dense <- suppressWarnings(causal_information_rate(
    case$model, case$sim$x, max_lag = max_lag, n_lag = max_lag + 1L
  ))
  fit <- suppressWarnings(aci(
    case$model, case$sim$x, dt = case$sim$dt,
    max_lag = max_lag, n_lag = max_lag + 1L
  ))
  expect_equal(fit@lead_time, as.numeric(dense) * case$sim$dt,
               tolerance = 1e-8)
})

test_that("an unresolvable anchor is censored rather than extrapolated", {
  case <- .cir_quadrature_case()
  ## Anchor 580 of 600 has twenty future steps, so a 120-step window cannot be
  ## closed. Oracle: the objective range is the integral of a non-negative
  ## profile over the lags the record supports, normalised by that profile's
  ## peak, so it is bounded above by the largest resolved lag -- twenty steps.
  ## Carrying the last divergence forward across the missing hundred lags
  ## breaks that bound.
  res <- suppressWarnings(causal_information_rate(
    case$model, case$sim$x, eval_points = 580L, max_lag = 120L
  ))
  expect_true(attr(res, "censored"))
  expect_lte(as.numeric(res), 20)

  ## An interior anchor with the same window is fully resolved.
  ok <- suppressWarnings(causal_information_rate(
    case$model, case$sim$x, eval_points = 300L, max_lag = 120L
  ))
  expect_false(attr(ok, "censored"))
})

test_that("a censored lead-time cannot carry a grounded verdict", {
  case <- .cir_quadrature_case()
  mechanism <- list(
    verified_on = as.Date("2026-01-01"),
    source = "the simulated data-generating process"
  )
  fit <- suppressWarnings(aci(
    case$model, case$sim$x, dt = case$sim$dt,
    eval_points = 580L, max_lag = 120L, mechanism = mechanism
  ))
  expect_true(fit@censored)
  expect_identical(fit@grounding, "[unverified]")
  expect_identical(fit@grounding_reason, "censored_horizon")
  expect_output(print(fit), "censored")
})

test_that("an unconverged lead-time cannot carry a grounded verdict", {
  ## KM-D1 (orchestra fitness audit, 2026-08-25): the convergence check
  ## (test above, "a grid too coarse to have converged warns") only warned;
  ## it never touched aci()'s grounding token, so a caller reading `grounding`
  ## alone -- rather than parsing warnings -- saw "grounded" on a lead-time
  ## known not to have converged. This mirrors the censoring degrade already
  ## in place and rides on the same `converged` attribute wave 2 added to
  ## causal_information_rate().
  case <- .cir_quadrature_case()
  mechanism <- list(
    verified_on = as.Date("2026-01-01"),
    source = "the simulated data-generating process"
  )
  fit <- suppressWarnings(aci(
    case$model, case$sim$x, dt = case$sim$dt,
    max_lag = 120L, n_lag = 12L, mechanism = mechanism
  ))
  expect_false(isTRUE(attr(
    suppressWarnings(causal_information_rate(
      case$model, case$sim$x, max_lag = 120L, n_lag = 12L
    )),
    "converged"
  )))
  expect_identical(fit@grounding, "[unverified]")
  expect_identical(fit@grounding_reason, "not_converged")
})

test_that("an unconverged aci_fit carries the orchestra abstention contract classes", {
  case <- .cir_quadrature_case()
  fit <- suppressWarnings(aci(
    case$model, case$sim$x, dt = case$sim$dt,
    max_lag = 120L, n_lag = 12L
  ))
  expect_true(inherits(fit, "kalmix_abstention"))
  expect_true(inherits(fit, "orchestra_refusal"))

  ## A fully resolved grid has nothing to abstain on for this reason: it is
  ## not stamped with the convergence-abstention class.
  resolved <- suppressWarnings(aci(
    case$model, case$sim$x, dt = case$sim$dt,
    max_lag = 120L, n_lag = 121L
  ))
  expect_false(inherits(resolved, "kalmix_abstention"))
})

test_that("the exact and efficient functionals bracket a non-monotone profile", {
  ## Oracle: the layer-cake identity. Integrating the divergence profile and
  ## integrating its superlevel-set measure over the tolerance give the same
  ## number; the published subjective range takes the LAST crossing rather than
  ## the measure of the superlevel set, so the exact functional is at least the
  ## efficient one, with equality exactly when the profile is monotone.
  lag <- 0:20
  ## A monotone profile: the two functionals coincide.
  mono <- exp(-lag / 5)
  eps <- 10^seq(-6, 0.5, length.out = 129L)
  a <- kalmix:::.cir_reduce(lag, mono, functional = "objective")
  b <- kalmix:::.cir_reduce(lag, mono, functional = "objective_exact",
                            epsilon = eps)
  expect_true(a$monotone)
  ## Equal up to the one-grid-step advance in the subjective range's
  ## definition, which is all the discretisation of a continuous identity can
  ## promise.
  expect_gte(b$value, a$value)
  expect_lt(b$value - a$value, 1)

  ## A profile that rises before it decays: the efficient form is strictly
  ## smaller, and the monotonicity flag says so.
  bump <- mono
  bump[2:4] <- bump[2:4] + 0.6
  c2 <- kalmix:::.cir_reduce(lag, bump, functional = "objective")
  d2 <- kalmix:::.cir_reduce(lag, bump, functional = "objective_exact",
                             epsilon = eps)
  expect_false(c2$monotone)
  expect_lt(c2$value, d2$value)
})
