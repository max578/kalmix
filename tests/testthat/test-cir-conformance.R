## Conformance of kalmix's ACI core against aciR (biometryhub/ACI), a second
## implementation of Andreou, Chen and Bollt (2026) whose numerical core is
## graded against the method authors' MATLAB reference. See setup-conformance.R
## for exactly where that oracle's independence lies, and why the sibling kernR
## cross-check does not carry it.

test_that("the causal metric matches aciR::aci_metric on the dyad model", {
  skip_if_not_installed("aciR")
  ## Oracle: aciR::aci_metric(), evaluated on aciR's own dyad filter and
  ## smoother moments, so the only thing under test is kalmix's closed form.
  ## aciR's metric is graded to a maximum absolute error of 4.57e-14 against
  ## the authors' MATLAB on this model (aciR oracle manifest, fixture `dyad`).
  oracle <- .aciR_dyad_oracle(n = 301L, window = 51L)
  metric <- aciR::aci_metric(oracle$filter, oracle$smoother)
  kl <- kalmix:::.gaussian_relative_entropy

  native <- vapply(
    seq_along(metric),
    function(t) {
      kl(
        oracle$smoother$mean[t], matrix(oracle$smoother$cov[t]),
        oracle$filter$mean[t], matrix(oracle$filter$cov[t])
      )
    },
    numeric(1L)
  )
  ## The tolerance absorbs the positive-definiteness jitter kalmix's Cholesky
  ## helper adds to the reference covariance; the measured maximum on this
  ## fixture is of order 1e-8.
  expect_equal(native, metric, tolerance = 1e-6)
})

test_that("the objective range reproduces aciR::aci_cir()'s objective", {
  skip_if_not_installed("aciR")
  ## Oracle: aciR::aci_cir()$objective at each anchor. The divergence profile
  ## handed to kalmix's reducer is computed entirely by aciR's own online
  ## smoother and metric, so the graded quantity is the reduction -- the
  ## functional the audit finding KM-01 is about. The paper normalises the
  ## integral by the PEAK of the profile, which equals D(0) only when the
  ## profile decreases with lag; on this model it does not.
  oracle <- .aciR_dyad_oracle()
  for (i in seq_along(oracle$window)) {
    prof <- oracle$profiles[[i]]
    expect_equal(
      kalmix:::.cir_trapezoid(prof$lag, prof$divergence),
      oracle$range$objective[i],
      tolerance = 5e-3
    )
  }
})

test_that("the subjective range reproduces aciR::aci_cir()'s subjective range", {
  skip_if_not_installed("aciR")
  ## Oracle: aciR::aci_cir()$subjective, the LAST elapsed time at which the
  ## divergence exceeds each threshold (Andreou, Chen and Bollt 2026, eq. 8),
  ## reported on aciR's own 129-point logarithmic threshold grid. This is an
  ## exact identity on a shared profile, not a quadrature, so the tolerance is
  ## at floating-point level.
  oracle <- .aciR_dyad_oracle()
  for (i in seq_along(oracle$window)) {
    prof <- oracle$profiles[[i]]
    expect_equal(
      kalmix:::.cir_subjective_range(
        prof$lag, prof$divergence, oracle$epsilon
      ),
      as.numeric(oracle$range$subjective[, i]),
      tolerance = 1e-10
    )
  }
})

test_that("objective_exact reproduces aciR::aci_cir()'s objective_exact", {
  skip_if_not_installed("aciR")
  ## Oracle: aciR::aci_cir()$objective_exact -- the range obtained by averaging
  ## the subjective ranges over the whole threshold grid, which the source
  ## paper defines and of which the integral form graded above is a
  ## computationally efficient UNDERESTIMATE. The two are different
  ## functionals, equal only on a monotone profile; the dyad profile is not
  ## monotone, so the gap is the finding KM-01 makes visible.
  ## The tolerance is looser than for `objective` because the threshold
  ## quadrature is itself grid-dependent: aciR's default 129-point grid and the
  ## reference's 513-point grid do not give the same number, and neither is
  ## more nearly correct than the other.
  oracle <- .aciR_dyad_oracle()
  for (i in seq_along(oracle$window)) {
    prof <- oracle$profiles[[i]]
    exact <- kalmix:::.cir_reduce(
      prof$lag, prof$divergence,
      functional = "objective_exact", epsilon = oracle$epsilon
    )
    expect_equal(exact$value, oracle$range$objective_exact[i],
                 tolerance = 5e-2)
    ## The efficient form is the underestimate, on every anchor.
    expect_lt(oracle$range$objective[i], oracle$range$objective_exact[i])
    expect_lt(
      kalmix:::.cir_trapezoid(prof$lag, prof$divergence), exact$value
    )
  }
})

test_that("a profile that has not decayed by its last lag is flagged censored", {
  skip_if_not_installed("aciR")
  ## Oracle: aciR's censoring contract -- "a reported time close to the end of
  ## the record cannot be resolved ... its range is right-censored, and the
  ## truncated value is a lower bound" (aci_cir documentation). Truncating a
  ## resolved profile before it decays must flip the flag; the full profile,
  ## which does decay, must not.
  oracle <- .aciR_dyad_oracle(window = 51L)
  prof <- oracle$profiles[[1L]]
  keep <- seq_len(20L)
  short <- kalmix:::.cir_reduce(prof$lag[keep], prof$divergence[keep])
  full <- kalmix:::.cir_reduce(prof$lag, prof$divergence)
  expect_true(short$censored)
  expect_false(full$censored)
  ## A censored range is a lower bound on the resolved one.
  expect_lt(short$value, full$value)
})
