## Independent-oracle fixtures for the ACI conformance suite.
##
## acir (biometryhub/ACI, public; renamed from aciR on 2026-08-28) is a second
## implementation of Andreou, Chen and Bollt (2026) whose numerical core is
## graded cell-by-cell against the method authors' own MATLAB reference
## (Andreou, github.com/marandmath/ACI_code, commit 733c49f, MIT), with the
## reference fixtures, their provenance and the measured agreement shipped in
## `system.file("evidence", "register.csv", package = "acir")`.
##
## Be exact about where the independence lies. acir shares a maintainer with
## kalmix, so agreement between the two R packages by itself would establish
## consistency and nothing more. What makes acir an oracle is the chain its
## numbers carry back to the authors' own code: the values graded below are
## values kalmix did not produce, computed by an implementation that must
## reproduce the reference's MATLAB output to a published tolerance before it
## ships. The sibling `kernR::relative_entropy()` cross-check has no such
## chain -- it is the same textbook closed form in the same estate -- so it
## grades transcription, not correspondence.
##
## acir works in continuous time on conditional Gaussian nonlinear systems and
## kalmix in discrete time on state-space models, so the two filters are
## different objects that agree only as `dt` tends to zero. What IS directly
## comparable, and what the audit findings KM-01 to KM-03 are about, is the
## reduction of a divergence-versus-lag profile to a causal influence range.
## The fixtures below therefore build the divergence profile with acir's OWN
## online smoother and metric, and hand it to kalmix's reducer.

## Build the dyad divergence profiles and acir's reported ranges.
##
## The nonlinear dyad model with intermittent extreme events is the model the
## authors' `dyad_interaction_model.m` implements and the one acir's `dyad`
## fixtures grade. The rename replaced the closure-based `aci_dyad_components()`
## step with a model object that carries its own coefficients, and split the
## old `aci_cir()` into an `aci_range()` call per functional: `method =
## "l1_linf"` is the old `objective`, `method = "exact"` with `quadrature =
## "matlab_eps_grid"` is the old `objective_exact` on the same threshold grid.
## `aci_online_smoother()` is now `aci_online()`.
.acir_dyad_oracle <- function(n = 601L,
                              dt = 0.001,
                              seed = 11L,
                              window = c(51L, 201L)) {
  model <- acir::aci_dyad_model()
  sim <- acir::aci_simulate(model, seed = seed, t_end = (n - 1) * dt, dt = dt)
  filt <- acir::aci_filter(
    model, sim, init = list(mean = model$meta$ic_default$y0, cov = 0.1)
  )
  smooth <- acir::aci_smoother(model, sim, filter = filt)
  epsilon <- 10^seq(-6, 0.5, length.out = 129L)
  table <- acir::lag_table(model, sim, mode = "forward", filter = filt)
  range <- acir::aci_range(
    table, method = "l1_linf", anchors = window, epsilon = epsilon
  )
  range_exact <- acir::aci_range(
    table, method = "exact", quadrature = "matlab_eps_grid",
    epsilon_grid = epsilon, anchors = window
  )
  full <- acir::aci_online(model, sim, lag = Inf, filter = filt)
  ## The divergence profile at each anchor, every value of it computed by
  ## acir: the relative entropy of the fully informed posterior from the
  ## posterior informed only to each later observation.
  profiles <- lapply(window, function(j) {
    lags <- seq.int(0L, n - j)
    divergence <- vapply(
      lags,
      function(lag) {
        sm <- acir::aci_online(model, sim, lag = lag, filter = filt)
        ## acir's own argument order is (more informed, less informed): the
        ## fully informed posterior first, the lag-limited one second.
        acir::aci_metric_pair(
          mu_p = full$mean[j], R_p = full$cov[j],
          mu_q = sm$mean[j], R_q = sm$cov[j],
          decompose = FALSE
        )
      },
      numeric(1L)
    )
    list(lag = lags * dt, divergence = divergence)
  })
  list(
    x = sim$obs$x, dt = dt, filter = filt, smoother = smooth,
    range = range, range_exact = range_exact, window = window,
    profiles = profiles, epsilon = epsilon
  )
}
