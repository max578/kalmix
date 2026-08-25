## Independent-oracle fixtures for the ACI conformance suite.
##
## aciR (biometryhub/ACI, public) is a second implementation of Andreou, Chen
## and Bollt (2026) whose numerical core is graded cell-by-cell against the
## method authors' own MATLAB reference (Andreou, github.com/marandmath/
## ACI_code, commit 733c49f, MIT), with the reference fixtures, their
## provenance and the measured agreement shipped in
## `system.file("extdata", "oracle-manifest.yml", package = "aciR")`.
##
## Be exact about where the independence lies. aciR shares a maintainer with
## kalmix, so agreement between the two R packages by itself would establish
## consistency and nothing more. What makes aciR an oracle is the chain its
## numbers carry back to the authors' own code: the values graded below are
## values kalmix did not produce, computed by an implementation that must
## reproduce the reference's MATLAB output to a published tolerance before it
## ships. The sibling `kernR::relative_entropy()` cross-check has no such
## chain -- it is the same textbook closed form in the same estate -- so it
## grades transcription, not correspondence.
##
## aciR works in continuous time on conditional Gaussian nonlinear systems and
## kalmix in discrete time on state-space models, so the two filters are
## different objects that agree only as `dt` tends to zero. What IS directly
## comparable, and what the audit findings KM-01 to KM-03 are about, is the
## reduction of a divergence-versus-lag profile to a causal influence range.
## The fixtures below therefore build the divergence profile with aciR's OWN
## online smoother and metric, and hand it to kalmix's reducer.

## Build the dyad divergence profiles and aciR's reported ranges.
##
## The nonlinear dyad model with intermittent extreme events is the model the
## authors' `dyad_interaction_model.m` implements and the one aciR's `dyad` and
## `cir` oracle fixtures grade.
.aciR_dyad_oracle <- function(n = 601L,
                              dt = 0.001,
                              seed = 11L,
                              window = c(51L, 201L)) {
  model <- aciR::aci_dyad_model()
  sim <- aciR::aci_simulate(model, n = n, dt = dt, seed = seed)
  comp <- aciR::aci_dyad_components(sim$x, model$parameters)
  filt <- aciR::aci_filter(
    sim$x, comp, dt = dt, mu0 = model$y0, R0 = 0.1
  )
  smooth <- aciR::aci_smoother(sim$x, comp, dt = dt, filt = filt)
  range <- aciR::aci_cir(
    sim$x, comp, dt = dt, filt = filt, window = window
  )
  full <- aciR::aci_online_smoother(
    sim$x, comp, dt = dt, filt = filt, lag = Inf
  )
  ## The divergence profile at each anchor, every value of it computed by
  ## aciR: the relative entropy of the fully informed posterior from the
  ## posterior informed only to each later observation.
  profiles <- lapply(window, function(j) {
    lags <- seq.int(0L, n - j)
    divergence <- vapply(
      lags,
      function(lag) {
        sm <- aciR::aci_online_smoother(
          sim$x, comp, dt = dt, filt = filt, lag = lag
        )
        aciR::aci_metric(
          list(mean = sm$mean[j], cov = sm$cov[j]),
          list(mean = full$mean[j], cov = full$cov[j])
        )
      },
      numeric(1L)
    )
    list(lag = lags * dt, divergence = divergence)
  })
  list(
    x = sim$x, dt = dt, filter = filt, smoother = smooth,
    range = range, window = window, profiles = profiles,
    epsilon = range$epsilon
  )
}
