## Self-contained OSSE fixtures for the Assimilative Causal Inference suite.
##
## A synthetic recharge-oscillator ENSO model with a KNOWN hidden-state
## data-generating process is the primary (self-consistency) oracle: because we
## simulate the truth, we can assert the ACI read-out behaves as it must. The
## hidden state z = (T, h) is the SST anomaly T and the thermocline-depth
## anomaly h evolving under the noise-sustained recharge oscillator of Jin
## (1997) -- complex eigenvalues give an irregular, roughly four-year ENSO-like
## oscillation. The OBSERVED series is a scalar that warm T drives through a
## coupling kappa (a crop-water-stress index in the grain anchor); kappa = 0 is
## the null under which the observation carries no information about the hidden
## cause. This mirrors the orchestra's conceptual ACI ENSO validation but is
## kalmix's own discrete-time, ssm-shaped build.

## Simulate the recharge oscillator and the coupled observation. The continuous
## drift dz = M z dt is discretised by the Euler-Maruyama transition
## A = I + M dt, so the kalmix ssm built from `amat` matches the DGP exactly.
.sim_enso_osse <- function(n = 600L,
                           kappa = 0.6,
                           dt = 1 / 12,
                           obs_sd = 0.5,
                           seed = 1L) {
  growth_rate <- 0
  coupling_gh <- 1.6
  coupling_ah <- 1.6
  damping_h <- 0.2
  state_sd_t <- 0.35
  state_sd_h <- 0.25
  drift <- matrix(
    c(growth_rate, coupling_gh, -coupling_ah, -damping_h),
    nrow = 2L, byrow = TRUE
  )
  amat <- diag(2L) + drift * dt

  withr::with_seed(seed, {
    z <- matrix(0, nrow = n, ncol = 2L)
    x <- numeric(n)
    for (t in 2:n) {
      noise <- c(
        stats::rnorm(1L, sd = state_sd_t * sqrt(dt)),
        stats::rnorm(1L, sd = state_sd_h * sqrt(dt))
      )
      z[t, ] <- as.numeric(amat %*% z[t - 1L, ]) + noise
      x[t] <- kappa * z[t, 1L] + stats::rnorm(1L, sd = obs_sd)
    }
    list(
      x = x, z = z, dt = dt, kappa = kappa, amat = amat,
      state_sd_t = state_sd_t, state_sd_h = state_sd_h, obs_sd = obs_sd
    )
  })
}

## Build the kalmix ssm that matches the known OSSE data-generating process.
## The observation operator maps T (the first hidden component) to the observed
## series through the coupling; a tiny floor is used at the null so the
## observation matrix is non-degenerate while still carrying no information.
.enso_osse_model <- function(sim, kappa = sim$kappa) {
  kappa <- if (kappa == 0) 1e-6 else kappa
  ssm(
    transition = sim$amat,
    observation = matrix(c(kappa, 0), nrow = 1L),
    state_cov = diag(c(
      sim$state_sd_t^2 * sim$dt,
      sim$state_sd_h^2 * sim$dt
    )),
    obs_cov = sim$obs_sd^2,
    init_state = c(0, 0),
    init_cov = diag(c(1, 1))
  )
}
