## Shared fixtures for the kalmix test suite.

## A mean-reverting AR(1) / Ornstein-Uhlenbeck series with known persistence.
.sim_ou <- function(n = 1000, ar = 0.85, seed = 1L) {
  withr::with_seed(seed, as.numeric(stats::arima.sim(list(ar = ar), n = n)))
}

## A cointegrated pair: b is a random walk, a = beta * b + u, where the
## spread u is a mean-reverting AR(1) -- the realistic cointegration case.
.sim_pair <- function(n = 400, beta = 2, ar = 0.8, seed = 1L) {
  withr::with_seed(seed, {
    b <- cumsum(stats::rnorm(n))
    u <- as.numeric(stats::arima.sim(list(ar = ar), n = n))
    a <- beta * b + u
    list(a = a, b = b, beta = beta)
  })
}

## A univariate local-level (random-walk-plus-noise) state-space model and a
## realisation from it, with the latent level returned so recovery can be
## checked against the known truth.
.sim_local_level <- function(n = 200,
                             state_sd = 0.3,
                             obs_sd = 1,
                             seed = 1L) {
  withr::with_seed(seed, {
    level <- cumsum(stats::rnorm(n, sd = state_sd))
    y <- level + stats::rnorm(n, sd = obs_sd)
    list(y = y, level = level, state_sd = state_sd, obs_sd = obs_sd)
  })
}

## A two-state Gaussian hidden Markov model and a realisation, returning the
## true state path so the filter and Viterbi can be checked against it.
.sim_hmm <- function(n = 150,
                    means = c(0, 4),
                    sd = 1,
                    stay = 0.95,
                    seed = 1L) {
  k <- length(means)
  trans <- matrix((1 - stay) / (k - 1L), nrow = k, ncol = k)
  diag(trans) <- stay
  withr::with_seed(seed, {
    state <- integer(n)
    state[1L] <- 1L
    for (t in 2:n) {
      state[t] <- sample.int(k, size = 1L, prob = trans[state[t - 1L], ])
    }
    y <- stats::rnorm(n, mean = means[state], sd = sd)
    list(y = y, state = state, means = means, sd = sd, transition = trans)
  })
}
