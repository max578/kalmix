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
