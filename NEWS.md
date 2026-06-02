# kalmix 0.0.0.9000

## New features

* First scaffold: a standalone spread-modelling core plus the dormant
  mixture-Kalman adapter that commissions `proxymix::gmm_filter()`.
* **New: `spread_series()`** -- builds a cointegration spread from two price
  series via an ordinary-least-squares hedge ratio.
* **New: `ou_fit()`** -- fits the Ornstein-Uhlenbeck (mean-reverting)
  dynamics of a spread by autoregression, returning a `spread_model` with
  the mean-reversion speed, long-run mean, volatility, and half-life.
* **New: `spread_model` S7 class** with a `print()` method and `half_life()`.
* **New (dormant): `kalmix_filter()`** -- the interface for filtering a
  regime-switching spread state through a finite Gaussian-mixture Kalman
  recursion. Capability-probed against `proxymix::gmm_filter()`; errors
  cleanly with a commissioning notice until that primitive ships
  (constellation Invariant 6, dormant-but-correct).
