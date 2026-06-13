# kalmix 0.1.0

This release turns kalmix from a thin pairs-trading scaffold into a general
state-space toolkit, and relicences it as open-source MIT. Mean-reverting
spread estimation is retained as one documented application of the general
machinery.

## Breaking changes

* The package is **relicensed from proprietary to MIT** (`LICENSE`,
  `LICENSE.md`, and `DESCRIPTION` `License: MIT + file LICENSE`). kalmix is now
  a public, open-source member of its package family.
* `kalmix_filter()` is no longer a dormant adapter that errors until
  `proxymix::gmm_filter()` ships. It now runs a real native regime filter and
  works standalone, with a new signature:

  ```r
  # before (0.0.0.9000): errored unless proxymix shipped a primitive
  kalmix_filter(model, y, n_regimes = 2)

  # after (0.1.0.9000): runs the native engine; proxymix is optional
  kalmix_filter(model, y, n_regimes = 2, vol_multipliers = NULL,
                engine = "native")
  ```

  The `...` argument is replaced by explicit `vol_multipliers` and `engine`
  arguments.

## New features

* **General linear-Gaussian state space.** A new `ssm()` constructor builds a
  vector-state, optionally time-varying, linear-Gaussian state-space model as
  an S7 object. `kalman_filter()` runs the closed-form forward recursion of
  Kalman (1960) and `rts_smoother()` the backward Rauch-Tung-Striebel (1965)
  smoother, returning `kalman_fit` and `rts_fit` objects with the filtered,
  predicted and smoothed states, innovations, and model log-likelihood.
* **Regime-switching mixture filter.** `mixture_filter()` filters a series that
  switches between a finite set of linear-Gaussian regimes through a native
  Gaussian pseudo-Bayesian (GPB1) recursion, returning a `regime_fit` with the
  per-step posterior over the active regime. It is self-contained; an optional
  `engine = "proxymix"` path is honoured only if that package is installed.
* **Hidden Markov models.** `hmm()` builds a discrete-state Gaussian hidden
  Markov model; `hmm_filter()` runs the forward-backward recursion for filtered
  and smoothed state posteriors and the log-likelihood, and `hmm_viterbi()`
  returns the jointly most likely state path.
* **Change-point detection.** `detect_changepoint()` finds shifts in the mean
  of a Gaussian sequence by binary segmentation on the likelihood-ratio
  statistic, with an empirically calibrated default penalty (`3 log n`) that
  holds the false-detection rate at or below six per cent on homogeneous series
  while retaining full power against a one-standard-deviation shift.
* **Single-case interrupted-time-series causal contrast.** `its_causal()`
  estimates the causal effect of an intervention on one series by forecasting a
  state-space counterfactual from the pre-intervention segment and contrasting
  the observed post-intervention path against it. The counterfactual's noise
  scale is fitted by maximum likelihood so the effect interval is
  well-calibrated -- a genuinely null effect reports an interval that straddles
  zero.
* **Pairs-trading application retained.** `spread_series()`, `ou_fit()`, the
  `spread_model` class and `half_life()` are unchanged; `kalmix_filter()` now
  expands a fitted spread into volatility regimes for the general
  `mixture_filter()` engine.

## Minor improvements and fixes

* The package title and description are reframed from pairs-trading-only to a
  general state-space toolkit.
* `DESCRIPTION` gains `URL:` and `BugReports:` fields and a `Collate:` order so
  the S7 classes load before the functions that reference them.
