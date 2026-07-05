# kalmix (development version)

## New features

* **Heavy-tailed (Student-t) observation model.** `ssm()` gains an
  `obs_family` argument (`"gaussian"`, the default, or `"student_t"`) and an
  `obs_df` degrees-of-freedom argument. A Student-t family runs a deterministic
  outlier-robust Kalman filter (the variational fixed point of Agamennoni et
  al. 2012 over the Gaussian scale-mixture representation of the t), so an
  outlier is down-weighted rather than allowed to distort the state. The new
  `estimate_obs_df()` selects the degrees of freedom by profile likelihood --
  the grid includes `Inf`, and a likelihood-ratio parsimony margin keeps
  light-tailed data Gaussian. `innovation_diagnostics()` is now family-aware:
  under a Student-t model the scale and family-fit tests read off the t
  probability-integral transform, so the adequacy guard can *certify* a causal
  read-out on a genuinely heavy-tailed series (climate extremes, financial
  returns) where the Gaussian model can only abstain. The Gaussian path is
  unchanged, and the t model collapses to it as the degrees of freedom grow.

* **Model-adequacy grounding for assimilative causal inference.**
  `innovation_diagnostics()` runs the standard state-space residual battery --
  a Ljung-Box whiteness test, a normalised-innovation-squared scale test and a
  Jarque-Bera normality test on the standardised one-step-ahead innovations --
  and reduces it to an adequacy verdict. `aci()` now grounds its read-out
  through it: a model the innovations reject abstains with a `[unverified]`
  token (`grounding`, `grounding_reason` and `adequacy` are new `aci_fit`
  fields), and a verdict is labelled `"grounded"` only when an adequate model
  carries declared, dated mechanism provenance through the new `mechanism`
  argument. The grounding tokens follow the orchestra's provenance vocabulary,
  so passing the diagnostics establishes self-consistency without, on its own,
  claiming external verification.

* **Adaptive online smoother for the causal information rate.**
  `causal_information_rate()` gains an `engine` argument. The default
  `"online"` engine computes the objective rate of Andreou, Chen and Bollt
  (2026, eq. 9) through a fixed-point smoother recursion that grows each
  anchor's smoothed estimate forward in place, so the lead-time costs one
  filter-smoother pass with no re-smoothing and now applies to time-varying
  models. The previous re-smoothing construction is retained as
  `engine = "expanding"` and cross-checks the online engine to numerical
  precision.

## Bug fixes

* **`its_causal()` average-effect intervals were understated.** The augmented
  running-sum covariance recursion dropped the off-diagonal block of its
  process-noise covariance, losing the correlation between the state noise and
  the running sum it also drives. Average-effect intervals now widen slightly
  (the local-level null case understated the summed-counterfactual variance by
  about 4 per cent), so a borderline effect near an interval endpoint can
  change verdict. Point estimates, pointwise counterfactual intervals and the
  cumulative-effect path are unchanged. The recursion is now pinned to a
  direct Monte-Carlo simulation oracle in the test suite, and the null
  coverage test tightened from a one-sided floor to a two-sided band around
  the nominal 95 per cent.

* The NOAA Oceanic Nino Index external-oracle test now skips cleanly when the
  `curl` package is absent, so `R CMD check` is clean on hosts without it
  (the offline guard previously errored rather than skipping).

# kalmix 0.2.0

This development cycle adds the time-resolved assimilative-causal-inference
(ACI) path on top of the Phase 1 state-space core.

## New features

* **Assimilative causal inference.** `aci()` reads a time-resolved causal
  information series off the Kalman filter and Rauch-Tung-Striebel smoother of
  an `ssm`: the per-step relative entropy of the smoothing distribution of the
  latent state from its filtering distribution (Andreou, Chen and Bollt 2026,
  eq. 7), the information the future of the series adds to the present estimate
  of the state. It returns an `aci_fit` carrying the causal information series,
  its time-average, and -- by default -- the objective causal information rate
  and the implied decision lead-time.
* **Causal information rate and lead-time.** `causal_information_rate()`
  computes the threshold-free objective causal information rate of Andreou,
  Chen and Bollt (2026, eqs. 8--9) by integrating an expanding-future-window
  divergence profile, returning a decision lead-time in the time units of the
  series.
* **Native Gaussian relative entropy.** The smoother-vs-filter divergence is
  computed natively in closed form, so the ACI path is self-contained. When the
  orchestra's `kernR` package is installed its `relative_entropy()` is used as
  an independent oracle in the test suite; the comparison skips when `kernR` is
  absent, so the package builds and checks cleanly either way (`kernR` is
  `Suggests`-only).
* **Validation.** A synthetic recharge-oscillator ENSO model with a known
  hidden-state data-generating process is the primary self-consistency oracle
  (the causal information is zero under the null and positive when coupled, the
  hidden state is recovered from the observed effect alone, and the lead-time
  grows with observation noise); the real NOAA Oceanic Nino Index is the
  external oracle, run by an online-gated test against a recorded fixture. A new
  *Assimilative causal inference with kalmix* vignette walks through both.


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
