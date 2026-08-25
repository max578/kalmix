# kalmix (development version)

## Breaking changes

* **The objective causal influence range is normalised by the peak of the
  divergence profile, not by its value at lag zero.** This is the
  normalisation of Andreou, Chen and Bollt (2026), and the two agree only
  when the profile decreases with lag, which on a real record it often does
  not. Reported lead-times change, in either direction; on the package's own
  dyad conformance fixture the old normaliser was 58 per cent high at one
  anchor.
* **`causal_information_rate()` evaluates every lag in the window by
  default** (`n_lag = NULL`), where it previously used a twelve-point grid.
  A twelve-point trapezoid over a sharply decaying profile over-estimates
  the integral; on the package's ENSO fixture the old default was about
  9 per cent high against the fully resolved quadrature. Pass an explicit
  `n_lag` for the former behaviour.
* **An anchor whose future window runs past the end of the record is
  right-censored rather than extrapolated.** The divergence at the last
  resolvable lag was previously carried forward as a constant over every
  remaining lag, which inflated the integral without warning; at step 580 of
  a 600-step series with a 120-step window that produced a range of 32.7
  steps from twenty steps of data. The integral now stops at the last
  resolved lag and the result is flagged.

## New features

* `causal_information_rate()` gains `functional`, choosing between the
  efficient integral form of the range (`"objective"`, the default, which
  the source paper gives as an underestimate) and the threshold-averaged
  form (`"objective_exact"`), and `epsilon` for the tolerance grid the
  latter averages over.
* `causal_information_rate()` gains `margin` and `tol`, and returns its
  scalar with `censored`, `censored_fraction`, `monotone` and `n_lag`
  attributes. A thinned quadrature grid is compared with the same reduction
  on every other point and warns when the two disagree by more than `tol`.
* `aci()` gains `n_lag`, `engine` and `functional`, so the flagship verb can
  refine and re-specify the lead-time it reports; previously neither the
  grid nor the engine could be reached through it.
* `aci_fit` gains `censored` and `monotone` properties, `print()` marks a
  censored lead-time as a lower bound, and a censored anchor degrades the
  grounding token to `[unverified]` with reason `"censored_horizon"`.
* `causal_information_rate()` returns a `converged` attribute alongside
  `censored`: `TRUE`/`FALSE` when the convergence check ran, `NA` when it
  did not apply (a resolved grid, or `tol = Inf`). `aci()` reads it: a
  lead-time whose quadrature has not converged now also degrades the
  grounding token to `[unverified]`, with reason `"not_converged"`, and the
  returned `aci_fit` carries the `orchestra_refusal` / `kalmix_abstention`
  classes (the same convention `ou_fit()`'s refusal now uses, below) so a
  cross-member caller can recognise the lead-time as unreliable without
  parsing the warning. Previously the convergence check only warned and left
  the grounding token untouched.

## Bug fixes

* The documented functional now matches the computed one. The roxygen
  previously asserted an integration-by-parts identity between the
  threshold-averaged and integral forms that holds only for a decreasing
  profile, while defining the subjective range as the *first* lag below the
  tolerance where the source paper takes the *last* lag above it.
* `ou_fit()`'s refusal on a non-mean-reverting series now carries the
  `orchestra_refusal` / `kalmix_refusal` condition classes, so it can be
  caught and recognised as a structured refusal by name rather than by
  parsing its message; the refusal itself -- what is refused and when -- is
  unchanged.
* `aci()` now also carries the `orchestra_refusal` / `kalmix_abstention`
  classes when the model itself fails the innovation-adequacy diagnostics
  (`grounding_reason = "model_inadequate"`) or when there is too little data
  to judge adequacy at all (`"insufficient_data"`) -- previously only a
  censored or unconverged lead-time was classed, so an inadequate model
  whose lead-time happened to be reliable (or was not requested,
  `lead_time = FALSE`) returned an un-classed `aci_fit` carrying only the
  `"[unverified]"` grounding string. The ordinary, expected
  `"mechanism_unverified"` case -- an adequate model with no declared
  provenance -- is unaffected and still carries no abstention class; only
  the genuine adequacy failures now do. What is refused is unchanged, only
  how it is signalled.

## Testing

* New conformance suite (`test-cir-conformance.R`) grading the ACI core
  against `aciR` (`biometryhub/ACI`), a second implementation of the same
  method whose numerical core is graded cell-by-cell against the method
  authors' MATLAB reference; that chain back to the authors' own code, not
  authorship, is what makes it an oracle here. The
  causal metric, the objective range, the subjective range at 129
  tolerances and the exact range are each graded on the nonlinear dyad
  model. `aciR` joins `Suggests`; the tests skip when it is absent.
* New quadrature suite (`test-cir-quadrature.R`) grading the lag grid, the
  convergence warning, the censoring contract and the grounding degradation,
  and now also the convergence-driven grounding degradation and the
  `orchestra_refusal`/`kalmix_abstention` classes on an unconverged
  `aci_fit`. `test-spread.R` grades the `orchestra_refusal`/`kalmix_refusal`
  classes on `ou_fit()`'s refusal.
* The check workflow installs `aciR` from GitHub rather than dropping it from
  `Suggests`, so the conformance grade runs in the clean room on all five
  platforms; `kernR` is still dropped, being private. `testthat`'s floor rises
  to 3.1.5 for `expect_no_warning()`.

## Documentation

* Every figure-producing chunk across the five vignettes now carries a
  `fig.cap` and at least one prose sentence interpreting the figure
  (previously none of the ten figures were captioned).
* *Assimilative causal inference with kalmix* states which range functional
  the reported lead-time is, that the two engines cross-check the smoother
  recursion and not the quadrature, and how a censored anchor is reported.
* *Change-points and hidden states with kalmix* gains the governing
  likelihood-ratio statistic for `detect_changepoint()`'s binary
  segmentation, cited to Scott and Knott (1974); the reference list now
  attributes binary segmentation to Scott and Knott rather than to
  Killick, Fearnhead and Eckley (2012), whose linear-cost PELT algorithm
  is noted as a "see also" alternative rather than the implemented method.
* *Change-points and hidden states with kalmix* gains the governing
  likelihood-ratio statistic for `detect_changepoint()`'s binary
  segmentation, cited to Scott and Knott (1974); the reference list now
  attributes binary segmentation to Scott and Knott rather than to
  Killick, Fearnhead and Eckley (2012), whose linear-cost PELT algorithm
  is noted as a "see also" alternative rather than the implemented method.

# kalmix 0.5.0 (2026-07-05)

## New features

* **Intercepts (control inputs).** `ssm()` gains `state_intercept` (`c_t`)
  and `obs_intercept` (`d_t`) in the standard
  Durbin-Koopman form, static or time-varying, defaulting to zero. They
  thread through both filter families, the regime machinery and
  `ssm_fit()`. A stationary state with a non-zero long-run mean is now a
  native one-liner (`state_intercept = theta * (1 - phi)`), and the
  intercept-form regime filter reproduces `kalmix_filter()`'s shipped
  deviation-form path exactly, so the two mechanisms validate each other.

* **Regime (Kim) smoother.** The new `mixture_smoother()` runs the Kim
  (1994) backward pass over the GPB1 forward recursion, returning smoothed
  regime probabilities and collapsed state moments. With identical regimes
  it reproduces `rts_smoother()` exactly; on a known switching path its
  regime attribution is at least as accurate as the filter's.

* **Assimilative causal inference on regime models.** `aci()` and
  `causal_information_rate()` accept a list of `ssm` regimes (with the
  chain's `transition` and `init_prob`): the causal information compares
  the Kim smoother's collapsed posteriors with the GPB1 filter's, the
  verdict grounds through the mixture-PIT diagnostics, and the lead-time
  uses the expanding-window construction. Both densities are the model
  class's own Gaussian collapses, an approximation the documentation
  states.

* **New vignette.** *The Elliott-van der Hoek-Malcolm spread model with
  kalmix* reproduces the founding 2005 model with the general machinery
  and demonstrates why the paper filters: least squares on the noisy
  observed spread attenuates the autoregression and misprices the
  half-life several-fold, while `ssm_fit()` recovers the true dynamics.

# kalmix 0.4.0 (2026-07-05)

## New features

* **Missing observations.** A fully `NA` observation (row) is now a missing
  observation everywhere: `kalman_filter()` runs a prediction-only step
  under both observation families, `mixture_filter()` carries the regime
  chain prior through the gap, the HMM verbs use a flat emission, and
  `innovation_diagnostics()` tests the complete steps. Gappy series
  (market holidays, sensor dropouts, index gaps) fit, filter, smooth and
  certify directly. Partially missing multivariate rows are refused with a
  clear error. The gappy filter and smoother are pinned to a direct
  joint-Gaussian conditioning oracle at every step, and the gappy
  log-likelihood to the joint density of the observed entries.

* **Maximum-likelihood estimation.** The new `ssm_fit()` estimates the
  process covariance, the observation covariance and (optionally) the
  transition matrix of an `ssm` by maximising the prediction-error
  log-likelihood over log-Cholesky factors, from a deterministic
  multi-start ladder, with delta-method standard errors. The returned
  `ssm_mle` carries a fitted `ssm` ready for every verb. Fitted variances
  agree with `stats::StructTS` on nested models and recover known
  parameters in simulation within their standard errors.

* **Regime certification.** `mixture_filter()` and `kalmix_filter()` now
  expose the collapsed one-step innovations and predictive covariances on
  `regime_fit`, plus (univariate) the exact probability integral transform
  of each observation under its mixture predictive.
  `innovation_diagnostics()` accepts a `regime_fit` and tests the PIT
  (reported as `obs_family = "gaussian_mixture"`), so a regime-switching
  model can be certified directly; the collapsed-Gaussian reference would
  wrongly fail a well-specified mixture, whose innovations are
  heavy-tailed relative to one Gaussian by construction. The full `aci()`
  read-out on a regime fit awaits a regime smoother and says so.

# kalmix 0.3.0 (2026-07-05)

## Breaking changes

* **The optional `proxymix` engine is removed.** `mixture_filter()` and
  `kalmix_filter()` no longer take an `engine` argument, and `proxymix` has
  left `Suggests`. The delegated path assumed a regime-switching signature
  that `proxymix::gmm_filter()` does not provide (it is a mixture-prior
  Gaussian-sum filter), so requesting that engine could never run the
  delegated call. The package's native Gaussian pseudo-Bayesian (GPB1)
  recursion, always the default, is now the only engine. Callers who passed
  `engine = "native"` can simply drop the argument; results are unchanged.

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
  argument. Passing the diagnostics establishes self-consistency without, on
  its own, claiming external verification.

* **Adaptive online smoother for the causal information rate.**
  `causal_information_rate()` gains an `engine` argument. The default
  `"online"` engine computes the objective rate of Andreou, Chen and Bollt
  (2026, eq. 9) through a fixed-point smoother recursion that grows each
  anchor's smoothed estimate forward in place, so the lead-time costs one
  filter-smoother pass with no re-smoothing and now applies to time-varying
  models. The previous re-smoothing construction is retained as
  `engine = "expanding"` and cross-checks the online engine to numerical
  precision.

## Package infrastructure

* A README entry surface, a three-OS by three-R
  continuous-integration check matrix, an explicit pre-1.0 API-stability
  policy (`API_STABILITY.md` in the repository), and a `CITATION` that reads
  the version from the package metadata instead of pinning it by hand.

## Bug fixes

* **The causal-influence-range profile now matches its published
  definition.** The expanding-window divergence behind
  `causal_information_rate()` and the `aci()` lead-time measured the lagged
  estimate from the complete smoother; the published definition integrates
  over the complete smoother (the reverse relative-entropy direction), and
  its normaliser is then exactly the per-step causal information at the
  anchor. Both engines were corrected together, so lead-times shift slightly
  in general; the two engines still agree to numerical tolerance and the
  recorded ENSO sanity bounds are unchanged. A dated equation-grounding
  record now accompanies the package sources.

* **`kalmix_filter()` now tracks a spread with a non-zero long-run mean.** The
  per-regime state-space models carried no intercept term, so the filtered
  state was contracted toward zero rather than toward the fitted long-run
  mean: on a spread centred at 5 the filtered mean sat roughly 6 per cent low.
  The regime models are now the Ornstein-Uhlenbeck deviation form -- the
  filter runs on `y - theta` and the returned state estimates are shifted
  back to the observed scale, for both the mixture path and the
  single-regime Kalman path. Regime probabilities, innovations and the
  log-likelihood are unaffected by the reformulation, and for a spread whose
  fitted long-run mean is near zero the change is negligible.

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
  computed natively in closed form, so the ACI path is self-contained. When
  the `kernR` package is installed its `relative_entropy()` is used as
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
