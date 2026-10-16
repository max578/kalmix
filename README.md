# kalmix <img src="man/figures/logo.png" align="right" height="139" alt="kalmix logo" />

<!-- badges: start -->
[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![R-CMD-check](https://github.com/max578/kalmix/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/max578/kalmix/actions/workflows/R-CMD-check.yaml)
[![Codecov test coverage](https://codecov.io/gh/max578/kalmix/graph/badge.svg)](https://app.codecov.io/gh/max578/kalmix)
[![r-universe](https://max578.r-universe.dev/badges/kalmix)](https://max578.r-universe.dev/kalmix)
[![pkgdown](https://github.com/max578/kalmix/actions/workflows/pkgdown.yaml/badge.svg)](https://max578.github.io/kalmix/)
<!-- badges: end -->

A compact, dependency-light R toolkit for state-space modelling: exact
linear-Gaussian filtering and smoothing, regime switching, hidden Markov
models, change-point detection, and causal read-outs built on the smoother.
Imports only `S7`, `stats` and `cli`.

## What it does

- **State-space core.** `ssm()` declares a linear-Gaussian model (time-varying
  matrices supported, Gaussian or heavy-tailed Student-t observations);
  `kalman_filter()` and `rts_smoother()` run the closed-form forward and
  backward passes; `innovation_diagnostics()` checks model adequacy on the
  standardised innovations.
- **Regimes and discrete states.** `mixture_filter()` is a native Gaussian
  pseudo-Bayesian (GPB1) regime-switching filter; `hmm()`, `hmm_filter()` and
  `hmm_viterbi()` cover Gaussian hidden Markov models; `detect_changepoint()`
  finds mean shifts by penalised likelihood-ratio segmentation.
- **Causal read-outs.** `its_causal()` is a single-case interrupted-time-series
  contrast on the Kalman forecast, and `aci()` implements assimilative causal
  inference (Andreou, Chen and Bollt 2026): the information the future of a
  series carries about the latent state, with an objective decision lead-time
  from `causal_information_rate()` and adequacy-gated abstention.
- **Pairs trading.** The founding application: `spread_series()`, `ou_fit()`
  and `kalmix_filter()` estimate a mean-reverting spread and track its
  volatility regimes, after Elliott, van der Hoek and Malcolm (2005).

## Example

```r
library(kalmix)

## A mean-reverting spread with a calm and a turbulent half.
set.seed(1)
calm <- as.numeric(arima.sim(list(ar = 0.7), n = 150, sd = 0.4))
wild <- as.numeric(arima.sim(list(ar = 0.7), n = 150, sd = 2))
spread <- c(calm, wild)

model <- ou_fit(spread)             # OU / AR(1) mean-reversion fit
half_life(model)                    # trading half-life, in steps

fit <- kalmix_filter(model, spread) # two-regime mixture Kalman filter
tail(round(fit@regime_prob, 2))     # posterior probability of the wild regime
```

The vignettes walk through each verb family: *Getting started: state-space
models with kalmix*, *Pairs trading with kalmix*, *Change-points and hidden
states with kalmix* and *Assimilative causal inference with kalmix*.

## Related packages

- [KFAS](https://CRAN.R-project.org/package=KFAS) filters, smooths,
  forecasts and simulates multivariate state-space models with Gaussian,
  Poisson, binomial, negative binomial and gamma observations.
- [depmixS4](https://CRAN.R-project.org/package=depmixS4) fits hidden
  Markov models to mixed categorical and continuous time series.
- [changepoint](https://CRAN.R-project.org/package=changepoint) finds
  single and multiple changepoints in mean and variance with many
  frequentist and non-parametric methods.
- [CausalImpact](https://CRAN.R-project.org/package=CausalImpact)
  estimates the effect of an intervention on a time series with Bayesian
  structural time-series models.

Each of these covers one of kalmix's areas in more depth. kalmix puts the
linear-Gaussian filter, regime switching, hidden Markov models and
change-points behind one small interface with three dependencies, and adds
assimilative causal inference and the pairs-trading spread model.

## Installation

```r
remotes::install_github("max578/kalmix", build_vignettes = TRUE)
```

## Lineage and credits

The pairs-trading application implements and extends the spread model of
Elliott, van der Hoek and Malcolm (2005), *Pairs trading*, Quantitative
Finance 5(3), 271-276. The package is developed with Johannes van der Hoek.
The assimilative-causal-inference path follows Andreou, Chen and Bollt (2026),
*Assimilative causal inference*, Nature Communications 17, 1854.

## Contributing

Bug reports and suggestions are welcome as
[GitHub issues](https://github.com/max578/kalmix/issues).

## Citation

```r
citation("kalmix")
```

## Licence

MIT. See `LICENSE`.
