# changepoint.R -- likelihood-ratio change-point detection in the mean.
#
# detect_changepoint() finds shifts in the mean of a Gaussian sequence. The
# single-change test is the classic likelihood-ratio (equivalently CUSUM)
# statistic: for every candidate split it compares the log-likelihood of one
# common mean against two segment means, and the location maximising the
# improvement is the candidate change-point. Multiple change-points are found
# by binary segmentation -- recursively splitting each segment whose best split
# clears a penalty -- which keeps the search exact-per-segment and transparent.
# A penalised criterion (BIC by default) guards against spurious splits.

#' A change-point detection result
#'
#' An S7 object holding the output of [detect_changepoint()]: the detected
#' change-point locations, the per-segment means, and the penalised criterion
#' the search used.
#'
#' @param changepoints An integer vector of the detected change-point indices
#'   (the last index of each segment but the last); empty if no change is
#'   supported.
#' @param n_segments Integer scalar — the number of segments,
#'   `length(changepoints) + 1`.
#' @param segment_means A numeric vector of the fitted mean of each segment.
#' @param penalty Numeric scalar — the per-split penalty that was applied.
#' @param method Character scalar — the detection method used.
#'
#' @returns An S7 object of class `changepoint_fit`.
#' @family change-point
#' @export
changepoint_fit <- S7::new_class(
  name = "changepoint_fit",
  package = "kalmix",
  properties = list(
    changepoints = S7::class_integer,
    n_segments = S7::class_integer,
    segment_means = S7::class_double,
    penalty = S7::class_double,
    method = S7::class_character
  )
)

#' Detect change-points in the mean of a sequence
#'
#' Finds shifts in the mean of a Gaussian sequence by binary segmentation on
#' the single-change-point likelihood-ratio statistic. At each stage the
#' segment is split at the location whose two-mean fit most improves the
#' Gaussian fit over a single-mean fit; the split is kept only when the
#' improvement clears a penalty, otherwise the segment is declared homogeneous.
#'
#' The single-change statistic at split \eqn{\tau} is the drop in the Gaussian
#' twice-negative-profile-log-likelihood, computed in closed form from the
#' segment sums and sums of squares,
#' \deqn{n \log \hat\sigma^2_{1:n} - \tau \log \hat\sigma^2_{1:\tau} -
#'       (n - \tau)\log \hat\sigma^2_{\tau+1:n},}
#' the variance being estimated within each candidate segment.
#'
#' The penalty matters more than the statistic. Under the null of one
#' homogeneous segment the best-split statistic is far from zero -- on a
#' length-100 series its null mean is near `7` and its 95th percentile near
#' `13` -- so the plain Schwarz penalty \eqn{\log(n)} (about `4.6` at that
#' length) over-segments badly. The default `"mbic"` penalty here is an
#' empirically calibrated \eqn{3\log(n)}, strengthened in the spirit of the
#' modified Bayesian information criterion of Zhang and Siegmund (2007) to
#' account for the change location being estimated (their criterion is a
#' different functional; the name is kept for the shared motivation). On
#' homogeneous Gaussian series it was calibrated to a false-detection rate at
#' or below `6\%` from `n = 100` upward while retaining full power against a
#' one-standard-deviation mean shift. The weaker `"sic"` (\eqn{\log(n)},
#' alias `"bic"`) penalty and any numeric value are offered for callers who
#' want to tune sensitivity deliberately, accepting more false positives.
#'
#' @param x A numeric vector — the sequence to scan.
#' @param penalty Numeric scalar, `"mbic"`, `"sic"`, or `"bic"`. A number is
#'   used directly as the per-split penalty on the statistic; `"mbic"` (the
#'   default) uses \eqn{3\log(n)} and `"sic"` (alias `"bic"`) uses
#'   \eqn{\log(n)}. Larger penalties yield fewer change-points.
#' @param min_segment Integer scalar — the shortest admissible segment length,
#'   so a single outlier is not declared its own segment. Defaults to `5`.
#'
#' @returns A [changepoint_fit].
#' @family change-point
#' @references
#' Scott, A. J. and Knott, M. (1974). A cluster analysis method for grouping
#' means in the analysis of variance. *Biometrics*, 30(3), 507--512. (Binary
#' segmentation.) Killick, R., Fearnhead, P. and Eckley, I. A. (2012). Optimal
#' detection of changepoints with a linear computational cost. *Journal of the
#' American Statistical Association*, 107(500), 1590--1598.
#'
#' Zhang, N. R. and Siegmund, D. O. (2007). A modified Bayes information
#' criterion with applications to the analysis of comparative genomic
#' hybridization data. *Biometrics*, 63(1), 22--32. (The motivation for
#' strengthening the per-split penalty; the `"mbic"` value here is an
#' empirically calibrated \eqn{3\log(n)}, not their functional.)
#' @export
#' @examples
#' set.seed(1)
#' x <- c(rnorm(50, 0), rnorm(50, 4), rnorm(50, 1))
#' fit <- detect_changepoint(x)
#' fit@changepoints
detect_changepoint <- function(x, penalty = "mbic", min_segment = 5L) {
  .check_changepoint_input(x, min_segment)
  n <- length(x)
  pen <- .resolve_changepoint_penalty(penalty, n)

  cps <- sort(.binary_segment(x, 0L, n, pen, as.integer(min_segment)))
  bounds <- c(0L, cps, n)
  seg_means <- vapply(seq_len(length(bounds) - 1L), function(i) {
    mean(x[(bounds[i] + 1L):bounds[i + 1L]])
  }, numeric(1L))

  changepoint_fit(
    changepoints = as.integer(cps),
    n_segments = length(cps) + 1L,
    segment_means = seg_means,
    penalty = pen,
    method = "binary-segmentation"
  )
}

#' Recursive binary segmentation of one interval
#'
#' Finds the best single split of `x[(lo + 1):hi]`, keeps it when the
#' twice-log-likelihood-ratio clears the penalty, and recurses into the two
#' halves. Returns the change-point indices (in the coordinates of the full
#' series) found within this interval.
#'
#' @param x The full numeric series.
#' @param lo Integer scalar — the index before the interval start (zero-based
#'   boundary).
#' @param hi Integer scalar — the last index of the interval.
#' @param pen Numeric scalar — the per-split penalty.
#' @param min_segment Integer scalar — the shortest admissible segment.
#'
#' @returns An integer vector of change-point indices within the interval.
#' @noRd
#' @keywords internal
.binary_segment <- function(x, lo, hi, pen, min_segment) {
  seg_len <- hi - lo
  if (seg_len < 2L * min_segment) {
    return(integer(0L))
  }
  split <- .best_split(x, lo, hi, min_segment)
  if (is.null(split) || split$stat <= pen) {
    return(integer(0L))
  }
  c(
    .binary_segment(x, lo, split$tau, pen, min_segment),
    split$tau,
    .binary_segment(x, split$tau, hi, pen, min_segment)
  )
}

#' Best single split of an interval by the likelihood ratio
#'
#' Scans every admissible split of `x[(lo + 1):hi]` and returns the location
#' and twice-log-likelihood-ratio of the split that most improves the Gaussian
#' fit, computed in closed form from the running segment sums and sums of
#' squares. Returns `NULL` when no admissible split exists.
#'
#' @param x The full numeric series.
#' @param lo,hi The interval boundaries (zero-based `lo`, inclusive `hi`).
#' @param min_segment Integer scalar — the shortest admissible segment.
#'
#' @returns A list with `tau` (the split index) and `stat` (the statistic), or
#'   `NULL`.
#' @noRd
#' @keywords internal
.best_split <- function(x, lo, hi, min_segment) {
  seg <- x[(lo + 1L):hi]
  n <- length(seg)
  total_ss <- .cost_gauss(sum(seg), sum(seg^2), n)

  taus <- seq.int(min_segment, n - min_segment)
  if (length(taus) == 0L) {
    return(NULL)
  }
  csum <- cumsum(seg)
  csum_sq <- cumsum(seg^2)

  stats <- vapply(taus, function(k) {
    left <- .cost_gauss(csum[k], csum_sq[k], k)
    right <- .cost_gauss(csum[n] - csum[k], csum_sq[n] - csum_sq[k], n - k)
    total_ss - left - right
  }, numeric(1L))

  best <- which.max(stats)
  list(tau = lo + taus[best], stat = stats[best])
}

#' Twice the maximised Gaussian log-likelihood cost of a segment
#'
#' Returns \eqn{n \log(\hat\sigma^2)} (the term that drives the
#' likelihood-ratio) computed from a segment's sum and sum of squares, so the
#' segmentation never re-touches the raw data. A guard keeps a degenerate
#' zero-variance segment finite.
#'
#' @param s Numeric scalar — the segment sum.
#' @param ss Numeric scalar — the segment sum of squares.
#' @param n Integer scalar — the segment length.
#'
#' @returns Numeric scalar.
#' @noRd
#' @keywords internal
.cost_gauss <- function(s, ss, n) {
  var_hat <- (ss - s^2 / n) / n
  n * log(max(var_hat, .Machine$double.eps))
}

#' Resolve the change-point penalty argument to a number
#'
#' Maps the named penalties to their length-dependent values -- `"mbic"` to
#' \eqn{3\log(n)} and `"sic"` (alias `"bic"`) to \eqn{\log(n)} -- and passes a
#' numeric value through after checking it is a single non-negative finite
#' number.
#'
#' @param penalty The `penalty` argument to [detect_changepoint()].
#' @param n Integer scalar — the series length.
#'
#' @returns Numeric scalar — the per-split penalty.
#' @noRd
#' @keywords internal
.resolve_changepoint_penalty <- function(penalty, n) {
  if (is.character(penalty)) {
    pen <- switch(penalty,
      mbic = 3 * log(n),
      sic = log(n),
      bic = log(n),
      cli::cli_abort('`penalty` must be "mbic", "sic", "bic", or a number.')
    )
    return(pen)
  }
  pen <- as.numeric(penalty)
  if (length(pen) != 1L || !is.finite(pen) || pen < 0) {
    cli::cli_abort(
      '`penalty` must be "mbic", "sic", "bic", or a non-negative number.'
    )
  }
  pen
}

#' Validate the input to detect_changepoint()
#'
#' @param x The series.
#' @param min_segment The shortest admissible segment.
#'
#' @returns `NULL`, invisibly; called for its error side effect.
#' @noRd
#' @keywords internal
.check_changepoint_input <- function(x, min_segment) {
  if (!is.numeric(x) || !is.null(dim(x))) {
    cli::cli_abort("`x` must be a numeric vector.")
  }
  if (anyNA(x)) {
    cli::cli_abort("`x` must not contain missing values.")
  }
  if (length(min_segment) != 1L || min_segment < 2L) {
    cli::cli_abort("`min_segment` must be a single integer of at least two.")
  }
  if (length(x) < 2L * min_segment) {
    cli::cli_abort(
      "`x` must have at least {2L * as.integer(min_segment)} observations for `min_segment = {min_segment}`."
    )
  }
  invisible(NULL)
}

#' @export
S7::method(print, changepoint_fit) <- function(x, ...) {
  cat(sprintf(
    "<changepoint_fit>: %d change-point%s, %d segment%s (%s)\n",
    length(x@changepoints), if (length(x@changepoints) == 1L) "" else "s",
    x@n_segments, if (x@n_segments == 1L) "" else "s", x@method
  ))
  if (length(x@changepoints) > 0L) {
    cat("  at indices :", paste(x@changepoints, collapse = ", "), "\n")
  }
  cat("  segment means :", paste(round(x@segment_means, 3), collapse = ", "),
      "\n")
  invisible(x)
}
