# find_zoi.R

#' Detect zones of interest in an intensity trace
#'
#' Identifies the zones of interest (ZOIs) in a smoothed extracted-ion
#' chromatogram: each maximal run of consecutive scans whose intensity exceeds
#' the noise cutoff, together with the index of the apex (highest-intensity
#' point) within the run. ZOIs are the starting points for measuring mass
#' tolerance, peak height, and peak width.
#'
#' A run of length one, an isolated point above the cutoff, is returned as a
#' single-point ZOI. The minimum number of collected masses is applied later,
#' during feature measurement, because following even a single-point ZOI's mass
#' trace into the neighbouring scans can still gather enough values to measure.
#'
#' @param intensity Numeric vector of intensities ordered by scan, typically the
#'   smoothed extracted-ion chromatogram from [smooth_intensity()]. Must be
#'   finite: no `NA`, `NaN`, or infinite values.
#' @param cutoff Single non-negative number giving the noise cutoff, typically
#'   the `cutoff` element returned by [estimate_noise()]. A scan is included when
#'   its intensity is strictly greater than this value.
#'
#' @return A data frame with one row per zone of interest and integer columns
#'   `start`, `end` (first and last scan index of the run) and `apex` (scan
#'   index of the highest point in the run). A trace with no points above the
#'   cutoff returns a zero-row data frame.
#'
#' @examples
#' find_zoi(c(2, 1, 50, 80, 40, 1, 90, 30, 1, 100), cutoff = 10)
#'
#' @export
find_zoi <- function(intensity, cutoff) {
  check_numeric_vector(intensity, "intensity")
  check_nonneg_number(cutoff, "cutoff")
  zoi <- zoi_indices(intensity, cutoff)
  data.frame(
    start = zoi$start,
    end = zoi$end,
    apex = zoi$apex
  )
}

#' Locate zones of interest as bare index vectors
#'
#' The unvalidated core of [find_zoi()]. The measurement loop calls this once
#' per mass bin, tens of thousands of times per file, and only ever reads the
#' index vectors. So it skips the argument checks, the caller having already
#' validated the trace, and the `data.frame()` wrapper that the loop would
#' immediately take apart again.
#'
#' @param intensity Numeric vector of per-scan intensities.
#' @param cutoff Single non-negative noise cutoff.
#'
#' @return A list of three integer vectors, `start`, `end`, and `apex`.
#' @noRd
zoi_indices <- function(intensity, cutoff) {
  above <- which(intensity > cutoff)
  if (length(above) == 0L) {
    return(list(start = integer(0), end = integer(0), apex = integer(0)))
  }

  gaps <- which(diff(above) != 1L)
  run_start <- c(1L, gaps + 1L)
  run_end <- c(gaps, length(above))
  start <- above[run_start]
  end <- above[run_end]
  apex <- start + vapply(
    seq_along(start),
    function(j) which.max(intensity[start[j]:end[j]]) - 1L,
    integer(1)
  )

  list(start = start, end = end, apex = apex)
}
