# smooth_intensity.R

#' Smooth a chromatographic intensity trace
#'
#' Applies a triangular (Bartlett) weighted moving average to an intensity
#' vector. This is the smoothing step used when building extracted-ion
#' chromatograms prior to zone-of-interest detection. The window is symmetric,
#' with `half_window` points on each side of the centre; weights rise linearly
#' towards the centre (`1, 2, ..., half_window + 1, ..., 2, 1`) and are
#' renormalised at the two ends so that each position is averaged only over the
#' points that actually exist. A `half_window` of `0` returns the input unchanged.
#'
#' @param intensity Numeric vector of intensities ordered by scan (retention
#'   time). Must be finite: no `NA`, `NaN`, or infinite values.
#' @param half_window Single non-negative integer giving the number of points
#'   on each side of the smoothing window. `0` (the default) applies no
#'   smoothing.
#'
#' @return A numeric vector the same length as `intensity`.
#'
#' @examples
#' smooth_intensity(c(10, 100, 40, 80, 20), half_window = 1)
#' smooth_intensity(c(10, 100, 40, 80, 20), half_window = 0)
#'
#' @export
smooth_intensity <- function(intensity, half_window = 0L) {
  check_numeric_vector(intensity, "intensity")
  check_nonneg_integer(half_window, "half_window")
  half_window <- as.integer(half_window)

  if (length(intensity) <= 1L || half_window == 0L) {
    return(intensity)
  }

  n <- half_window
  kernel <- c(seq_len(n), n + 1L, rev(seq_len(n)))
  padded <- c(rep(0, n), intensity, rep(0, n))
  mask <- c(rep(0, n), rep(1, length(intensity)), rep(0, n))
  numerator <- stats::filter(
    padded,
    kernel,
    sides = 2
  )
  denominator <- stats::filter(
    mask,
    kernel,
    sides = 2
  )
  keep <- (n + 1L):(n + length(intensity))
  as.numeric(numerator[keep] / denominator[keep])
}
