# flag_isolated_zoi.R

#' Flag zones of interest that have no close neighbour
#'
#' Identifies the "clean" zones of interest within a mass bin: those far enough
#' in retention time from adjacent zones that they can be matched unambiguously
#' across samples. Only clean zones are used to estimate instrument mass and
#' retention-time shift, since a zone with a close neighbour risks being matched
#' to the wrong peak in another file.
#'
#' A zone qualifies when the gap to its neighbour is strictly greater than
#' `min_gap`. A bin containing a single zone always qualifies.
#'
#' The neighbour set is every zone detected in the bin, including zones that
#' later fail the mass-count or ppm filters in `zoi_features()`: a nearby peak
#' interferes with matching whether or not it yields a usable measurement. Pass
#' every zone detected in the bin, then apply the measurement filters afterwards.
#'
#' @param apex_rt Numeric vector of apex retention times for the zones in one
#'   mass bin, in non-decreasing order. Must be finite. `find_zoi()` returns scan
#'   indices rather than times, so index your retention times by its `apex`
#'   column; its zones are already ordered by scan.
#' @param min_gap Single non-negative number giving the minimum separation, in
#'   the same units as `apex_rt`, for a zone to count as isolated. Defaults to
#'   `300` (five minutes, for retention times in seconds).
#' @param legacy_isolation Single `TRUE` or `FALSE`. When `FALSE` (the default)
#'   a zone must be separated from its neighbours on both sides, the criterion
#'   the paper describes. When `TRUE` it reproduces the original, which tests
#'   every zone but the last against only the *following* zone, and the last
#'   against only the preceding one. [paramounter()] sets this from [pm_config]'s
#'   `legacy`.
#'
#' @return A logical vector the same length as `apex_rt`, `TRUE` where the zone
#'   is isolated.
#'
#' @noRd
flag_isolated_zoi <- function(
    apex_rt,
    min_gap = 300,
    legacy_isolation = FALSE
) {
  check_numeric_vector(apex_rt, "apex_rt")
  check_nonneg_number(min_gap, "min_gap")
  check_flag(legacy_isolation, "legacy_isolation")

  n <- length(apex_rt)
  if (n == 0L) {
    return(logical(0))
  }
  if (n > 1L && any(diff(apex_rt) < 0)) {
    stop(
      "`apex_rt` must be in non-decreasing order.",
      call. = FALSE
    )
  }

  gaps <- diff(apex_rt) > min_gap
  forward <- c(gaps, TRUE)
  backward <- c(TRUE, gaps)

  if (legacy_isolation) {
    out <- forward
    out[n] <- backward[n]
    out
  } else {
    forward & backward
  }
}
