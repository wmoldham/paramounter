# flag_isolate_zoi.R

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
#' later fail the mass-count or ppm filters in [zoi_features()]: a nearby peak
#' interferes with matching whether or not it yields a usable measurement. Apply
#' this function to all apex retention times from [find_zoi()] and combine the
#' result with the measurement filters afterwards.
#'
#' @param apex_rt Numeric vector of apex retention times for the zones of
#'   interest in one mass bin, in non-decreasing order (as returned by
#'   [find_zoi()], whose zones are ordered by scan). Must be finite.
#' @param min_gap Single non-negative number giving the minimum separation, in
#'   the same units as `apex_rt`, for a zone to count as isolated. Defaults to
#'   `300` (five minutes, for retention times in seconds).
#' @param legacy_isolation Single `TRUE` or `FALSE`. When `TRUE` (the default)
#'   the original Paramounter behaviour is reproduced: every zone except the
#'   last is tested only against the *following* zone, and the last zone only
#'   against the preceding one, so a zone immediately after a close neighbour is
#'   still treated as clean. When `FALSE` a zone must be separated from
#'   neighbours on both sides, which is the criterion described in the paper.
#'   This toggle exists only to allow comparison against the published values
#'   and is expected to be removed once reproduction is confirmed.
#'
#' @return A logical vector the same length as `apex_rt`, `TRUE` where the zone
#'   is isolated.
#'
#' @examples
#' # the middle zone sits 10 s after its neighbour but 690 s before the next
#' flag_isolated_zoi(c(100, 110, 800))
#' flag_isolated_zoi(c(100, 110, 800), legacy_isolation = FALSE)
#'
#' @export
flag_isolated_zoi <- function(
    apex_rt,
    min_gap = 300,
    legacy_isolation = TRUE
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
