# collect_zoi_masses.R

#' Follow a mass trace outward from a zone-of-interest apex
#'
#' Starting at the apex of a zone of interest, walks left and then right across
#' adjacent scans, collecting at each scan the m/z value nearest the reference
#' m/z. The walk in a given direction stops when it reaches the end of the
#' trace, a scan with no signal in the mass bin, a scan whose mean intensity
#' falls below the noise cutoff, or (once more than `min_points` values have been
#' collected and the smoothed intensity has begun rising again, indicating an
#' adjacent coeluting peak) a scan whose m/z is a Tukey-fence outlier relative to
#' the values gathered so far. The collected m/z values are used downstream to
#' measure mass tolerance, and the scan boundaries to measure peak width.
#'
#' @param mz_list List with one element per scan; each element is the numeric
#'   vector of m/z values falling in the mass bin at that scan (possibly empty).
#' @param int_list List parallel to `mz_list` giving the corresponding
#'   intensities; element `i` must have the same length as `mz_list[[i]]`.
#' @param eic Numeric vector of per-scan intensities for the bin (the smoothed
#'   extracted-ion chromatogram). Must be finite and the same length as
#'   `mz_list`. Used to detect the rise into an adjacent peak.
#' @param cutoff Single non-negative noise cutoff. A scan is included while its
#'   mean bin intensity is at or above this value.
#' @param apex Single integer scan index of the ZOI apex (the starting point),
#'   typically the `apex` column from [find_zoi()].
#' @param reference_mz Single finite reference m/z, typically the m/z of the most
#'   intense point at the apex. The nearest value to this is taken at each scan.
#' @param min_points Single non-negative integer. The outlier test is applied
#'   only once more than this many values have been collected (default `5`).
#' @param fence Single non-negative multiplier for the inter-quartile range in
#'   the Tukey outlier test (default `1.5`).
#'
#' @return A named list with elements `masses` (the collected m/z values,
#'   beginning with `reference_mz`), `left_edge` and `right_edge` (the leftmost
#'   and rightmost included scan indices; both equal `apex` when no neighbouring
#'   scans are included).
#'
#' @examples
#' mz_list <- list(
#'   numeric(0), 220.0102, 220.0098, 220.0100, 220.0101, 220.0103, numeric(0)
#' )
#' int_list <- list(numeric(0), 1000, 3000, 5000, 3000, 1000, numeric(0))
#' eic <- c(0, 1000, 3000, 5000, 3000, 1000, 0)
#' collect_zoi_masses(mz_list, int_list, eic, cutoff = 500, apex = 4, reference_mz = 220.0100)
#'
#' @export
collect_zoi_masses <- function(
    mz_list,
    int_list,
    eic,
    cutoff,
    apex,
    reference_mz,
    min_points = 5L,
    fence = 1.5
) {
  n <- length(eic)
  if (!is.numeric(eic)) {
    stop(
      "`eic` must be a numeric vector.",
      call. = FALSE
    )
  }
  if (any(!is.finite(eic))) {
    stop(
      "`eic` must not contain NA, NaN, or infinite values.",
      call. = FALSE
    )
  }
  if (!is.list(mz_list)) {
    stop(
      "`mz_list` must be a list.",
      call. = FALSE
    )
  }
  if (!is.list(int_list)) {
    stop(
      "`int_list` must be a list.",
      call. = FALSE
    )
  }
  if (length(mz_list) != n || length(int_list) != n) {
    stop(
      "`mz_list` and `int_list` must have the same length as `eic`.",
      call. = FALSE
    )
  }
  if (!all(lengths(mz_list) == lengths(int_list))) {
    stop(
      "Each element of `mz_list` and `int_list` must have matching length.",
      call. = FALSE
    )
  }
  if (!is.numeric(cutoff) ||
      length(cutoff) != 1L ||
      !is.finite(cutoff) ||
      cutoff < 0) {
    stop(
      "`cutoff` must be a single non-negative number.",
      call. = FALSE
    )
  }
  if (!is.numeric(apex) ||
      length(apex) != 1L ||
      !is.finite(apex) ||
      apex != as.integer(apex) ||
      apex < 1 ||
      apex > n) {
    stop(
      "`apex` must be a single integer index within `eic`.",
      call. = FALSE
    )
  }
  if (!is.numeric(reference_mz) ||
      length(reference_mz) != 1L ||
      !is.finite(reference_mz)) {
    stop(
      "`reference_mz` must be a single finite number.",
      call. = FALSE
    )
  }
  if (!is.numeric(min_points) ||
      length(min_points) != 1L ||
      !is.finite(min_points) ||
      min_points < 0 ||
      min_points != as.integer(min_points)) {
    stop(
      "`min_points` must be a single non-negative integer.",
      call. = FALSE
    )
  }
  if (!is.numeric(fence) ||
      length(fence) != 1L ||
      !is.finite(fence) ||
      fence < 0) {
    stop(
      "`fence` must be a single non-negative number.",
      call. = FALSE
    )
  }
  apex <- as.integer(apex)

  walk <- function(masses, dir) {
    ind <- apex + dir
    repeat {
      if (ind < 1L || ind > n) break
      if (length(mz_list[[ind]]) == 0L) break
      if (mean(int_list[[ind]]) < cutoff) break
      candidates <- mz_list[[ind]]
      nearest <- if (length(candidates) == 1L) {
        candidates
      } else {
        candidates[which.min(abs(candidates - reference_mz))]
      }
      masses <- c(masses, nearest)
      if (eic[ind] > eic[ind - dir] && length(masses) > min_points) {
        quart <- stats::quantile(
          masses,
          c(0.25, 0.75),
          names = FALSE,
          type = 7
        )
        spread <- quart[2] - quart[1]
        latest <- masses[length(masses)]
        if (latest < quart[1] - fence * spread ||
            latest > quart[2] + fence * spread) {
          break
        }
      }
      ind <- ind + dir
    }
    list(masses = masses, edge = ind - dir)
  }

  left <- walk(reference_mz, -1L)
  right <- walk(left$masses, 1L)
  list(
    masses = right$masses,
    left_edge = left$edge,
    right_edge = right$edge
  )
}
