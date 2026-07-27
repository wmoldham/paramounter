# collect_zoi_masses.R

#' Follow a mass trace outward from a zone-of-interest apex
#'
#' Starting at the apex of a zone of interest, walks left and then right across
#' adjacent scans, collecting at each scan the m/z value nearest the reference
#' m/z. The walk in a given direction stops at the end of the trace, at a scan
#' with no signal in the mass bin, or at a scan whose mean intensity falls below
#' the noise cutoff.
#'
#' It also stops when the m/z it collects is a Tukey-fence outlier against the
#' values gathered so far, which marks an adjacent coeluting peak. That test
#' applies only once more than `min_points` values have been collected and the
#' smoothed intensity has begun rising again.
#'
#' The collected m/z values are used downstream to measure mass tolerance, and
#' the scan boundaries to measure peak width.
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
  check_list(mz_list, "mz_list")
  check_list(int_list, "int_list")
  check_numeric_vector(eic, "eic")
  n <- length(eic)
  if (length(mz_list) != n || length(int_list) != n) {
    stop(
      "`mz_list` and `int_list` must have the same length as `eic`.",
      call. = FALSE
    )
  }
  check_matching_lengths(mz_list, int_list)
  check_nonneg_number(cutoff, "cutoff")
  check_index(apex, "apex", n, "eic")
  check_finite_number(reference_mz, "reference_mz")
  check_nonneg_integer(min_points, "min_points")
  check_nonneg_number(fence, "fence")
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
