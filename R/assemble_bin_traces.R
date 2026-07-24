# assemble_bin_traces.R

#' Reconstruct one bin's per-scan traces
#'
#' Expands a single bin from [bin_peaks()] into the full-length per-scan
#' structures the measurement functions consume: a list of the m/z values at
#' each scan, the matching intensities, and the raw extracted-ion chromatogram
#' (the mean intensity in the bin at each scan, zero where the bin is empty).
#' Scans with no peak in the bin give an empty vector and a zero EIC value.
#'
#' This is done one bin at a time rather than for all bins at once, so the
#' full-length structures never accumulate across the whole m/z axis.
#'
#' @param bin A single bin element from the `bins` list returned by
#'   [bin_peaks()]: a list with parallel numeric vectors `scan`, `mz`, and
#'   `intensity`, where `scan` holds integer scan indices within `[1, n_scans]`
#'   and `mz`, `intensity` are finite.
#' @param n_scans Single positive integer giving the total number of scans, the
#'   `n_scans` element returned by [bin_peaks()].
#'
#' @return A list with `mz_list` and `int_list` (each length `n_scans`, the m/z
#'   values and intensities in the bin at each scan) and `eic` (a length
#'   `n_scans` numeric vector of the per-scan mean intensity, zero where the bin
#'   is empty). The raw EIC is suitable for [smooth_intensity()], and `mz_list`
#'   and `int_list` for [collect_zoi_masses()].
#'
#' @examples
#' binned <- bin_peaks(
#'   list(c(100.01, 100.06), c(100.02, 100.07)),
#'   list(c(10, 20), c(40, 50)),
#'   bin_width = 0.05,
#'   mz_range = c(100.00, 100.10)
#' )
#' assemble_bin_traces(binned$bins[[1]], binned$n_scans)
#'
#' @export
assemble_bin_traces <- function(bin, n_scans) {
  if (!is.list(bin) || !all(c("scan", "mz", "intensity") %in% names(bin))) {
    stop(
      "`bin` must be a list with `scan`, `mz`, and `intensity` elements.",
      call. = FALSE
    )
  }
  check_count(n_scans, "n_scans")
  scan <- bin$scan
  mz <- bin$mz
  intensity <- bin$intensity
  if (!is.numeric(scan) || !is.numeric(mz) || !is.numeric(intensity)) {
    stop(
      "`scan`, `mz`, and `intensity` must be numeric.",
      call. = FALSE
    )
  }
  if (length(scan) != length(mz) || length(scan) != length(intensity)) {
    stop(
      "`scan`, `mz`, and `intensity` must have the same length.",
      call. = FALSE
    )
  }
  if (length(scan) > 0L) {
    if (any(!is.finite(scan)) ||
        any(scan != as.integer(scan)) ||
        any(scan < 1) ||
        any(scan > n_scans)) {
      stop(
        "`scan` must contain integer indices within [1, n_scans].",
        call. = FALSE
      )
    }
    if (any(!is.finite(mz)) || any(!is.finite(intensity))) {
      stop(
        "`mz` and `intensity` must be finite.",
        call. = FALSE
      )
    }
  }

  mz_list <- rep(list(numeric(0)), n_scans)
  int_list <- rep(list(numeric(0)), n_scans)
  eic <- numeric(n_scans)
  if (length(scan) > 0L) {
    by_mz <- split(mz, scan)
    by_int <- split(intensity, scan)
    present <- as.integer(names(by_mz))
    mz_list[present] <- by_mz
    int_list[present] <- by_int
    eic[present] <- vapply(by_int, mean, numeric(1))
  }
  list(mz_list = mz_list, int_list = int_list, eic = eic)
}
