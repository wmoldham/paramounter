# bin_peaks.R

#' Assign mass-spectral peaks to m/z bins
#'
#' Slices the m/z axis into fixed-width bins and groups every peak from every
#' scan into its bin, the first step of extracting per-bin chromatograms. Each
#' peak is placed in the half-open bin `[low, high)` that contains it, matching
#' the original Paramounter binning. Bins that receive no peaks are omitted.
#'
#' Unlike the original implementation, which rescans every peak for every bin,
#' this assigns each peak to its bin once, so the cost grows with the number of
#' peaks rather than with bins times scans. Reconstruct a single bin's per-scan
#' traces with [assemble_bin_traces()].
#'
#' The bin edges are `seq(mz_range[1], mz_range[2], by = bin_width)`. A peak
#' lying on or above the final edge falls outside every bin and is dropped; with
#' the default range this affects only the single highest peak, as in the
#' original.
#'
#' @param mz_list List with one element per scan, each a numeric vector of the
#'   m/z values in that scan (possibly empty). Elements must be finite.
#' @param int_list List parallel to `mz_list` giving the corresponding
#'   intensities; element `i` must have the same length as `mz_list[[i]]`.
#' @param bin_width Single positive bin width in m/z units (default `0.05`).
#' @param mz_range Optional length-2 numeric `c(low, high)` with `low < high`
#'   giving the m/z range to bin over. When `NULL` (the default) the range spans
#'   the smallest and largest observed peak, as in the original.
#'
#' @return A list with `n_scans` (the number of scans) and `bins`, a list with
#'   one element per non-empty bin in ascending m/z order. Each bin element is a
#'   list with `bin` (the bin index), `mz_low` and `mz_high` (its edges), and the
#'   parallel vectors `scan`, `mz`, and `intensity` holding that bin's peaks in
#'   scan order.
#'
#' @examples
#' mz <- list(c(100.01, 100.06, 100.12), c(100.02, 100.07), numeric(0))
#' intensity <- list(c(10, 20, 30), c(40, 50), numeric(0))
#' bin_peaks(mz, intensity, bin_width = 0.05, mz_range = c(100.00, 100.15))
#'
#' @export
bin_peaks <- function(mz_list, int_list, bin_width = 0.05, mz_range = NULL) {
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
  if (length(mz_list) != length(int_list)) {
    stop(
      "`mz_list` and `int_list` must have the same length.",
      call. = FALSE
    )
  }
  if (!all(lengths(mz_list) == lengths(int_list))) {
    stop(
      "Each element of `mz_list` and `int_list` must have matching length.",
      call. = FALSE
    )
  }
  n_scans <- length(mz_list)
  mz <- unlist(mz_list, use.names = FALSE)
  intensity <- unlist(int_list, use.names = FALSE)
  if (is.null(mz)) {
    mz <- numeric(0)
  }
  if (is.null(intensity)) {
    intensity <- numeric(0)
  }
  if (!is.numeric(mz) || any(!is.finite(mz))) {
    stop(
      "All elements of `mz_list` must be finite numeric vectors.",
      call. = FALSE
    )
  }
  if (!is.numeric(intensity) || any(!is.finite(intensity))) {
    stop(
      "All elements of `int_list` must be finite numeric vectors.",
      call. = FALSE
    )
  }
  check_positive_number(bin_width, "bin_width")
  if (!is.null(mz_range)) {
    if (!is.numeric(mz_range) ||
        length(mz_range) != 2L ||
        any(!is.finite(mz_range)) ||
        mz_range[1] >= mz_range[2]) {
      stop(
        "`mz_range` must be a length-2 numeric vector with `mz_range[1] < mz_range[2]`.",
        call. = FALSE
      )
    }
  }

  if (length(mz) == 0L) {
    return(list(n_scans = n_scans, bins = list()))
  }
  scans <- rep(seq_len(n_scans), lengths(mz_list))
  if (is.null(mz_range)) {
    mz_range <- c(min(mz), max(mz))
  }
  edges <- seq(mz_range[1], mz_range[2], by = bin_width)
  if (length(edges) < 2L) {
    return(list(n_scans = n_scans, bins = list()))
  }

  bin <- findInterval(mz, edges)
  keep <- bin >= 1L & bin <= (length(edges) - 1L)
  scans <- scans[keep]
  mz <- mz[keep]
  intensity <- intensity[keep]
  bin <- bin[keep]

  by_bin <- split(seq_along(bin), bin)
  bins <- lapply(names(by_bin), function(bin_name) {
    rows <- by_bin[[bin_name]]
    b <- as.integer(bin_name)
    list(
      bin = b,
      mz_low = edges[b],
      mz_high = edges[b + 1L],
      scan = scans[rows],
      mz = mz[rows],
      intensity = intensity[rows]
    )
  })
  list(n_scans = n_scans, bins = bins)
}
