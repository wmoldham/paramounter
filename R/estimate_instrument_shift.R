# estimate_instrument_shift.R

#' Estimate instrument mass and retention-time shift
#'
#' Summarises the matched features from [match_zoi_across_files()] into the two
#' shift distributions that set the cross-file alignment tolerances: for each
#' feature, the spread of its m/z across files (`max - min`) and the spread of
#' its retention time across files. The maxima of these two distributions set the
#' alignment windows for the target software. Trimming and the maximum are
#' applied downstream, alongside the other distributions, so this function
#' returns the untrimmed per-feature values.
#'
#' @param matched List of two numeric matrices `mz` and `rt`, as returned by
#'   [match_zoi_across_files()]: one row per matched feature, one column per
#'   file, with matching dimensions and at least two columns. Values must be
#'   finite.
#'
#' @return A named list with elements `mass_shift` and `rt_shift`, each a numeric
#'   vector with one value per matched feature giving the across-file range of
#'   m/z and of retention time respectively. Both are empty when `matched` has no
#'   rows.
#'
#' @examples
#' matched <- list(
#'   mz = rbind(c(200.000, 200.003), c(300.000, 300.002)),
#'   rt = rbind(c(500, 505), c(800, 803))
#' )
#' estimate_instrument_shift(matched)
#'
#' @export
estimate_instrument_shift <- function(matched) {
  if (!is.list(matched) || !all(c("mz", "rt") %in% names(matched))) {
    stop(
      "`matched` must be a list with `mz` and `rt` elements.",
      call. = FALSE
    )
  }
  mz <- matched$mz
  rt <- matched$rt
  if (!is.matrix(mz) || !is.numeric(mz)) {
    stop(
      "`matched$mz` must be a numeric matrix.",
      call. = FALSE
    )
  }
  if (!is.matrix(rt) || !is.numeric(rt)) {
    stop(
      "`matched$rt` must be a numeric matrix.",
      call. = FALSE
    )
  }
  if (!identical(dim(mz), dim(rt))) {
    stop(
      "`matched$mz` and `matched$rt` must have the same dimensions.",
      call. = FALSE
    )
  }
  if (ncol(mz) < 2L) {
    stop(
      "`matched` must contain at least two files (columns).",
      call. = FALSE
    )
  }
  if (nrow(mz) > 0L && (any(!is.finite(mz)) || any(!is.finite(rt)))) {
    stop(
      "`matched$mz` and `matched$rt` must not contain NA, NaN, or infinite values.",
      call. = FALSE
    )
  }

  row_span <- function(m) {
    if (nrow(m) == 0L) {
      numeric(0)
    } else {
      apply(m, 1, function(row) diff(range(row)))
    }
  }
  list(
    mass_shift = row_span(mz),
    rt_shift = row_span(rt)
  )
}
