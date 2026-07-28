# point_estimates.R

#' Shared point estimates from the parameter distributions
#'
#' Computes the summary values the software translators share, guarding empty
#' distributions with `NA`. The peak-width bounds apply the original's wide-peak
#' adjustment, using the 5%-trimmed-mean height-to-width ratio.
#'
#' The wide-peak branch fires when the widest measured peak exceeds 35 seconds and
#' the height-to-width ratio exceeds 515. It halves the upper bound, which keeps
#' the search from chasing peaks far wider than the bulk of the data. In the
#' original it also drops the lower bound to zero, and that part does not
#' survive contact with a real dataset: `CentWave` turns `peakwidth` into wavelet
#' scales as roughly `peakwidth / mean(diff(scantime)) / 2`, so a lower bound of
#' zero asks it to look for peaks of no width, which cannot be told from a spike.
#' Under `legacy = FALSE` the lower bound therefore stays data-driven in both
#' branches; `legacy = TRUE` returns the zero, which is what Table S-7 published.
#'
#' The `ratio > 515` guard is left as the original wrote it. It compares an
#' intensity-per-second quantity against an absolute constant calibrated on
#' Bruker Q-TOF counts, so on instruments with much larger intensities it is
#' effectively always true and the branch turns on `hi > 35` alone. That is worth
#' knowing, but re-deriving a scale-free replacement would need more than one
#' dataset, and once the lower bound is fixed the branch only halves an upper
#' bound that an outlier drove up.
#'
#' @param distributions Named list of one numeric vector per quantity, the
#'   `distributions` property of a [universal_parameters]. Quantities with no
#'   measurements are empty vectors.
#' @param legacy Single `TRUE` or `FALSE`. When `TRUE` the wide-peak branch
#'   returns a lower bound of zero, reproducing the original.
#' @return A named list of estimates.
#' @noRd
point_estimates <- function(distributions, legacy = FALSE) {
  d <- distributions
  present <- function(x) length(x) > 0L
  if (present(d$width_seconds) && present(d$height)) {
    w <- mean(d$width_seconds, trim = 0.05)
    h <- mean(d$height, trim = 0.05)
    ratio <- h / w
    lo <- min(d$width_seconds)
    hi <- max(d$width_seconds)
    wide <- hi > 35 && ratio > 515
    peakwidth <- c(
      if (wide && legacy) 0 else ceiling(lo) + 4,
      if (wide) (ceiling(hi) + 7) / 2 else ceiling(hi) + 5
    )
  } else {
    peakwidth <- c(NA_real_, NA_real_)
  }
  list(
    max_ppm = if (present(d$ppm)) ceiling(max(d$ppm)) else NA_real_,
    max_mz_diff = if (present(d$mz_diff)) ceiling(max(d$mz_diff) * 100) / 100 else NA_real_,
    min_noise = if (present(d$noise)) floor(min(d$noise)) else NA_real_,
    min_peak_scan = if (present(d$width_scans)) floor(min(d$width_scans)) else NA_real_,
    min_peak_height = if (present(d$height)) floor(min(d$height)) else NA_real_,
    min_sn = if (present(d$sn)) max(3, min(d$sn)) else NA_real_,
    peakwidth = peakwidth,
    max_mass_shift = if (present(d$mass_shift)) max(d$mass_shift) else NA_real_,
    max_rt_shift = if (present(d$rt_shift)) max(d$rt_shift) else NA_real_
  )
}

#' Require that the estimates a translator needs were measured
#'
#' Each translator needs a different subset of `point_estimates()`, and any of
#' them is `NA` when the underlying distribution was never measured. Report all
#' the missing ones rather than failing on the first, so a user with an empty run
#' sees everything at once.
#'
#' @param est The list returned by `point_estimates()`.
#' @param required Named character vector. Names are the labels to report,
#'   values the estimates they come from. Order is preserved in the message.
#' @param software Software name, used in the error message.
#'
#' @return `est`, invisibly.
#' @noRd
require_estimates <- function(est, required, software) {
  absent <- vapply(required, function(field) anyNA(est[[field]]), logical(1))
  if (any(absent)) {
    stop(
      sprintf(
        "Cannot derive %s parameters: no measurements for %s.",
        software,
        paste(names(required)[absent], collapse = ", ")
      ),
      call. = FALSE
    )
  }
  invisible(est)
}

#' Assemble a translator's parameter/value table
#'
#' Shared by [to_msdial()] and [to_mzmine()], which differ only in their labels
#' and the values they compute. The alignment tolerances are appended only when
#' the instrument shifts were measured, which needs two or more files.
#'
#' @param parameter,value Parallel vectors: `parameter` character, `value`
#'   numeric.
#' @param est The list returned by `point_estimates()`.
#' @param shift_labels Length-two character vector naming the mass- and
#'   retention-time alignment tolerances, or `NULL` to omit them.
#'
#' @return A data frame with columns `parameter` and `value`.
#' @noRd
parameter_table <- function(parameter, value, est, shift_labels = NULL) {
  if (!is.null(shift_labels) && !is.na(est$max_mass_shift)) {
    parameter <- c(parameter, shift_labels)
    value <- c(value, est$max_mass_shift, est$max_rt_shift / 60)
  }
  data.frame(
    parameter = parameter,
    value = round(value, 3),
    stringsAsFactors = FALSE
  )
}

#' Write a table to CSV when a path is given
#'
#' @param table A data frame.
#' @param file Optional path to write to. `NULL` writes nothing.
#'
#' @return `table`.
#' @noRd
write_table_maybe <- function(table, file) {
  if (!is.null(file)) {
    utils::write.csv(table, file, row.names = FALSE)
  }
  table
}
