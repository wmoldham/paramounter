# point_estimates.R

#' Shared point estimates from the parameter distributions
#'
#' Computes the summary values the software translators share, guarding empty
#' distributions with `NA`. The peak-width bounds apply the original's wide-peak
#' adjustment, using the 5%-trimmed-mean height-to-width ratio.
#'
#' @param distributions The `distributions` list of a [universal_parameters].
#' @return A named list of estimates.
#' @noRd
point_estimates <- function(distributions) {
  d <- distributions
  present <- function(x) length(x) > 0L
  if (present(d$width_seconds) && present(d$height)) {
    w <- mean(d$width_seconds, trim = 0.05)
    h <- mean(d$height, trim = 0.05)
    ratio <- h / w
    lo <- min(d$width_seconds)
    hi <- max(d$width_seconds)
    peakwidth <- if (hi > 35 && ratio > 515) {
      c(0, (ceiling(hi) + 7) / 2)
    } else {
      c(ceiling(lo) + 4, ceiling(hi) + 5)
    }
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

#' Require that the estimates a translator needs were actually measured
#'
#' Each translator needs a different subset of `point_estimates()`, and any of
#' them is `NA` when the underlying distribution was never measured. Report all
#' the missing ones at once rather than failing on the first, so a user with an
#' empty run learns everything that is wrong in one go.
#'
#' @param est The list returned by `point_estimates()`.
#' @param required Named character vector mapping the label to report to the
#'   name of the estimate it comes from. Order is preserved in the message.
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
#' and expressions. The alignment tolerances are appended only when the
#' instrument shifts were measured, which needs two or more files.
#'
#' @param parameter,value Parallel vectors of setting names and values.
#' @param est The list returned by `point_estimates()`.
#' @param shift_labels Length-two character vector naming the mass- and
#'   retention-time alignment tolerances, or `NULL` to never append them.
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

#' Return a table, writing it to CSV first when a path was given
#'
#' @param table A data frame.
#' @param file Optional path.
#'
#' @return `table`.
#' @noRd
write_table_maybe <- function(table, file) {
  if (!is.null(file)) {
    utils::write.csv(table, file, row.names = FALSE)
  }
  table
}
