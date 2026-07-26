# to_msdial.R

#' @noRd
msdial_table <- function(params) {
  est <- point_estimates(params@distributions)
  missing <- character(0)
  if (is.na(est$max_mz_diff)) missing <- c(missing, "mz_diff")
  if (is.na(est$min_peak_height)) missing <- c(missing, "height")
  if (is.na(est$min_peak_scan)) missing <- c(missing, "width_scans")
  if (length(missing) > 0L) {
    stop(
      sprintf("Cannot derive MS-DIAL parameters: no measurements for %s.", paste(missing, collapse = ", ")),
      call. = FALSE
    )
  }
  parameter <- c(
    "mass accuracy: MS1 tolerance (Da)",
    "peak detection: minimum peak height",
    "peak detection: mass slice width (Da)",
    "peak detection: minimum peak width (scans)"
  )
  value <- c(est$max_mz_diff, est$min_peak_height, est$max_mz_diff, est$min_peak_scan)
  if (!is.na(est$max_mass_shift)) {
    parameter <- c(parameter, "alignment: MS1 tolerance (Da)", "alignment: retention time tolerance (min)")
    value <- c(value, est$max_mass_shift, est$max_rt_shift / 60)
  }
  data.frame(parameter = parameter, value = round(value, 3), stringsAsFactors = FALSE)
}

#' Translate universal parameters to MS-DIAL settings
#'
#' Converts a [universal_parameters] object into the MS-DIAL parameter values,
#' following the Paramounter method: the MS1 tolerance and mass slice width are
#' the rounded-up maximum absolute mass tolerance (Da); the minimum peak height
#' is the rounded-down minimum height; the minimum peak width is the rounded-down
#' minimum width in scans; and, when instrument shifts are available (two or more
#' files), the alignment MS1 and retention-time tolerances are the maximum mass
#' shift and the maximum retention-time shift in minutes.
#'
#' @param params A [universal_parameters] object from [paramounter()].
#' @param file Optional path; when given, the table is also written there as CSV.
#'
#' @return A data frame with columns `parameter` and `value`, the MS-DIAL
#'   settings to enter (returned invisibly is not used; the value is returned so
#'   it can be inspected or written).
#'
#' @examples
#' \dontrun{
#' params <- paramounter(files)
#' to_msdial(params, file = "msdial_parameters.csv")
#' }
#'
#' @export
to_msdial <- function(params, file = NULL) {
  if (!S7::S7_inherits(params, universal_parameters)) {
    stop("`params` must be a universal_parameters object.", call. = FALSE)
  }
  table <- msdial_table(params)
  if (!is.null(file)) {
    utils::write.csv(table, file, row.names = FALSE)
  }
  table
}
