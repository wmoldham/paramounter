# to_msdial.R

#' @noRd
msdial_table <- function(params) {
  est <- point_estimates(params@distributions)
  require_estimates(
    est,
    c(
      mz_diff = "max_mz_diff",
      height = "min_peak_height",
      width_scans = "min_peak_scan"
    ),
    "MS-DIAL"
  )
  parameter_table(
    parameter = c(
      "mass accuracy: MS1 tolerance (Da)",
      "peak detection: minimum peak height",
      "peak detection: mass slice width (Da)",
      "peak detection: minimum peak width (scans)"
    ),
    value = c(
      est$max_mz_diff, est$min_peak_height, est$max_mz_diff, est$min_peak_scan
    ),
    est = est,
    shift_labels = c(
      "alignment: MS1 tolerance (Da)",
      "alignment: retention time tolerance (min)"
    )
  )
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
  check_s7(params, universal_parameters, "params")
  write_table_maybe(msdial_table(params), file)
}
