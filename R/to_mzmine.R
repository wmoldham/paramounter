# to_mzmine.R

#' @noRd
mzmine_table <- function(params) {
  est <- point_estimates(params@distributions)
  require_estimates(
    est,
    c(
      noise = "min_noise",
      width_scans = "min_peak_scan",
      height = "min_peak_height",
      ppm = "max_ppm",
      "width_seconds/height" = "peakwidth"
    ),
    "MZmine"
  )
  parameter_table(
    parameter = c(
      "mass detection: noise level",
      "ADAP chromatogram builder: min group size (scans)",
      "ADAP chromatogram builder: group intensity threshold",
      "ADAP chromatogram builder: min highest intensity",
      "ADAP chromatogram builder: m/z tolerance (ppm)",
      "chromatogram deconvolution: peak duration min (min)",
      "chromatogram deconvolution: peak duration max (min)"
    ),
    value = c(
      est$min_noise, est$min_peak_scan, est$min_noise, est$min_peak_height,
      est$max_ppm, est$peakwidth[1] / 60, est$peakwidth[2] / 60
    ),
    est = est,
    shift_labels = c(
      "alignment: m/z tolerance (Da)",
      "alignment: RT tolerance (min)"
    )
  )
}

#' Translate universal parameters to MZmine settings
#'
#' Converts a [universal_parameters] object into the MZmine parameter values,
#' following the Paramounter method. The mass-detection noise level and the ADAP
#' group intensity threshold are the rounded-down minimum noise. The minimum
#' group size is the rounded-down minimum width in scans, and the minimum highest
#' intensity the rounded-down minimum height. The m/z tolerance is the rounded-up
#' maximum relative tolerance in ppm, and the peak-duration range is the derived
#' peak-width bounds in minutes.
#'
#' When the instrument shifts are available, which needs two or more files, the
#' alignment m/z tolerance is the maximum mass shift in Da and the alignment
#' retention-time tolerance is the maximum retention-time shift in minutes.
#'
#' @param params A [universal_parameters] object from [paramounter()].
#' @param file Optional path; when given, the table is also written there as CSV.
#'
#' @return A data frame with columns `parameter` and `value`, the MZmine
#'   settings to enter.
#'
#' @examples
#' \dontrun{
#' params <- paramounter(files)
#' to_mzmine(params, file = "mzmine_parameters.csv")
#' }
#'
#' @export
to_mzmine <- function(params, file = NULL) {
  check_s7(params, universal_parameters, "params")
  write_table_maybe(mzmine_table(params), file)
}
