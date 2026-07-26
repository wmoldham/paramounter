# to_mzmine.R

#' @noRd
mzmine_table <- function(params) {
  est <- point_estimates(params@distributions)
  missing <- character(0)
  if (is.na(est$min_noise)) missing <- c(missing, "noise")
  if (is.na(est$min_peak_scan)) missing <- c(missing, "width_scans")
  if (is.na(est$min_peak_height)) missing <- c(missing, "height")
  if (is.na(est$max_ppm)) missing <- c(missing, "ppm")
  if (anyNA(est$peakwidth)) missing <- c(missing, "width_seconds/height")
  if (length(missing) > 0L) {
    stop(
      sprintf("Cannot derive MZmine parameters: no measurements for %s.", paste(missing, collapse = ", ")),
      call. = FALSE
    )
  }
  parameter <- c(
    "mass detection: noise level",
    "ADAP chromatogram builder: min group size (scans)",
    "ADAP chromatogram builder: group intensity threshold",
    "ADAP chromatogram builder: min highest intensity",
    "ADAP chromatogram builder: m/z tolerance (ppm)",
    "chromatogram deconvolution: peak duration min (min)",
    "chromatogram deconvolution: peak duration max (min)"
  )
  value <- c(
    est$min_noise, est$min_peak_scan, est$min_noise, est$min_peak_height,
    est$max_ppm, est$peakwidth[1] / 60, est$peakwidth[2] / 60
  )
  if (!is.na(est$max_mass_shift)) {
    parameter <- c(parameter, "alignment: m/z tolerance (Da)", "alignment: RT tolerance (min)")
    value <- c(value, est$max_mass_shift, est$max_rt_shift / 60)
  }
  data.frame(parameter = parameter, value = round(value, 3), stringsAsFactors = FALSE)
}

#' Translate universal parameters to MZmine settings
#'
#' Converts a [universal_parameters] object into the MZmine parameter values,
#' following the Paramounter method: the mass-detection noise level and the ADAP
#' group intensity threshold are the rounded-down minimum noise; the minimum
#' group size is the rounded-down minimum width in scans; the minimum highest
#' intensity is the rounded-down minimum height; the m/z tolerance is the
#' rounded-up maximum relative tolerance (ppm); the peak-duration range is the
#' derived peak-width bounds in minutes; and, when instrument shifts are
#' available, the alignment m/z and retention-time tolerances are the maximum
#' mass shift and maximum retention-time shift in minutes.
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
  if (!S7::S7_inherits(params, universal_parameters)) {
    stop("`params` must be a universal_parameters object.", call. = FALSE)
  }
  table <- mzmine_table(params)
  if (!is.null(file)) {
    utils::write.csv(table, file, row.names = FALSE)
  }
  table
}
