# paramounter.R

#' Measure universal LC-MS parameters from data files
#'
#' Runs the complete Paramounter analysis over a set of centroided
#' mass-spectrometry data files.
#'
#' Each file is read and measured in turn. Its peaks are binned, and every mass
#' bin's noise, zones of interest, mass tolerance, peak width, signal-to-noise,
#' and height are determined. The per-file results are then pooled and trimmed
#' into a single [universal_parameters] object. Given two or more files, the
#' estimated instrument mass and retention-time shifts are included too.
#'
#' Translate that object to software-specific settings with [to_xcms()],
#' [to_msdial()], or [to_mzmine()].
#'
#' @param files Character vector of paths to the data files to analyse (mzML,
#'   mzXML, or CDF). At least two files are needed to estimate instrument shifts;
#'   with a single file the shift distributions are empty.
#' @param config A [pm_config] giving the analysis settings (default
#'   `pm_config()`: the published numeric settings, with the corrections to the
#'   original method applied). Pass `pm_config(legacy = TRUE)` to reproduce the
#'   original end to end instead.
#' @param reader Function used to read one file into per-scan traces, called as
#'   `reader(file)` and returning a list with `mz`, `intensity`, and `rtime`.
#'   Defaults to [read_ms_data()]; override it to read a format the default
#'   reader does not handle, or to supply pre-processed data.
#'
#' @return A [universal_parameters] object holding the measured parameter
#'   distributions, their summary, and the files and configuration used.
#'
#' @examples
#' \dontrun{
#' files <- list.files("demo_data", pattern = "\\.mzXML$", full.names = TRUE)
#' params <- paramounter(files)
#' params@summary
#' }
#'
#' @export
paramounter <- function(files, config = pm_config(), reader = read_ms_data) {
  if (!is.character(files) || length(files) == 0L || any(is.na(files))) {
    stop("`files` must be a non-empty character vector of file paths.", call. = FALSE)
  }
  missing <- !file.exists(files)
  if (any(missing)) {
    stop(
      sprintf("File(s) not found: %s", paste(files[missing], collapse = ", ")),
      call. = FALSE
    )
  }
  check_s7(config, pm_config, "config")
  if (!is.function(reader)) {
    stop("`reader` must be a function.", call. = FALSE)
  }

  per_file <- lapply(files, function(file) {
    traces <- reader(file)
    measure_file(traces$mz, traces$intensity, traces$rtime, config)
  })
  aggregate_files(per_file, files, config)
}
