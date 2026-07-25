# read_ms_data.R

#' Read one LC-MS data file
#'
#' Reads a single mass-spectrometry data file and extracts the per-scan peak
#' data the analysis works from: for each scan at the requested MS level, the m/z
#' values and their intensities, together with the scan retention times. Reading
#' uses the Spectra package with its mzR backend, so any format mzR supports
#' (mzML, mzXML, CDF) can be read.
#'
#' The analysis expects centroided spectra; profile data can be read but should
#' be centroided first for meaningful results. One file is read per call; the
#' orchestrator reads each of several files in turn.
#'
#' @param file Single path to an existing mass-spectrometry data file.
#' @param ms_level Single positive integer giving the MS level to extract
#'   (default `1`).
#'
#' @return A named list with `mz` and `intensity` (each a list holding one
#'   numeric vector per scan, parallel to one another) and `rtime` (a numeric
#'   vector of scan retention times). All three have one entry per scan at the
#'   requested MS level.
#'
#' @examples
#' \dontrun{
#' f <- system.file("microtofq/MM14.mzML", package = "msdata")
#' traces <- read_ms_data(f)
#' length(traces$rtime)
#' }
#'
#' @export
read_ms_data <- function(file, ms_level = 1L) {
  if (!is.character(file) || length(file) != 1L || is.na(file)) {
    stop(
      "`file` must be a single file path.",
      call. = FALSE
    )
  }
  if (!file.exists(file)) {
    stop(
      sprintf("File not found: %s", file),
      call. = FALSE
    )
  }
  check_count(ms_level, "ms_level")
  ms_level <- as.integer(ms_level)

  spectra <- Spectra::Spectra(file, source = Spectra::MsBackendMzR())
  spectra <- Spectra::filterMsLevel(spectra, ms_level)
  if (length(spectra) == 0L) {
    stop(
      sprintf("No MS level %d spectra found in %s", ms_level, file),
      call. = FALSE
    )
  }

  list(
    mz = as.list(Spectra::mz(spectra)),
    intensity = as.list(Spectra::intensity(spectra)),
    rtime = Spectra::rtime(spectra)
  )
}
