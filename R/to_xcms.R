# to_xcms.R

#' Derive xcms parameter values from universal parameters
#'
#' The numeric core of [to_xcms()], separated so it can be tested without xcms
#' installed. Reads the trimmed distributions from a [universal_parameters]
#' object and applies the Paramounter conversions.
#'
#' @param params A [universal_parameters] object.
#' @param sample_groups Optional sample-group vector, one per file.
#'
#' @return A named list of the computed values.
#' @noRd
xcms_values <- function(params, sample_groups = NULL) {
  d <- params@distributions
  required <- c("ppm", "noise", "width_seconds", "width_scans", "sn", "height")
  empty <- required[vapply(d[required], length, integer(1)) == 0L]
  if (length(empty) > 0L) {
    stop(
      sprintf("Cannot derive xcms parameters: no measurements for %s.", paste(empty, collapse = ", ")),
      call. = FALSE
    )
  }
  legacy <- params@config@legacy
  n_files <- length(params@files)

  ppm <- ceiling(max(d$ppm))
  minnoise <- floor(min(d$noise))
  w <- mean(d$width_seconds, trim = 0.05)
  h <- mean(d$height, trim = 0.05)
  ratio <- h / w
  min_pw <- min(d$width_seconds)
  max_pw <- max(d$width_seconds)
  if (max_pw > 35 && ratio > 515) {
    peakwidth <- c(0, (ceiling(max_pw) + 7) / 2)
  } else {
    peakwidth <- c(ceiling(min_pw) + 4, ceiling(max_pw) + 5)
  }
  minpeakscan <- floor(min(d$width_scans))
  snthresh <- max(3, min(d$sn))
  bw <- if (!legacy && length(d$rt_shift) > 0L) max(d$rt_shift) else 5
  binSize <- if (length(d$mass_shift) > 0L) max(d$mass_shift) else 0.012
  if (is.null(sample_groups)) {
    sample_groups <- rep(1L, n_files)
  }

  list(
    ppm = ppm,
    peakwidth = peakwidth,
    snthresh = snthresh,
    prefilter = c(minpeakscan, minnoise),
    mzdiff = -0.01,
    noise = minnoise,
    integrate = 1L,
    bw = bw,
    binSize = binSize,
    minFraction = 0.5,
    minSamples = 1L,
    maxFeatures = 100L,
    obiwarp_binSize = 1,
    sample_groups = sample_groups
  )
}

#' Translate universal parameters to xcms settings
#'
#' Converts a [universal_parameters] object into the parameter objects for an
#' xcms workflow: a `CentWaveParam` for chromatographic peak detection, a
#' `PeakDensityParam` for correspondence (grouping), and an `ObiwarpParam` for
#' retention-time alignment. The conversions follow the Paramounter method:
#' `ppm` is the rounded-up maximum relative tolerance; `noise` and the prefilter
#' intensity are the rounded-down minimum noise; the prefilter peak count is the
#' rounded-down minimum peak width in scans; `snthresh` is the minimum
#' signal-to-noise floored at 3; and the peak-width range is derived from the
#' peak-width (seconds) distribution, with the wide-peak adjustment from the
#' original. Grouping `binSize` is the maximum mass shift; grouping `bw` is the
#' maximum retention-time shift, unless `config@legacy` is `TRUE`, which
#' reproduces the original's fixed `bw = 5`.
#'
#' With fewer than two files the instrument shifts cannot be measured, so `bw`
#' and `binSize` fall back to `5` and `0.012` and a warning is issued.
#'
#' @param params A [universal_parameters] object from [paramounter()].
#' @param sample_groups Optional vector assigning each file to a sample group,
#'   one entry per file, used for the `PeakDensityParam`. Defaults to a single
#'   group; set it to your experimental design.
#'
#' @return A named list with `chrom_peaks` (a `CentWaveParam`), `group` (a
#'   `PeakDensityParam`), and `retention` (an `ObiwarpParam`), ready for
#'   `findChromPeaks()`, `groupChromPeaks()`, and `adjustRtime()`.
#'
#' @examples
#' \dontrun{
#' params <- paramounter(files)
#' xp <- to_xcms(params)
#' data <- findChromPeaks(data, xp$chrom_peaks)
#' data <- groupChromPeaks(data, xp$group)
#' data <- adjustRtime(data, xp$retention)
#' }
#'
#' @export
to_xcms <- function(params, sample_groups = NULL) {
  if (!S7::S7_inherits(params, universal_parameters)) {
    stop("`params` must be a universal_parameters object.", call. = FALSE)
  }
  if (!is.null(sample_groups) && length(sample_groups) != length(params@files)) {
    stop("`sample_groups` must have one entry per file.", call. = FALSE)
  }
  if (length(params@distributions$mass_shift) == 0L) {
    warning(
      "Instrument shifts unavailable (fewer than two files); grouping `bw` and `binSize` fall back to 5 and 0.012.",
      call. = FALSE
    )
  }

  v <- xcms_values(params, sample_groups)
  list(
    chrom_peaks = xcms::CentWaveParam(
      ppm = v$ppm,
      peakwidth = v$peakwidth,
      snthresh = v$snthresh,
      prefilter = v$prefilter,
      mzdiff = v$mzdiff,
      noise = v$noise,
      integrate = v$integrate
    ),
    group = xcms::PeakDensityParam(
      sampleGroups = v$sample_groups,
      bw = v$bw,
      binSize = v$binSize,
      minFraction = v$minFraction,
      minSamples = v$minSamples,
      maxFeatures = v$maxFeatures
    ),
    retention = xcms::ObiwarpParam(binSize = v$obiwarp_binSize)
  )
}
