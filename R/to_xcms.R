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
  est <- point_estimates(
    params@distributions,
    params@config@legacy,
    params@config@tolerance_quantile
  )
  require_estimates(
    est,
    c(
      ppm = "max_ppm",
      noise = "min_noise",
      width_scans = "min_peak_scan",
      sn = "min_sn",
      "width_seconds/height" = "peakwidth"
    ),
    "xcms"
  )
  legacy <- params@config@legacy
  bw <- if (!legacy && !is.na(est$max_rt_shift)) est$max_rt_shift else 5
  binSize <- if (!is.na(est$max_mass_shift)) est$max_mass_shift else 0.012
  if (is.null(sample_groups)) {
    sample_groups <- rep(1L, length(params@files))
  }
  list(
    ppm = est$max_ppm,
    peakwidth = est$peakwidth,
    snthresh = est$min_sn,
    prefilter = c(est$min_peak_scan, est$min_noise),
    mzdiff = -0.01,
    noise = est$min_noise,
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
#' When the instrument shifts were never measured, `bw` and `binSize` fall back
#' to `5` and `0.012` and a warning is issued. That happens with a single file,
#' and also when several files share no zone that matched across all of them.
#'
#' @param params A [universal_parameters] object from [paramounter()].
#' @param sample_groups Optional vector assigning each file of the **experiment
#'   you are about to process** to a sample group, one entry per file, passed to
#'   the `PeakDensityParam`. This is normally longer than `params@files`, because
#'   the point of measuring a handful of representative injections is that you do
#'   not have to measure all of them; its length is not checked against them.
#'
#'   The default puts every *measured* file in one group, which is only right
#'   when you are processing exactly the files you measured. Set it explicitly
#'   otherwise.
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
  check_s7(params, universal_parameters, "params")
  # deliberately not checked against length(params@files): the files measured and
  # the files to be processed are different sets, which is the point of measuring
  # a few representative injections
  if (!is.null(sample_groups) && length(sample_groups) == 0L) {
    stop("`sample_groups` must have at least one entry.", call. = FALSE)
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
