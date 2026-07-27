# measure_file.R

#' Measure universal-parameter contributions from one file
#'
#' Runs the full per-file measurement chain over one file's extracted per-scan
#' peak data: it bins the peaks, and for every mass bin with signal it estimates
#' the noise floor and locates the zones of interest, measuring each zone's mass
#' tolerance, peak width, signal-to-noise, and height. This is the per-file step
#' of the analysis, called once per file by the orchestrator; [read_ms_data()]
#' supplies its input and the aggregation across files consumes its output.
#'
#' The relative-tolerance (ppm) cutoff is not applied here. Every zone clearing
#' the mass-count bounds is returned, so its ppm value can enter the pooled
#' distribution the cutoff is derived from; the cutoff is applied once, across
#' all files, during aggregation. Each returned zone carries an `isolated` flag
#' from [flag_isolated_zoi()] marking whether it can seed cross-file shift
#' estimation.
#'
#' The mass-trace outlier settings (`min_points`, `fence` of
#' [collect_zoi_masses()]) are left at their defaults and are not exposed on
#' [pm_config]; they can be added there later if they need to be tunable.
#'
#' @param mz List with one numeric vector of m/z values per scan, as in the `mz`
#'   element returned by [read_ms_data()].
#' @param intensity List of intensities parallel to `mz`.
#' @param rtime Numeric vector of scan retention times, one per scan (the same
#'   length as `mz`).
#' @param config A [pm_config] giving the analysis settings (default
#'   `pm_config()`).
#'
#' @return A named list with `noise` (a numeric vector holding the noise cutoff
#'   of each mass bin that contained signal) and `zoi` (a data frame with one row
#'   per measured zone and columns `ppm`, `mz_diff`, `width_seconds`,
#'   `width_scans`, `sn`, `height`, `reference_mz`, `apex_rt`, and `isolated`).
#'
#' @noRd
measure_file <- function(mz, intensity, rtime, config = pm_config()) {
  if (!is.list(mz) || !is.list(intensity)) {
    stop("`mz` and `intensity` must be lists.", call. = FALSE)
  }
  if (length(mz) != length(intensity)) {
    stop("`mz` and `intensity` must have the same length.", call. = FALSE)
  }
  check_numeric_vector(rtime, "rtime")
  if (length(rtime) != length(mz)) {
    stop("`rtime` must have one value per scan (the same length as `mz`).", call. = FALSE)
  }
  check_s7(config, pm_config, "config")

  binned <- bin_peaks(mz, intensity, bin_width = config@bin_width)
  n_scans <- binned$n_scans
  noise_values <- numeric(0)
  rows <- list()

  for (bin in binned$bins) {
    traces <- assemble_bin_traces(bin, n_scans)
    if (all(traces$eic == 0)) {
      next
    }
    smoothed <- smooth_intensity(traces$eic, config@smooth_half_window)
    noise <- estimate_noise(smoothed, config@noise_block_size, config@noise_sd_factor)
    noise_values <- c(noise_values, noise$cutoff)
    zois <- find_zoi(smoothed, noise$cutoff)
    if (nrow(zois) == 0L) {
      next
    }
    isolated <- flag_isolated_zoi(
      rtime[zois$apex],
      min_gap = config@isolation_gap,
      legacy_isolation = config@legacy
    )
    for (z in seq_len(nrow(zois))) {
      apex <- zois$apex[z]
      apex_peaks <- traces$int_list[[apex]]
      if (length(apex_peaks) == 0L) {
        next
      }
      reference_mz <- traces$mz_list[[apex]][which.max(apex_peaks)]
      walk <- collect_zoi_masses(
        traces$mz_list,
        traces$int_list,
        smoothed,
        noise$cutoff,
        apex,
        reference_mz
      )
      features <- zoi_features(
        walk$masses,
        reference_mz,
        walk$left_edge,
        walk$right_edge,
        apex,
        rtime,
        smoothed[apex],
        noise$noise_mean,
        noise$noise_sd,
        noise$cutoff,
        sd_range = config@mass_sd_range,
        ppm_cutoff = Inf,
        min_masses = config@min_masses,
        max_masses = config@max_masses,
        legacy_scan_count = config@legacy
      )
      if (is.null(features)) {
        next
      }
      rows[[length(rows) + 1L]] <- data.frame(
        ppm = features$ppm,
        mz_diff = features$mz_diff,
        width_seconds = features$width_seconds,
        width_scans = features$width_scans,
        sn = features$sn,
        height = features$height,
        reference_mz = reference_mz,
        apex_rt = rtime[apex],
        isolated = isolated[z]
      )
    }
  }

  zoi <- if (length(rows)) {
    do.call(rbind, rows)
  } else {
    data.frame(
      ppm = numeric(0),
      mz_diff = numeric(0),
      width_seconds = numeric(0),
      width_scans = integer(0),
      sn = numeric(0),
      height = numeric(0),
      reference_mz = numeric(0),
      apex_rt = numeric(0),
      isolated = logical(0)
    )
  }
  list(noise = noise_values, zoi = zoi)
}
