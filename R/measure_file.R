# measure_file.R

# The columns of the per-zone result, as length-one prototypes. Declared once so
# the empty case matches the assembled case. Each prototype doubles as the
# `vapply` FUN.VALUE that type-checks its column as it is built.
ZOI_COLUMNS <- list(
  ppm = NA_real_,
  mz_diff = NA_real_,
  width_seconds = NA_real_,
  width_scans = NA_integer_,
  sn = NA_real_,
  height = NA_real_,
  reference_mz = NA_real_,
  apex_rt = NA_real_,
  isolated = NA
)

#' Assemble the collected per-zone measurements into a data frame
#'
#' Zones are collected as bare lists and turned into columns once, at the end.
#' There are tens of thousands of zones per file, and building a one-row
#' `data.frame()` for each is slow: every call deparses its own column names.
#'
#' @param rows List of per-zone named lists, each holding one value per column
#'   of `ZOI_COLUMNS`.
#'
#' @return A data frame with one row per zone and the columns of `ZOI_COLUMNS`,
#'   with zero rows when `rows` is empty.
#' @noRd
assemble_zoi <- function(rows) {
  columns <- lapply(names(ZOI_COLUMNS), function(column) {
    if (length(rows) == 0L) {
      ZOI_COLUMNS[[column]][0]
    } else {
      vapply(rows, `[[`, ZOI_COLUMNS[[column]], column)
    }
  })
  names(columns) <- names(ZOI_COLUMNS)
  data.frame(columns, stringsAsFactors = FALSE)
}

#' Measure universal-parameter contributions from one file
#'
#' Runs the full measurement chain over one file's extracted per-scan peak data.
#' It bins the peaks. Then, for every mass bin with signal, it estimates the
#' noise floor and locates the zones of interest, measuring each zone's mass
#' tolerance, peak width, signal-to-noise, and height.
#'
#' [read_ms_data()] supplies its input; `aggregate_files()` consumes its output.
#'
#' The relative-tolerance (ppm) cutoff is not applied here. Every zone clearing
#' the mass-count bounds is returned, so its ppm value can enter the pooled
#' distribution the cutoff is derived from; the cutoff is applied once, across
#' all files, during aggregation. Each returned zone carries an `isolated` flag
#' from [flag_isolated_zoi()] marking whether it can seed cross-file shift
#' estimation.
#'
#' The mass-trace outlier settings, `min_points` and `fence` of
#' [collect_zoi_masses()], are left at their defaults and are not exposed on
#' [pm_config].
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
  # both accumulators are filled by index rather than grown with c()/append.
  # There is one noise value per bin and one row per zone, tens of thousands of
  # each, so growing them reallocates on nearly every iteration.
  n_bins <- length(binned$bins)
  noise_values <- numeric(n_bins)
  n_noise <- 0L
  rows <- vector("list", n_bins)
  n_rows <- 0L

  for (bin in binned$bins) {
    traces <- assemble_bin_traces(bin, n_scans)
    if (all(traces$eic == 0)) {
      next
    }
    smoothed <- smooth_intensity(traces$eic, config@smooth_half_window)
    noise <- estimate_noise(smoothed, config@noise_block_size, config@noise_sd_factor)
    n_noise <- n_noise + 1L
    noise_values[n_noise] <- noise$cutoff
    zois <- zoi_indices(smoothed, noise$cutoff)
    n_zoi <- length(zois$apex)
    if (n_zoi == 0L) {
      next
    }
    isolated <- flag_isolated_zoi(
      rtime[zois$apex],
      min_gap = config@isolation_gap,
      legacy_isolation = config@legacy
    )
    for (z in seq_len(n_zoi)) {
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
      n_rows <- n_rows + 1L
      if (n_rows > length(rows)) {
        length(rows) <- 2L * length(rows)
      }
      rows[[n_rows]] <- list(
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

  list(
    noise = noise_values[seq_len(n_noise)],
    zoi = assemble_zoi(rows[seq_len(n_rows)])
  )
}
