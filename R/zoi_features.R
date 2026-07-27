# zoi_features.R

#' Measure universal parameters for a single zone of interest
#'
#' Converts the mass trace and scan boundaries of one zone of interest (as
#' produced by [collect_zoi_masses()]) into the per-ZOI measurements that feed
#' the universal parameters: relative mass tolerance (ppm), absolute mass
#' tolerance (Da), peak width in seconds and in scans, signal-to-noise ratio,
#' and peak height. Two validity filters are applied first, matching the
#' original Paramounter workflow: the mass vector must contain between
#' `min_masses` and `max_masses` values, and the relative mass tolerance must
#' fall below `ppm_cutoff`. A zone that fails either filter yields no measurement
#' and `NULL` is returned.
#'
#' The relative mass tolerance is `sd_range * sd(masses) / reference_mz * 1e6`
#' and the absolute tolerance is `sd_range * sd(masses)`, both representing the
#' mass fluctuation across the peak at the chosen confidence (default `sd_range`
#' of 2, roughly 95%). The signal-to-noise ratio is
#' `(apex_intensity - noise_mean) / noise_sd`; when `noise_sd` is zero it is
#' infinite, as in the original.
#'
#' @param masses Numeric vector of m/z values gathered across the peak, the
#'   `masses` element from [collect_zoi_masses()]. Must be finite.
#' @param reference_mz Single positive reference m/z used to convert the mass
#'   fluctuation to ppm.
#' @param left_edge,right_edge Single integer scan indices of the leftmost and
#'   rightmost points of the peak (`left_edge <= right_edge`), from
#'   [collect_zoi_masses()].
#' @param apex Single integer scan index of the apex, satisfying
#'   `left_edge <= apex <= right_edge`.
#' @param rtime Numeric vector of retention times per scan. Must be finite and
#'   long enough to contain `right_edge`.
#' @param apex_intensity Single non-negative intensity at the apex (the value of
#'   the smoothed extracted-ion chromatogram at `apex`).
#' @param noise_mean,noise_sd Single non-negative mean and standard deviation of
#'   the noise floor, from [estimate_noise()].
#' @param cutoff Single non-negative noise cutoff. When it is zero the
#'   signal-to-noise ratio falls back to `apex_intensity`.
#' @param sd_range Single positive multiplier on the mass standard deviation
#'   (default `2`).
#' @param ppm_cutoff Single positive upper bound on the relative mass tolerance;
#'   zones at or above it are discarded. Defaults to `Inf` (no filtering).
#' @param min_masses,max_masses Single positive integers bounding the number of
#'   collected masses for a zone to be measured (defaults `2` and `199`).
#' @param legacy_scan_count Single `TRUE` or `FALSE`. When `FALSE` (the default)
#'   `width_scans` is the corrected inclusive span `right_edge - left_edge + 1`.
#'   When `TRUE` it reproduces the original Paramounter implementation, which
#'   counts `right_edge - left_edge + 2` and so overstates the span by one; use
#'   that only to compare against the published values. Keep this in step with
#'   [pm_config]'s `legacy`, which is what sets it when the pipeline is run
#'   through [paramounter()].
#'
#' @return A named list with elements `reference_mz`, `apex_rt`, `ppm`,
#'   `mz_diff` (absolute tolerance, Da), `width_seconds`, `width_scans`, `sn`,
#'   and `height`; or `NULL` if the zone fails the mass-count or ppm filter.
#'
#' @examples
#' rtime <- seq(10, 28, by = 2)
#' masses <- c(220.0098, 220.0100, 220.0102, 220.0101, 220.0099)
#' zoi_features(
#'   masses,
#'   reference_mz = 220.0100,
#'   left_edge = 2,
#'   right_edge = 6,
#'   apex = 4,
#'   rtime = rtime,
#'   apex_intensity = 5000,
#'   noise_mean = 100,
#'   noise_sd = 30,
#'   cutoff = 120
#' )
#'
#' @export
zoi_features <- function(
    masses,
    reference_mz,
    left_edge,
    right_edge,
    apex,
    rtime,
    apex_intensity,
    noise_mean,
    noise_sd,
    cutoff,
    sd_range = 2,
    ppm_cutoff = Inf,
    min_masses = 2L,
    max_masses = 199L,
    legacy_scan_count = FALSE
) {
  check_numeric_vector(masses, "masses")
  check_numeric_vector(rtime, "rtime", allow_empty = FALSE)
  n_rt <- length(rtime)
  check_index(left_edge, "left_edge", n_rt, "rtime")
  check_index(right_edge, "right_edge", n_rt, "rtime")
  check_index(apex, "apex", n_rt, "rtime")
  if (left_edge > right_edge) {
    stop(
      "`left_edge` must not exceed `right_edge`.",
      call. = FALSE
    )
  }
  if (apex < left_edge || apex > right_edge) {
    stop(
      "`apex` must lie between `left_edge` and `right_edge`.",
      call. = FALSE
    )
  }
  check_positive_number(reference_mz, "reference_mz")
  check_nonneg_number(apex_intensity, "apex_intensity")
  check_nonneg_number(noise_mean, "noise_mean")
  check_nonneg_number(noise_sd, "noise_sd")
  check_nonneg_number(cutoff, "cutoff")
  check_positive_number(sd_range, "sd_range")
  check_positive_number(ppm_cutoff, "ppm_cutoff", allow_infinite = TRUE)
  check_count(min_masses, "min_masses")
  check_count(max_masses, "max_masses")
  if (min_masses > max_masses) {
    stop(
      "`min_masses` must not exceed `max_masses`.",
      call. = FALSE
    )
  }
  check_flag(legacy_scan_count, "legacy_scan_count")
  left_edge <- as.integer(left_edge)
  right_edge <- as.integer(right_edge)
  apex <- as.integer(apex)

  n_masses <- length(masses)
  if (n_masses < min_masses || n_masses > max_masses) {
    return(NULL)
  }
  mass_sd <- stats::sd(masses)
  ppm <- sd_range * mass_sd / reference_mz * 1e6
  if (!(ppm < ppm_cutoff)) {
    return(NULL)
  }

  list(
    reference_mz = reference_mz,
    apex_rt = rtime[apex],
    ppm = ppm,
    mz_diff = sd_range * mass_sd,
    width_seconds = rtime[right_edge] - rtime[left_edge],
    width_scans = right_edge - left_edge + if (legacy_scan_count) 2L else 1L,
    sn = if (cutoff > 0) (apex_intensity - noise_mean) / noise_sd else apex_intensity,
    height = apex_intensity
  )
}
