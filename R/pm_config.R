# pm_config.R

#' Configuration for a Paramounter analysis
#'
#' Bundles every setting that controls parameter estimation into one validated
#' object, which the analysis threads through to each step. The numeric settings
#' default to the published Paramounter values, so `pm_config()` returns the
#' standard configuration. The exception is `legacy`, which defaults to `FALSE`:
#' the package applies its corrections to the original method unless you ask it
#' not to.
#'
#' @param bin_width Positive m/z bin width for [bin_peaks()] (default `0.05`).
#' @param smooth_half_window Non-negative integer smoothing half-window for
#'   [smooth_intensity()] (default `0`, no smoothing).
#' @param noise_block_size Positive integer block size for the rank-based noise
#'   estimate in [estimate_noise()] (default `10`).
#' @param noise_sd_factor Non-negative multiplier on the noise-floor standard
#'   deviation in [estimate_noise()] (default `3`).
#' @param mass_sd_range Positive multiplier on the mass standard deviation when
#'   converting to a tolerance in [zoi_features()] (default `2`).
#' @param ppm_cutoff Positive relative-tolerance cutoff for discarding noisy
#'   zones, or `NULL` (the default) to determine it automatically from the data
#'   at `ppm_quantile`.
#' @param ppm_quantile Single number in `(0, 1]` giving the quantile used to set
#'   the ppm cutoff automatically when `ppm_cutoff` is `NULL` (default `0.95`).
#' @param min_masses,max_masses Positive integers bounding the number of masses a
#'   zone must contain to be measured (defaults `2` and `199`).
#' @param isolation_gap Non-negative minimum retention-time separation for a zone
#'   to count as isolated in [flag_isolated_zoi()] (default `300`).
#' @param match_mz_tol,match_rt_tol Positive m/z and retention-time half-windows
#'   for cross-file matching in [match_zoi_across_files()] (defaults `0.015` and
#'   `30`).
#' @param trim Single number in `(0, 1]` giving the fraction of each parameter
#'   distribution kept when trimming outliers (default `0.97`).
#' @param legacy Single `TRUE` or `FALSE`. When `FALSE` (the default) the
#'   pipeline uses the corrected behaviour throughout. Set `TRUE` to reproduce
#'   the original Paramounter end to end, including its known quirks: it drives
#'   the scan-count, isolation, first-match, and grouping-bandwidth choices
#'   together. Use it to compare against the published values, not for new
#'   analyses — see the README for what each choice changes.
#'
#' @return A `pm_config` object.
#'
#' @include validate.R
#'
#' @examples
#' pm_config()
#' pm_config(bin_width = 0.02, smooth_half_window = 4L, legacy = TRUE)
#' pm_config(ppm_cutoff = 20) # fixed ppm cutoff instead of automatic
#'
#' @export
pm_config <- new_class(
  "pm_config",
  properties = list(
    bin_width = new_property(
      class_any,
      default = 0.05,
      validator = function(value) as_message(check_positive_number(value, "bin_width"))
    ),
    smooth_half_window = new_property(
      class_any,
      default = 0L,
      validator = function(value) as_message(check_nonneg_integer(value, "smooth_half_window"))
    ),
    noise_block_size = new_property(
      class_any,
      default = 10L,
      validator = function(value) as_message(check_count(value, "noise_block_size"))
    ),
    noise_sd_factor = new_property(
      class_any,
      default = 3,
      validator = function(value) as_message(check_nonneg_number(value, "noise_sd_factor"))
    ),
    mass_sd_range = new_property(
      class_any,
      default = 2,
      validator = function(value) as_message(check_positive_number(value, "mass_sd_range"))
    ),
    ppm_cutoff = new_property(
      class_any,
      default = NULL,
      validator = function(value) {
        if (is.null(value)) {
          NULL
        } else {
          as_message(check_positive_number(value, "ppm_cutoff"))
        }
      }
    ),
    ppm_quantile = new_property(
      class_any,
      default = 0.95,
      validator = function(value) {
        as_message(check_scalar(
          value,
          "ppm_quantile",
          "a single number in (0, 1]",
          min = 0,
          min_strict = TRUE,
          max = 1
        ))
      }
    ),
    min_masses = new_property(
      class_any,
      default = 2L,
      validator = function(value) as_message(check_count(value, "min_masses"))
    ),
    max_masses = new_property(
      class_any,
      default = 199L,
      validator = function(value) as_message(check_count(value, "max_masses"))
    ),
    isolation_gap = new_property(
      class_any,
      default = 300,
      validator = function(value) as_message(check_nonneg_number(value, "isolation_gap"))
    ),
    match_mz_tol = new_property(
      class_any,
      default = 0.015,
      validator = function(value) as_message(check_positive_number(value, "match_mz_tol"))
    ),
    match_rt_tol = new_property(
      class_any,
      default = 30,
      validator = function(value) as_message(check_positive_number(value, "match_rt_tol"))
    ),
    trim = new_property(
      class_any,
      default = 0.97,
      validator = function(value) {
        as_message(check_scalar(
          value,
          "trim",
          "a single number in (0, 1]",
          min = 0,
          min_strict = TRUE,
          max = 1
        ))
      }
    ),
    legacy = new_property(
      class_any,
      default = FALSE,
      validator = function(value) as_message(check_flag(value, "legacy"))
    )
  ),
  validator = function(self) {
    if (is.numeric(self@min_masses) &&
        is.numeric(self@max_masses) &&
        length(self@min_masses) == 1L &&
        length(self@max_masses) == 1L &&
        isTRUE(self@min_masses > self@max_masses)) {
      "`min_masses` must not exceed `max_masses`."
    } else {
      NULL
    }
  }
)

#' @importFrom S7 method
method(print, pm_config) <- function(x, ...) {
  fmt <- function(v) if (is.null(v)) "auto" else format(v)
  cat("<pm_config>\n")
  cat("  binning / smoothing\n")
  cat(sprintf("    bin_width          %s\n", fmt(x@bin_width)))
  cat(sprintf("    smooth_half_window %s\n", fmt(x@smooth_half_window)))
  cat("  noise\n")
  cat(sprintf("    noise_block_size   %s\n", fmt(x@noise_block_size)))
  cat(sprintf("    noise_sd_factor    %s\n", fmt(x@noise_sd_factor)))
  cat("  mass tolerance\n")
  cat(sprintf("    mass_sd_range      %s\n", fmt(x@mass_sd_range)))
  cat(sprintf("    ppm_cutoff         %s\n", fmt(x@ppm_cutoff)))
  cat(sprintf("    ppm_quantile       %s\n", fmt(x@ppm_quantile)))
  cat("  zoi filters\n")
  cat(sprintf("    min_masses         %s\n", fmt(x@min_masses)))
  cat(sprintf("    max_masses         %s\n", fmt(x@max_masses)))
  cat("  isolation / matching\n")
  cat(sprintf("    isolation_gap      %s\n", fmt(x@isolation_gap)))
  cat(sprintf("    match_mz_tol       %s\n", fmt(x@match_mz_tol)))
  cat(sprintf("    match_rt_tol       %s\n", fmt(x@match_rt_tol)))
  cat("  aggregation / reproduction\n")
  cat(sprintf("    trim               %s\n", fmt(x@trim)))
  cat(sprintf("    legacy             %s\n", fmt(x@legacy)))
  invisible(x)
}
