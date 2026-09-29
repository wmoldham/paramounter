# pm_config.R

#' Build a property validator from a `check_*` function
#'
#' Every `pm_config` property validates the same way: run a `check_*` helper and
#' return the `NULL`-or-string form S7 expects. Without this, every property
#' carried its own copy of that closure.
#'
#' The property name is still passed as a string, because an S7 validator
#' receives only the value and has no way to learn which property it guards.
#'
#' @param check A `check_*` function taking `(value, name)`.
#' @param name The property name, used in the error message.
#'
#' @return A validator function suitable for `new_property()`.
#' @noRd
# Must stay above the pm_config roxygen block. A helper placed between that block
# and new_class() merges the two, and this @noRd then suppresses the class's own
# help page.
checked <- function(check, name) {
  force(check)
  force(name)
  function(value) as_message(check(value, name))
}

#' Configuration for a Paramounter analysis
#'
#' Bundles every setting that controls parameter estimation into one validated
#' object, which the analysis threads through to each step. The numeric settings
#' default to the published Paramounter values, so `pm_config()` returns the
#' standard configuration. The exception is `legacy`, which defaults to `FALSE`:
#' the package applies its corrections to the original method unless you ask it
#' not to.
#'
#' @param bin_width Single positive number. The width of the mass bins, in Da
#'   (default `0.05`).
#' @param smooth_half_window Single non-negative integer. The half-window of the
#'   moving average applied to each chromatogram (default `0`, no smoothing).
#' @param noise_block_size Single positive integer. The noise estimate walks the
#'   sorted intensities this many at a time (default `10`).
#' @param noise_sd_factor Single non-negative number. How many standard
#'   deviations above the running mean an intensity must sit to count as signal
#'   rather than noise (default `3`).
#' @param mass_sd_range Single positive number. The mass tolerance reported for a
#'   zone is this many standard deviations of its collected masses (default `2`).
#' @param ppm_cutoff Single positive number. Zones whose relative mass tolerance
#'   reaches this value are discarded as noise. `NULL`, the default, takes the
#'   cutoff from the data at `ppm_quantile` instead.
#' @param ppm_quantile Single number in `(0, 1]`. The quantile of the pooled ppm
#'   distribution used as the cutoff when `ppm_cutoff` is `NULL` (default
#'   `0.95`).
#' @param min_masses,max_masses Single positive integers bounding how many masses
#'   a zone must contain to be measured (defaults `2` and `199`).
#' @param isolation_gap Single non-negative number. The retention-time
#'   separation, in seconds, a zone needs from its neighbours to count as
#'   isolated (default `300`, five minutes).
#' @param match_mz_tol,match_rt_tol Single positive numbers. The half-windows for
#'   matching zones across files: m/z in Da, retention time in seconds (defaults
#'   `0.015` and `30`).
#' @param trim Single number in `(0, 1]`. The fraction of each *threshold*
#'   distribution kept — `noise`, `width_scans`, `sn` and `height`. The rest is
#'   dropped from one end, whichever is opposite the statistic that quantity is
#'   summarised by (default `0.97`). Under `legacy = TRUE` this applies to every
#'   distribution, as in the original.
#' @param tolerance_quantile Single number in `(0, 1]`. The position taken from
#'   each *tolerance* distribution — `ppm`, `mz_diff`, and the upper peak-width
#'   bound (default `0.95`). Ignored under `legacy = TRUE`.
#'
#'   The split matters because the two kinds of parameter fail differently. A
#'   loose threshold only lets more candidates through, and downstream filtering
#'   removes what does not reproduce. A loose tolerance changes how peaks are
#'   built, which cannot be undone later: on a high-dynamic-range Orbitrap
#'   dataset an over-wide `ppm` doubled the chromatographic peaks and returned
#'   *fewer* features. Lower this when [check_estimates()] shows a tolerance
#'   estimate sitting far above the bulk of its distribution.
#'
#'   This is taken directly from the untrimmed distribution rather than stacked
#'   on `trim`. The two would otherwise compose: trimming to `0.97` and then
#'   taking `0.95` of what remains lands at `0.9215`, not `0.95`.
#' @param tolerance_isolated Single `TRUE` or `FALSE`. When `TRUE`, the
#'   *tolerance* distributions — `ppm`, `mz_diff` and `width_seconds` — are
#'   measured only on isolated zones, those with no neighbour in the same mass bin
#'   within `isolation_gap` (default `FALSE`, every zone). Ignored under
#'   `legacy = TRUE`.
#' @param tolerance_min_scans Single positive integer. The *tolerance*
#'   distributions are measured only on zones at least this many scans wide
#'   (default `1L`, every zone). Ignored under `legacy = TRUE`.
#'
#'   The two settings answer the same problem and are meant to be used together.
#'   On some data most measured zones are not chromatographic peaks at all: on a
#'   polarity-switching HILIC run the median zone was 2 scans (1.2 s) wide, and
#'   only 230 of 27,654 were isolated. Those zones set the tolerances. `ppm` came
#'   out at 45 at the 95th percentile in every intensity decile, so no quantile
#'   rescues it, and the upper peak-width bound (27 s) fell below the median width
#'   of a real peak (~25 s), which made `CentWave` split broad peaks. Restricted to
#'   isolated zones at least 5 scans wide, the same file gave a 95th-percentile
#'   `ppm` of 2.7 and peak widths of 5-104 s. Isolation alone was not enough - the
#'   isolated set still had a 95th-percentile `ppm` of 30 - because an isolated
#'   2-scan blip is still a blip.
#'
#'   Only tolerances are restricted. The *threshold* distributions (`noise`,
#'   `sn`, `height`, `width_scans`) keep every zone, because a loose threshold is
#'   recoverable downstream and a loose tolerance is not; see
#'   `tolerance_quantile`. Use [check_estimates()] to see whether a tolerance is
#'   being set by the bulk of its distribution or by its tail.
#' @param legacy Single `TRUE` or `FALSE`. When `FALSE` (the default) the
#'   pipeline applies its corrections to the original method. When `TRUE` it
#'   reproduces the original end to end, setting the scan-count, isolation,
#'   first-match, and grouping-bandwidth behaviour together. Use it to compare
#'   against published values. The README tabulates what each one changes.
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
  # class_any throughout: the check_* helpers do the type work and give better
  # messages than S7's own class dispatch would
  properties = list(
    bin_width = new_property(
      class_any,
      default = 0.05,
      validator = checked(check_positive_number, "bin_width")
    ),
    smooth_half_window = new_property(
      class_any,
      default = 0L,
      validator = checked(check_nonneg_integer, "smooth_half_window")
    ),
    noise_block_size = new_property(
      class_any,
      default = 10L,
      validator = checked(check_count, "noise_block_size")
    ),
    noise_sd_factor = new_property(
      class_any,
      default = 3,
      validator = checked(check_nonneg_number, "noise_sd_factor")
    ),
    mass_sd_range = new_property(
      class_any,
      default = 2,
      validator = checked(check_positive_number, "mass_sd_range")
    ),
    # the only optional setting: NULL means "derive the cutoff from the data"
    ppm_cutoff = new_property(
      class_any,
      default = NULL,
      validator = checked(
        function(value, name) {
          if (!is.null(value)) check_positive_number(value, name)
        },
        "ppm_cutoff"
      )
    ),
    ppm_quantile = new_property(
      class_any,
      default = 0.95,
      validator = checked(check_unit_fraction, "ppm_quantile")
    ),
    min_masses = new_property(
      class_any,
      default = 2L,
      validator = checked(check_count, "min_masses")
    ),
    max_masses = new_property(
      class_any,
      default = 199L,
      validator = checked(check_count, "max_masses")
    ),
    isolation_gap = new_property(
      class_any,
      default = 300,
      validator = checked(check_nonneg_number, "isolation_gap")
    ),
    match_mz_tol = new_property(
      class_any,
      default = 0.015,
      validator = checked(check_positive_number, "match_mz_tol")
    ),
    match_rt_tol = new_property(
      class_any,
      default = 30,
      validator = checked(check_positive_number, "match_rt_tol")
    ),
    trim = new_property(
      class_any,
      default = 0.97,
      validator = checked(check_unit_fraction, "trim")
    ),
    tolerance_quantile = new_property(
      class_any,
      default = 0.95,
      validator = checked(check_unit_fraction, "tolerance_quantile")
    ),
    tolerance_isolated = new_property(
      class_any,
      default = FALSE,
      validator = checked(check_flag, "tolerance_isolated")
    ),
    tolerance_min_scans = new_property(
      class_any,
      default = 1L,
      validator = checked(check_count, "tolerance_min_scans")
    ),
    legacy = new_property(
      class_any,
      default = FALSE,
      validator = checked(check_flag, "legacy")
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
  cat(sprintf("    tolerance_quantile %s\n", fmt(x@tolerance_quantile)))
  cat(sprintf("    tolerance_isolated %s\n", fmt(x@tolerance_isolated)))
  cat(sprintf("    tolerance_min_scans %s\n", fmt(x@tolerance_min_scans)))
  cat(sprintf("    legacy             %s\n", fmt(x@legacy)))
  invisible(x)
}
