# point_estimates.R

#' Shared point estimates from the parameter distributions
#'
#' Computes the summary values the software translators share, guarding empty
#' distributions with `NA`. The peak-width bounds apply the original's wide-peak
#' adjustment, using the 5%-trimmed-mean height-to-width ratio.
#'
#' @param distributions The `distributions` list of a [universal_parameters].
#' @return A named list of estimates.
#' @noRd
point_estimates <- function(distributions) {
  d <- distributions
  present <- function(x) length(x) > 0L
  if (present(d$width_seconds) && present(d$height)) {
    w <- mean(d$width_seconds, trim = 0.05)
    h <- mean(d$height, trim = 0.05)
    ratio <- h / w
    lo <- min(d$width_seconds)
    hi <- max(d$width_seconds)
    peakwidth <- if (hi > 35 && ratio > 515) {
      c(0, (ceiling(hi) + 7) / 2)
    } else {
      c(ceiling(lo) + 4, ceiling(hi) + 5)
    }
  } else {
    peakwidth <- c(NA_real_, NA_real_)
  }
  list(
    max_ppm = if (present(d$ppm)) ceiling(max(d$ppm)) else NA_real_,
    max_mz_diff = if (present(d$mz_diff)) ceiling(max(d$mz_diff) * 100) / 100 else NA_real_,
    min_noise = if (present(d$noise)) floor(min(d$noise)) else NA_real_,
    min_peak_scan = if (present(d$width_scans)) floor(min(d$width_scans)) else NA_real_,
    min_peak_height = if (present(d$height)) floor(min(d$height)) else NA_real_,
    min_sn = if (present(d$sn)) max(3, min(d$sn)) else NA_real_,
    peakwidth = peakwidth,
    max_mass_shift = if (present(d$mass_shift)) max(d$mass_shift) else NA_real_,
    max_rt_shift = if (present(d$rt_shift)) max(d$rt_shift) else NA_real_
  )
}
