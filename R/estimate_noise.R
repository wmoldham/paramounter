# estimate_noise.R

#' Estimate the noise floor of an intensity trace
#'
#' Determines a mass-bin-specific noise cutoff by rank-based intensity sorting.
#' The non-zero intensities are sorted ascending and scanned in blocks: at each
#' block boundary `N`, the mean and standard deviation of the lowest `N` values
#' define a threshold `mean + sd_factor * sd`. The scan stops at the first
#' boundary whose next-highest value meets or exceeds that threshold, marking the
#' transition from noise to signal. The cutoff is the intensity at that boundary,
#' returned with the mean and standard deviation of the noise floor it defines
#' (both are needed downstream to compute signal-to-noise ratios).
#'
#' If `block_size` or fewer non-zero values are present, the trace is too short
#' to characterise a floor and the maximum is used as the cutoff. Setting
#' `block_size = 1` evaluates the threshold at every value, the exact definition
#' given in the Paramounter paper; the default of `10` evaluates every tenth
#' value, matching the original implementation and running substantially faster
#' on dense traces.
#'
#' Only positive intensities are considered; zeros and any non-positive values
#' are ignored.
#'
#' @param intensity Numeric vector of intensities for a single mass bin,
#'   typically the smoothed extracted-ion chromatogram from
#'   [smooth_intensity()]. Must be finite: no `NA`, `NaN`, or infinite values.
#' @param block_size Single positive integer giving the block width for the
#'   rank-based scan (default `10`).
#' @param sd_factor Single non-negative number multiplying the standard
#'   deviation in the noise threshold (default `3`).
#'
#' @return A named list with elements `cutoff` (the estimated noise cutoff
#'   intensity), `noise_mean` (mean of the noise floor at the cutoff), and
#'   `noise_sd` (standard deviation of the noise floor at the cutoff). A trace
#'   with no positive values returns zeros for all three.
#'
#' @examples
#' set.seed(1)
#' eic <- c(rep(0, 5), rlnorm(40, 6, 0.5), 1e5, 8e4)
#' estimate_noise(eic)
#'
#' @export
estimate_noise <- function(intensity, block_size = 10L, sd_factor = 3) {
  if (!is.numeric(intensity)) {
    stop(
      "`intensity` must be a numeric vector.",
      call. = FALSE
    )
  }
  if (any(!is.finite(intensity))) {
    stop(
      "`intensity` must not contain NA, NaN, or infinite values.",
      call. = FALSE
    )
  }
  if (!is.numeric(block_size) ||
      length(block_size) != 1L ||
      !is.finite(block_size) ||
      block_size < 1 ||
      block_size != as.integer(block_size)) {
    stop(
      "`block_size` must be a single positive integer.",
      call. = FALSE
    )
  }
  if (!is.numeric(sd_factor) ||
      length(sd_factor) != 1L ||
      !is.finite(sd_factor) ||
      sd_factor < 0) {
    stop(
      "`sd_factor` must be a single non-negative number.",
      call. = FALSE
    )
  }
  block_size <- as.integer(block_size)

  s <- sort(intensity[intensity > 0])
  n_signal <- length(s)
  if (n_signal == 0L) {
    return(list(cutoff = 0, noise_mean = 0, noise_sd = 0))
  }
  if (n_signal <= block_size) {
    return(list(
      cutoff = s[n_signal],
      noise_mean = mean(s),
      noise_sd = if (n_signal >= 2L) stats::sd(s) else 0
    ))
  }

  checkpoints <- seq.int(block_size, n_signal, by = block_size)
  cumulative_sum <- cumsum(s)[checkpoints]
  cumulative_sq <- cumsum(s^2)[checkpoints]
  means <- cumulative_sum / checkpoints
  variances <- (cumulative_sq - cumulative_sum^2 / checkpoints) / (checkpoints - 1)
  sds <- sqrt(pmax(variances, 0))
  thresholds <- means + sd_factor * sds

  has_next <- checkpoints < n_signal
  next_value <- ifelse(has_next, s[checkpoints + 1L], NA_real_)
  exceeded <- has_next & (next_value >= thresholds)

  first_break <- which(exceeded)[1]
  position <- if (is.na(first_break)) length(checkpoints) else first_break
  cutoff_index <- checkpoints[position]

  list(
    cutoff = s[cutoff_index],
    noise_mean = means[position],
    noise_sd = sds[position]
  )
}
