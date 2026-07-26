# aggregate_files.R

#' Trim outliers from one tail of a distribution
#'
#' Keeps the given fraction of the values, dropping the outliers from whichever
#' tail is not wanted: `remove = "high"` keeps the lowest `keep` fraction (for
#' quantities summarised by their maximum), `remove = "low"` keeps the highest
#' `keep` fraction (for quantities summarised by their minimum). Non-finite
#' values are dropped first. The kept count is `round(n * keep)`, matching the
#' original.
#'
#' @param x Numeric vector.
#' @param keep Fraction to keep in `(0, 1]`.
#' @param remove Which tail to drop, `"high"` or `"low"`.
#'
#' @return The kept values as a numeric vector.
#' @noRd
trim_distribution <- function(x, keep = 0.97, remove = c("high", "low")) {
  remove <- match.arg(remove)
  x <- x[is.finite(x)]
  n <- length(x)
  if (n == 0L) {
    return(numeric(0))
  }
  k <- round(n * keep)
  if (k < 1L) {
    return(numeric(0))
  }
  sort(x, decreasing = (remove == "low"))[seq_len(k)]
}

#' Determine the ppm cutoff from the pooled ppm distribution
#'
#' Reproduces the original's cutoff: sort the pooled ppm values, keep the lowest
#' `trim` fraction, and take the value at the `quantile` position of that trimmed
#' set. Returns `Inf` (no filtering) when there are no values.
#'
#' @param ppm Numeric vector of pooled ppm values.
#' @param trim Fraction kept when trimming outliers before the quantile.
#' @param quantile Quantile position within the trimmed set.
#'
#' @return A single cutoff value.
#' @noRd
compute_ppm_cutoff <- function(ppm, trim, quantile) {
  s <- sort(ppm[is.finite(ppm)])
  n <- length(s)
  if (n == 0L) {
    return(Inf)
  }
  trimmed <- s[seq_len(max(round(n * trim), 1L))]
  trimmed[max(round(length(trimmed) * quantile), 1L)]
}

#' Aggregate per-file measurements into universal parameters
#'
#' Combines the per-file outputs of the measurement step into the final
#' [universal_parameters]. It pools the per-bin noise and stacks the per-zone
#' tables across files; determines the relative-tolerance cutoff (the manual
#' `config@ppm_cutoff`, or automatically from the pooled ppm at
#' `config@ppm_quantile`); drops zones failing that cutoff; trims each quantity's
#' distribution (dropping the tail opposite its summary statistic, per
#' `config@trim`); and, when more than one file is given, derives each file's
#' clean zones (surviving and isolated) and runs [match_zoi_across_files()] then
#' [estimate_instrument_shift()] to obtain the mass- and retention-time-shift
#' distributions. Called by the orchestrator on the results of the per-file
#' measurement.
#'
#' @param per_file List with one element per file, each the list returned by the
#'   per-file measurement step (with numeric `noise` and a data frame `zoi`).
#' @param files Character vector of the file paths, one per element of
#'   `per_file`.
#' @param config The [pm_config] giving the analysis settings.
#'
#' @return A [universal_parameters] object.
#' @noRd
aggregate_files <- function(per_file, files, config = pm_config()) {
  if (!is.list(per_file) || length(per_file) == 0L) {
    stop("`per_file` must be a non-empty list.", call. = FALSE)
  }
  for (i in seq_along(per_file)) {
    pf <- per_file[[i]]
    if (!is.list(pf) ||
        !all(c("noise", "zoi") %in% names(pf)) ||
        !is.numeric(pf$noise) ||
        !is.data.frame(pf$zoi)) {
      stop(
        sprintf("`per_file[[%d]]` must be a list with numeric `noise` and a data frame `zoi`.", i),
        call. = FALSE
      )
    }
  }
  if (!is.character(files) || length(files) != length(per_file)) {
    stop("`files` must be a character vector, one per file.", call. = FALSE)
  }
  if (!S7::S7_inherits(config, pm_config)) {
    stop("`config` must be a pm_config object.", call. = FALSE)
  }

  noise_pooled <- unlist(lapply(per_file, `[[`, "noise"))
  zoi_pooled <- do.call(rbind, lapply(per_file, `[[`, "zoi"))

  cutoff <- if (is.null(config@ppm_cutoff)) {
    compute_ppm_cutoff(zoi_pooled$ppm, config@trim, config@ppm_quantile)
  } else {
    config@ppm_cutoff
  }
  kept <- zoi_pooled[zoi_pooled$ppm < cutoff, , drop = FALSE]
  trim <- config@trim

  distributions <- empty_distributions()
  distributions$ppm <- as.numeric(kept$ppm)
  distributions$mz_diff <- trim_distribution(kept$mz_diff, trim, "high")
  distributions$noise <- trim_distribution(noise_pooled, trim, "high")
  distributions$width_seconds <- trim_distribution(kept$width_seconds, trim, "high")
  distributions$width_scans <- as.numeric(trim_distribution(kept$width_scans, trim, "high"))
  distributions$sn <- trim_distribution(kept$sn, trim, "low")
  distributions$height <- trim_distribution(kept$height, trim, "low")

  if (length(per_file) > 1L) {
    clean <- lapply(per_file, function(pf) {
      z <- pf$zoi
      keep_clean <- z$ppm < cutoff & z$isolated
      data.frame(mz = z$reference_mz[keep_clean], rt = z$apex_rt[keep_clean])
    })
    matched <- match_zoi_across_files(
      clean,
      config@match_mz_tol,
      config@match_rt_tol,
      config@legacy
    )
    shifts <- estimate_instrument_shift(matched)
    distributions$mass_shift <- trim_distribution(shifts$mass_shift, trim, "high")
    distributions$rt_shift <- trim_distribution(shifts$rt_shift, trim, "high")
  }

  universal_parameters(distributions = distributions, files = files, config = config)
}
