# check_estimates.R

# Which point estimate each measured quantity turns into, which end of the
# distribution it is taken from, and whether it governs how a peak is built or
# only how many candidates pass. The `kind` is what tells a user whether a
# distant estimate is a problem or the intended design.
ESTIMATE_MAP <- list(
  ppm = list(field = "max_ppm", end = "upper", kind = "tolerance"),
  mz_diff = list(field = "max_mz_diff", end = "upper", kind = "tolerance"),
  width_seconds = list(field = "peakwidth", end = "upper", kind = "tolerance"),
  mass_shift = list(field = "max_mass_shift", end = "upper", kind = "tolerance"),
  rt_shift = list(field = "max_rt_shift", end = "upper", kind = "tolerance"),
  noise = list(field = "min_noise", end = "lower", kind = "threshold"),
  width_scans = list(field = "min_peak_scan", end = "lower", kind = "threshold"),
  height = list(field = "min_peak_height", end = "lower", kind = "threshold"),
  sn = list(field = "min_sn", end = "lower", kind = "threshold")
)

#' Check how far each point estimate sits from its distribution
#'
#' Every processing parameter is drawn from one end of a measured distribution.
#' When that distribution is tight the estimate sits close to the bulk and the
#' parameter is representative. When it has a heavy tail the estimate is set by
#' the worst-measured zones in the file rather than by what the instrument can
#' do, and the resulting parameter can be far too permissive.
#'
#' This reports the gap, so you can see it before spending an hour in a peak
#' picker. `ratio` is how many times further from the median the estimate sits:
#' `estimate / median` for parameters taken from the upper end, `median /
#' estimate` for the lower end. A ratio near 1 means the estimate is typical of
#' the data. A large ratio means it is not.
#'
#' Read the `kind` column before acting on a large ratio, because the two kinds
#' fail differently:
#'
#' * A **tolerance** governs how a peak is built. A loose one is not recoverable
#'   later — it widens features and splits peaks that should have merged. On a
#'   high-dynamic-range Orbitrap dataset a `ppm` of 46 against a median of 2
#'   found twice the chromatographic peaks and produced *fewer* features than a
#'   `ppm` of 10. Lower [pm_config]'s `tolerance_quantile` when you see this.
#' * A **threshold** governs how many candidates pass. A loose one only adds
#'   candidates, and downstream filtering removes what does not reproduce. A
#'   large ratio here is the method working as intended, and tightening it has
#'   been measured to destroy real features. Leave these alone unless you have
#'   evidence from your own data.
#'
#' @param params A [universal_parameters] object from [paramounter()].
#'
#' @return A data frame with one row per measured quantity and columns
#'   `quantity`, `kind` (`"tolerance"` or `"threshold"`), `estimate` (the value
#'   the translators will use), `median` (of the measured distribution), and
#'   `ratio`. Quantities with no measurements are omitted.
#'
#' @examples
#' \dontrun{
#' params <- paramounter(files)
#' check_estimates(params)
#' }
#'
#' @export
check_estimates <- function(params) {
  check_s7(params, universal_parameters, "params")
  est <- point_estimates(
    params@distributions,
    params@config@legacy,
    params@config@tolerance_quantile
  )

  measured <- names(ESTIMATE_MAP)[
    vapply(
      names(ESTIMATE_MAP),
      function(q) length(params@distributions[[q]]) > 0L,
      logical(1)
    )
  ]

  rows <- lapply(measured, function(quantity) {
    spec <- ESTIMATE_MAP[[quantity]]
    value <- est[[spec$field]]
    # peakwidth is the only two-element estimate; its upper bound is the one
    # drawn from width_seconds
    if (spec$field == "peakwidth") {
      value <- value[2L]
    }
    middle <- stats::median(params@distributions[[quantity]])
    ratio <- if (spec$end == "upper") value / middle else middle / value
    data.frame(
      quantity = quantity,
      kind = spec$kind,
      estimate = value,
      median = middle,
      ratio = ratio,
      stringsAsFactors = FALSE
    )
  })

  out <- do.call(rbind, rows)
  if (is.null(out)) {
    out <- data.frame(
      quantity = character(0),
      kind = character(0),
      estimate = numeric(0),
      median = numeric(0),
      ratio = numeric(0),
      stringsAsFactors = FALSE
    )
  }
  rownames(out) <- NULL
  out[order(out$kind, -out$ratio), , drop = FALSE]
}
