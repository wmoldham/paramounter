# universal_parameters.R

# The universal quantities measured by a Paramounter analysis, in canonical order.
PM_QUANTITIES <- c(
  "ppm", "mz_diff", "noise", "width_seconds", "width_scans",
  "sn", "height", "mass_shift", "rt_shift"
)

#' An empty set of universal-parameter distributions
#'
#' A named list with one empty numeric vector per quantity in `PM_QUANTITIES`,
#' the template the analysis fills in.
#'
#' @return A named list of empty numeric vectors.
#' @noRd
empty_distributions <- function() {
  stats::setNames(rep(list(numeric(0)), length(PM_QUANTITIES)), PM_QUANTITIES)
}

#' Universal LC-MS parameters measured from the data
#'
#' The result of a Paramounter analysis: the software-agnostic parameter
#' distributions measured across the input files, together with the settings and
#' the files they came from. Software-specific parameter sets (for xcms,
#' MS-DIAL, and MZmine) are produced from this object by the translation
#' functions, and the diagnostic report is drawn from its distributions.
#' Normally created by the analysis rather than by hand.
#'
#' The `summary` property is computed on access from `distributions`, so it
#' always reflects them: one row per quantity giving the count and the minimum,
#' maximum, mean, and median. The translation layer reads whichever statistic a
#' given software parameter requires.
#'
#' @param distributions Named list holding one finite numeric vector per
#'   universal quantity, after outlier trimming. Must contain exactly `ppm` and
#'   `mz_diff` (relative and absolute mass tolerance), `noise`, `width_seconds`
#'   and `width_scans` (peak width), `sn`, `height`, `mass_shift`, and
#'   `rt_shift`. A quantity with no measurements (for example the shifts when
#'   only one file was analysed) is an empty vector.
#' @param files Character vector of the files the parameters were measured from.
#' @param config The [pm_config] used for the analysis.
#'
#' @return A `universal_parameters` object. Its `summary` property is a data
#'   frame with columns `quantity`, `n`, `min`, `max`, `mean`, and `median`
#'
#' @include pm_config.R
#'
#' @examples
#' # normally produced by the analysis; an empty instance shows the structure
#' universal_parameters()@summary
#'
#' @export
universal_parameters <- new_class(
  "universal_parameters",
  properties = list(
    distributions = new_property(
      class_any,
      default = empty_distributions(),
      validator = function(value) {
        if (!is.list(value) || is.null(names(value))) {
          return("`distributions` must be a named list.")
        }
        if (!setequal(names(value), PM_QUANTITIES)) {
          return(sprintf(
            "`distributions` must have exactly these names: %s.",
            paste(PM_QUANTITIES, collapse = ", ")
          ))
        }
        for (nm in names(value)) {
          v <- value[[nm]]
          if (!is.numeric(v) || any(!is.finite(v))) {
            return(sprintf("`distributions$%s` must be a finite numeric vector.", nm))
          }
        }
        NULL
      }
    ),
    files = new_property(
      class_any,
      default = character(0),
      validator = function(value) {
        if (!is.character(value)) "`files` must be a character vector." else NULL
      }
    ),
    config = new_property(pm_config, default = pm_config()),
    summary = new_property(
      class_any,
      getter = function(self) {
        d <- self@distributions[PM_QUANTITIES]
        safe <- function(f) {
          vapply(d, function(v) if (length(v)) f(v) else NA_real_, numeric(1))
        }
        data.frame(
          quantity = PM_QUANTITIES,
          n = vapply(d, length, integer(1)),
          min = safe(min),
          max = safe(max),
          mean = safe(mean),
          median = safe(stats::median),
          row.names = NULL,
          stringsAsFactors = FALSE
        )
      }
    )
  )
)

method(print, universal_parameters) <- function(x, ...) {
  cat("<universal_parameters>\n")
  cat(sprintf(
    "  files: %d   config: <pm_config> legacy = %s\n",
    length(x@files),
    x@config@legacy
  ))
  cat("  parameters (from trimmed distributions):\n")
  s <- x@summary
  s[] <- lapply(s, function(col) if (is.numeric(col)) signif(col, 5) else col)
  print(s, row.names = FALSE)
  invisible(x)
}
