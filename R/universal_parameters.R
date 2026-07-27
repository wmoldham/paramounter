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
#' functions, and `plot()` draws its distributions. Normally created by the
#' analysis rather than by hand.
#'
#' The `summary` property is computed on access from `distributions`, so it
#' always reflects them: one row per quantity giving the count and the minimum,
#' maximum, mean, and median. The translation layer reads whichever statistic a
#' given software parameter requires.
#'
#' @param distributions Named list holding one finite numeric vector per
#'   universal quantity. Must contain exactly `ppm` and `mz_diff` (relative mass
#'   tolerance in ppm and absolute in Da), `noise`, `width_seconds` and
#'   `width_scans` (peak width in seconds and in scans), `sn`, `height`,
#'   `mass_shift` (Da), and `rt_shift` (seconds). A quantity with no
#'   measurements, such as the shifts when only one file was analysed, is an
#'   empty vector.
#' @param files Character vector of the files the parameters were measured from.
#' @param config The [pm_config] used for the analysis.
#'
#' @return A `universal_parameters` object. Its `summary` property is a data
#'   frame with columns `quantity`, `n`, `min`, `max`, `mean`, and `median`.
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
      # quoted for the same reason as `config` below, and additionally because
      # an evaluated list here is a *value* while the generated \usage holds a
      # `list(...)` *call*: they deparse identically but do not compare equal,
      # which R CMD check reports as a codoc mismatch it cannot render.
      default = quote(empty_distributions()),
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
    # quote() so the default is evaluated when an object is built rather than
    # when the class is defined. Constructing it here instead left roxygen
    # deparsing an already-built S7 object into `config = <object>`, which is
    # not valid R and so produced a malformed \usage section.
    config = new_property(pm_config, default = quote(pm_config())),
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

#' Plot the universal-parameter distributions
#'
#' Draws a diagnostic panel of histograms, one per measured quantity, with the
#' driving statistic marked — the minimum or maximum that feeds the software
#' parameters (both, for peak width in seconds). Quantities with no measurements,
#' such as the instrument shifts when only one file was analysed, are omitted.
#' The exact per-software parameter values come from [to_xcms()], [to_msdial()],
#' and [to_mzmine()].
#'
#' @param x A [universal_parameters] object.
#' @param y Ignored, and present only because the `plot` generic it specialises
#'   takes it. The panel is drawn entirely from `x`.
#' @param ... Ignored.
#'
#' @return `x`, invisibly.
#'
# No \usage block. roxygen2 renders an S7 method on an S4 generic as a bare
# `plot(x, y, ...)`, which R CMD check reads as documenting `plot` itself and
# wants an \alias{plot} for; the \S4method{} form it would accept cannot express
# the `paramounter::universal_parameters` signature, because the `::` breaks its
# \usage parser. The call is just `plot(params)`, shown in the examples.
#' @usage NULL
#'
#' @examples
#' \dontrun{
#' params <- paramounter(files)
#' plot(params)
#' }
#'
#' @export
method(plot, universal_parameters) <- function(x, y, ...) {
  d <- x@distributions
  meta <- list(
    ppm = list(label = "mass tolerance (ppm)", drive = "max"),
    mz_diff = list(label = "mass tolerance (Da)", drive = "max"),
    noise = list(label = "noise level", drive = "min"),
    width_seconds = list(label = "peak width (s)", drive = "range"),
    width_scans = list(label = "peak width (scans)", drive = "min"),
    sn = list(label = "signal / noise", drive = "min"),
    height = list(label = "peak height (log10)", drive = "min", log = TRUE),
    mass_shift = list(label = "mass shift (Da)", drive = "max"),
    rt_shift = list(label = "RT shift (s)", drive = "max")
  )
  canonical_order <- names(meta)
  present <- canonical_order[vapply(d[canonical_order], length, integer(1)) > 0L]
  if (length(present) == 0L) {
    stop("No non-empty distributions to plot.", call. = FALSE)
  }

  n <- length(present)
  ncols <- if (n == 1L) 1L else if (n <= 4L) 2L else 3L
  nrows <- ceiling(n / ncols)
  op <- graphics::par(mfrow = c(nrows, ncols), mar = c(4, 4, 3, 1))
  on.exit(graphics::par(op))

  for (q in present) {
    m <- meta[[q]]
    v <- d[[q]]
    plot_values <- if (isTRUE(m$log)) log10(v) else v
    graphics::hist(
      plot_values,
      main = m$label,
      xlab = "",
      ylab = "count",
      col = "grey85",
      border = "white"
    )
    mark <- function(value, text) {
      position <- if (isTRUE(m$log)) log10(value) else value
      graphics::abline(v = position, col = "firebrick", lty = 2, lwd = 2)
      graphics::mtext(text, side = 3, at = position, cex = 0.7, col = "firebrick")
    }
    if (m$drive == "max") {
      mark(max(v), sprintf("max %.4g", max(v)))
    } else if (m$drive == "min") {
      mark(min(v), sprintf("min %.4g", min(v)))
    } else {
      mark(min(v), sprintf("min %.4g", min(v)))
      mark(max(v), sprintf("max %.4g", max(v)))
    }
  }
  invisible(x)
}
