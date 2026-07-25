# validate.R

#' Validate a numeric vector argument
#'
#' @param x Value to check.
#' @param name Argument name, used in the error message.
#' @param allow_empty Whether a zero-length vector is acceptable.
#'
#' @return `x`, invisibly. Called for its side effect of raising an error.
#' @noRd
check_numeric_vector <- function(x, name, allow_empty = TRUE) {
  if (!is.numeric(x)) {
    stop(
      sprintf("`%s` must be a numeric vector.", name),
      call. = FALSE
    )
  }
  if (!allow_empty && length(x) == 0L) {
    stop(
      sprintf("`%s` must be a non-empty numeric vector.", name),
      call. = FALSE
    )
  }
  if (any(!is.finite(x))) {
    stop(
      sprintf("`%s` must not contain NA, NaN, or infinite values.", name),
      call. = FALSE
    )
  }
  invisible(x)
}

#' Validate a single numeric value against optional bounds
#'
#' The mechanical core behind the scalar validators. `label` completes the
#' sentence "`name` must be ...", so each caller controls its own wording.
#'
#' @param x Value to check.
#' @param name Argument name, used in the error message.
#' @param label Description of the requirement, e.g. `"a single positive integer"`.
#' @param min,max Inclusive bounds, unless `min_strict` is `TRUE`.
#' @param min_strict Whether `min` is exclusive.
#' @param integer Whether the value must be integer-valued.
#' @param allow_infinite Whether an infinite value is acceptable.
#'
#' @return `x`, invisibly.
#' @noRd
check_scalar <- function(
    x,
    name,
    label,
    min = -Inf,
    min_strict = FALSE,
    max = Inf,
    integer = FALSE,
    allow_infinite = FALSE
) {
  ok <- is.numeric(x) && length(x) == 1L && !is.na(x)
  if (ok && !allow_infinite) {
    ok <- is.finite(x)
  }
  if (ok && integer) {
    # round() rather than as.integer() so that huge values fail cleanly
    # instead of overflowing to NA with a warning
    ok <- is.finite(x) && abs(x) <= .Machine$integer.max && x == round(x)
  }
  if (ok) {
    ok <- if (min_strict) x > min else x >= min
  }
  if (ok) {
    ok <- x <= max
  }
  if (!ok) {
    stop(
      sprintf("`%s` must be %s.", name, label),
      call. = FALSE
    )
  }
  invisible(x)
}

#' @noRd
check_finite_number <- function(x, name) {
  check_scalar(x, name, "a single finite number")
}

#' @noRd
check_nonneg_number <- function(x, name) {
  check_scalar(
    x,
    name,
    "a single non-negative number",
    min = 0
  )
}

#' @noRd
check_positive_number <- function(x, name, allow_infinite = FALSE) {
  label <- if (allow_infinite) {
    "a single positive number (may be Inf)"
  } else {
    "a single positive number"
  }
  check_scalar(
    x,
    name,
    label,
    min = 0,
    min_strict = TRUE,
    allow_infinite = allow_infinite
  )
}

#' @noRd
check_nonneg_integer <- function(x, name) {
  check_scalar(
    x,
    name,
    "a single non-negative integer",
    min = 0,
    integer = TRUE
  )
}

#' @noRd
check_count <- function(x, name) {
  check_scalar(
    x,
    name,
    "a single positive integer",
    min = 1,
    integer = TRUE
  )
}

#' @noRd
check_index <- function(x, name, n, within) {
  check_scalar(
    x,
    name,
    sprintf("a single integer index within `%s`", within),
    min = 1,
    max = n,
    integer = TRUE
  )
}

#' @noRd
check_flag <- function(x, name) {
  if (!is.logical(x) || length(x) != 1L || is.na(x)) {
    stop(
      sprintf("`%s` must be a single TRUE or FALSE.", name),
      call. = FALSE
    )
  }
  invisible(x)
}

#' Convert a stop-based check into an S7 validator message
#'
#' Runs a `check_*` call and returns `NULL` if it passes or the error message as
#' a string if it fails, the form S7 property and class validators expect.
#'
#' @param expr A validation expression that errors on failure.
#'
#' @return `NULL` on success, or the error message string on failure.
#' @noRd
as_message <- function(expr) {
  tryCatch(
    {
      force(expr)
      NULL
    },
    error = function(e) conditionMessage(e)
  )
}
