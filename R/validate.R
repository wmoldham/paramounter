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
  # range() settles all three cases in one C-level pass and allocates only its
  # own length-two result: NA and NaN propagate into it, and an infinite value
  # becomes an infinite bound. `any(!is.finite(x))` needs two full-length
  # allocations instead, which showed up as a sixth of the per-file runtime
  # because the measurement loop revalidates the same traces for every mass bin.
  if (length(x) > 0L) {
    bounds <- range(x)
    if (!is.finite(bounds[1L]) || !is.finite(bounds[2L])) {
      stop(
        sprintf("`%s` must not contain NA, NaN, or infinite values.", name),
        call. = FALSE
      )
    }
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
check_unit_fraction <- function(x, name) {
  check_scalar(
    x,
    name,
    "a single number in (0, 1]",
    min = 0,
    min_strict = TRUE,
    max = 1
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

#' Validate that an argument is an instance of an S7 class
#'
#' The class name in the message is taken from the class itself, so the wording
#' cannot drift away from the class it checks.
#'
#' @param x Value to check.
#' @param class An S7 class object, such as [pm_config].
#' @param name Argument name, used in the error message.
#'
#' @return `x`, invisibly.
#' @noRd
check_s7 <- function(x, class, name) {
  if (!S7::S7_inherits(x, class)) {
    stop(
      sprintf("`%s` must be a %s object.", name, class@name),
      call. = FALSE
    )
  }
  invisible(x)
}

#' Validate a list argument
#'
#' @param x Value to check.
#' @param name Argument name, used in the error message.
#'
#' @return `x`, invisibly.
#' @noRd
check_list <- function(x, name) {
  if (!is.list(x)) {
    stop(
      sprintf("`%s` must be a list.", name),
      call. = FALSE
    )
  }
  invisible(x)
}

#' Validate that two parallel lists have element-wise matching lengths
#'
#' The per-scan m/z and intensity lists must line up element by element, not
#' just in overall length. The callers differ in how they check the outer
#' length — [bin_peaks()] compares the two lists, [collect_zoi_masses()]
#' compares both against the EIC — so only this inner check is shared.
#'
#' @param mz_list,int_list Parallel lists to compare.
#' @param mz_name,int_name Argument names, used in the error message.
#'
#' @return `mz_list`, invisibly.
#' @noRd
check_matching_lengths <- function(
    mz_list,
    int_list,
    mz_name = "mz_list",
    int_name = "int_list"
) {
  if (!all(lengths(mz_list) == lengths(int_list))) {
    stop(
      sprintf(
        "Each element of `%s` and `%s` must have matching length.",
        mz_name,
        int_name
      ),
      call. = FALSE
    )
  }
  invisible(mz_list)
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
