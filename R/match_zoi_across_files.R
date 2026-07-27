# match_zoi_across_files.R

#' Match clean zones of interest across files
#'
#' Aligns the isolated ("clean") zones of interest from several files so that
#' the same feature can be tracked between samples, which is the basis for
#' estimating instrument mass and retention-time shift. Matching is anchored on
#' the first file: for every subsequent file, each anchor zone is paired with a
#' zone falling inside a mass and retention-time window centred on the anchor's
#' own m/z and retention time. A feature is kept only if it can be matched in
#' every file.
#'
#' The window is always centred on the anchor file, never on a running consensus
#' or the previous file, so all files are compared against the same reference.
#' Supply the zones in the order they were detected. Row order changes the result
#' only when `legacy_first_match` is `TRUE`.
#'
#' A feature that is absent from any one file is dropped, so pass only files that
#' actually contain clean zones; an empty file causes every feature to be
#' discarded.
#'
#' @param zoi_list List of at least two data frames, one per file, each with
#'   numeric, finite `mz` and `rt` columns giving the reference m/z and apex
#'   retention time of that file's clean zones of interest, meaning the ones
#'   [flag_isolated_zoi()] marked. The first element is the anchor. Elements may
#'   be named to label the output columns.
#' @param mz_tol Single positive m/z half-window, in the same units as `mz`
#'   (default `0.015`). A candidate matches when its m/z is within `mz_tol` of
#'   the anchor's, inclusive.
#' @param rt_tol Single positive retention-time half-window, in the same units
#'   as `rt` (default `30`). Applied inclusively, like `mz_tol`.
#' @param legacy_first_match Single `TRUE` or `FALSE`. When `FALSE` (the
#'   default) the *nearest* candidate is taken, minimising the window-normalised
#'   distance `sqrt((dmz / mz_tol)^2 + (drt / rt_tol)^2)`, which is the
#'   principled choice when a window contains several candidates. When `TRUE`
#'   the original Paramounter behaviour is reproduced: the *first* candidate
#'   within the window (by row order) is taken, which can pair an anchor with a
#'   far candidate while a nearer one sits in the same window; use that only to
#'   compare against the published values. Keep this in step with [pm_config]'s
#'   `legacy`, which is what sets it when the pipeline is run through
#'   [paramounter()].
#'
#' @return A named list with two numeric matrices, `mz` and `rt`, each with one
#'   row per matched feature and one column per file (in the order of
#'   `zoi_list`, with column names taken from its names when present). Column one
#'   holds the anchor values. When nothing matches across all files, both
#'   matrices have zero rows.
#'
#' @examples
#' file1 <- data.frame(mz = c(200.000, 300.000), rt = c(500, 800))
#' file2 <- data.frame(mz = c(200.003, 300.002), rt = c(505, 803))
#' match_zoi_across_files(list(file1, file2))
#'
#' @export
match_zoi_across_files <- function(
    zoi_list,
    mz_tol = 0.015,
    rt_tol = 30,
    legacy_first_match = FALSE
) {
  if (!is.list(zoi_list) || is.data.frame(zoi_list)) {
    stop(
      "`zoi_list` must be a list of data frames.",
      call. = FALSE
    )
  }
  if (length(zoi_list) < 2L) {
    stop(
      "`zoi_list` must contain at least two files.",
      call. = FALSE
    )
  }
  for (i in seq_along(zoi_list)) {
    element <- zoi_list[[i]]
    if (!is.data.frame(element) || !all(c("mz", "rt") %in% names(element))) {
      stop(
        sprintf("`zoi_list[[%d]]` must be a data frame with `mz` and `rt` columns.", i),
        call. = FALSE
      )
    }
    check_numeric_vector(element$mz, sprintf("zoi_list[[%d]]$mz", i))
    check_numeric_vector(element$rt, sprintf("zoi_list[[%d]]$rt", i))
  }
  check_positive_number(mz_tol, "mz_tol")
  check_positive_number(rt_tol, "rt_tol")
  check_flag(legacy_first_match, "legacy_first_match")

  n_files <- length(zoi_list)
  anchor <- zoi_list[[1]]
  n_anchor <- nrow(anchor)
  matched_mz <- matrix(NA_real_, nrow = n_anchor, ncol = n_files)
  matched_rt <- matrix(NA_real_, nrow = n_anchor, ncol = n_files)
  matched_mz[, 1] <- anchor$mz
  matched_rt[, 1] <- anchor$rt
  alive <- rep(TRUE, n_anchor)

  if (n_anchor > 0L) {
    for (k in 2:n_files) {
      fk <- zoi_list[[k]]
      for (i in which(alive)) {
        anchor_mz <- anchor$mz[i]
        anchor_rt <- anchor$rt[i]
        in_window <- which(
          abs(fk$mz - anchor_mz) <= mz_tol &
            abs(fk$rt - anchor_rt) <= rt_tol
        )
        if (length(in_window) == 0L) {
          alive[i] <- FALSE
        } else {
          if (legacy_first_match) {
            pick <- in_window[1]
          } else {
            distance <- sqrt(
              (abs(fk$mz[in_window] - anchor_mz) / mz_tol)^2 +
                (abs(fk$rt[in_window] - anchor_rt) / rt_tol)^2
            )
            pick <- in_window[which.min(distance)]
          }
          matched_mz[i, k] <- fk$mz[pick]
          matched_rt[i, k] <- fk$rt[pick]
        }
      }
    }
  }

  file_names <- names(zoi_list)
  matched_mz <- matched_mz[alive, , drop = FALSE]
  matched_rt <- matched_rt[alive, , drop = FALSE]
  colnames(matched_mz) <- file_names
  colnames(matched_rt) <- file_names
  list(mz = matched_mz, rt = matched_rt)
}
