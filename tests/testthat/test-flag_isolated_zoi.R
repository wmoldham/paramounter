# test-flag_isolated_zoi.R

test_that("empty and single-zone bins are handled", {
  expect_identical(flag_isolated_zoi(numeric(0)), logical(0))
  expect_identical(flag_isolated_zoi(100), TRUE)
  expect_identical(flag_isolated_zoi(100, legacy_isolation = FALSE), TRUE)
})

test_that("well-separated zones are isolated and close zones are not", {
  expect_identical(flag_isolated_zoi(c(100, 800)), c(TRUE, TRUE))
  expect_identical(flag_isolated_zoi(c(100, 200)), c(FALSE, FALSE))
  expect_identical(
    flag_isolated_zoi(c(100, 800), legacy_isolation = FALSE),
    c(TRUE, TRUE)
  )
  expect_identical(
    flag_isolated_zoi(c(100, 200), legacy_isolation = FALSE),
    c(FALSE, FALSE)
  )
})

test_that("the separation must be strictly greater than min_gap", {
  expect_identical(flag_isolated_zoi(c(100, 400)), c(FALSE, FALSE))
  expect_identical(flag_isolated_zoi(c(100, 400.001)), c(TRUE, TRUE))
})

test_that("legacy mode tests forward only, corrected mode tests both sides", {
  # zone 2 is 10 s after zone 1 but 690 s before zone 3
  expect_identical(
    flag_isolated_zoi(c(100, 110, 800), legacy_isolation = TRUE),
    c(FALSE, TRUE, TRUE)
  )
  expect_identical(
    flag_isolated_zoi(c(100, 110, 800), legacy_isolation = FALSE),
    c(FALSE, FALSE, TRUE)
  )
})

test_that("a zone isolated on both sides qualifies in either mode", {
  expect_identical(
    flag_isolated_zoi(c(100, 800, 1600), legacy_isolation = FALSE),
    c(TRUE, TRUE, TRUE)
  )
})

test_that("min_gap is honoured", {
  expect_identical(
    flag_isolated_zoi(c(0, 50, 200), min_gap = 60, legacy_isolation = FALSE),
    c(FALSE, FALSE, TRUE)
  )
  expect_identical(
    flag_isolated_zoi(c(100, 101, 102), min_gap = 0, legacy_isolation = FALSE),
    c(TRUE, TRUE, TRUE)
  )
})

test_that("tied retention times are never isolated", {
  expect_identical(flag_isolated_zoi(c(100, 100), min_gap = 0), c(FALSE, FALSE))
  expect_silent(flag_isolated_zoi(c(100, 100, 900)))
})

test_that("legacy mode matches the original selection across random bins", {
  orig_isolated <- function(rt, gap = 300) {
    n <- length(rt)
    keep <- logical(n)
    for (z in seq_len(n)) {
      nxt <- if (z + 1 <= n) rt[z + 1] else NA_real_
      if (z == 1) {
        if (is.na(nxt)) keep[z] <- TRUE
        if (!is.na(nxt)) {
          if (rt[z] < (nxt - gap)) keep[z] <- TRUE
        }
      }
      if (z > 1) {
        if (is.na(nxt) && rt[z] > (rt[z - 1] + gap)) keep[z] <- TRUE
        if (!is.na(nxt)) {
          if (rt[z] < (nxt - gap)) keep[z] <- TRUE
        }
      }
    }
    keep
  }
  set.seed(13)
  for (rep in 1:1000) {
    n <- sample(1:12, 1)
    rt <- sort(cumsum(c(runif(1, 0, 100), runif(max(n - 1, 0), 5, 900))))[seq_len(n)]
    # orig_isolated implements the original's forward-only rule
    expect_identical(flag_isolated_zoi(rt, legacy_isolation = TRUE), orig_isolated(rt))
  }
})

test_that("invalid apex_rt errors", {
  expect_error(flag_isolated_zoi("a"), "numeric")
  expect_error(flag_isolated_zoi(list(1)), "numeric")
  expect_error(flag_isolated_zoi(c(1, NA)), "NA")
  expect_error(flag_isolated_zoi(c(1, Inf)), "NA")
  expect_error(flag_isolated_zoi(c(300, 100)), "non-decreasing")
})

test_that("invalid min_gap and legacy_isolation error", {
  expect_error(flag_isolated_zoi(c(1, 2), min_gap = -1), "non-negative")
  expect_error(flag_isolated_zoi(c(1, 2), min_gap = c(1, 2)), "single")
  expect_error(flag_isolated_zoi(c(1, 2), min_gap = NA_real_), "non-negative")
  expect_error(flag_isolated_zoi(c(1, 2), min_gap = Inf), "non-negative")
  expect_error(flag_isolated_zoi(c(1, 2), legacy_isolation = NA), "TRUE or FALSE")
  expect_error(flag_isolated_zoi(c(1, 2), legacy_isolation = 1), "TRUE or FALSE")
})
