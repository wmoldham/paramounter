# test-bin_peaks.R

test_that("peaks are grouped into their half-open bins in scan order", {
  b <- bin_peaks(
    list(c(100.01, 100.06, 100.12), c(100.02, 100.07), numeric(0)),
    list(c(10, 20, 30), c(40, 50), numeric(0)),
    bin_width = 0.05,
    mz_range = c(100.00, 100.15)
  )
  expect_equal(b$n_scans, 3L)
  expect_length(b$bins, 3L)
  expect_identical(b$bins[[1]]$scan, c(1L, 2L))
  expect_equal(b$bins[[1]]$mz, c(100.01, 100.02))
  expect_equal(b$bins[[1]]$intensity, c(10, 40))
  expect_equal(b$bins[[1]]$mz_low, 100.00)
  expect_equal(b$bins[[1]]$mz_high, 100.05)
  expect_equal(b$bins[[2]]$mz, c(100.06, 100.07))
  expect_equal(b$bins[[3]]$mz, 100.12)
  expect_identical(b$bins[[3]]$scan, 1L)
})

test_that("empty bins are omitted", {
  b <- bin_peaks(list(c(100.01, 100.02)), list(c(1, 2)), 0.05, c(100.00, 100.20))
  expect_length(b$bins, 1L)
  expect_equal(b$bins[[1]]$bin, 1L)
})

test_that("scans with no peaks give no bins", {
  b <- bin_peaks(list(numeric(0), numeric(0)), list(numeric(0), numeric(0)))
  expect_length(b$bins, 0L)
  expect_equal(b$n_scans, 2L)
})

test_that("the highest peak on the top edge is excluded, as in the original", {
  b <- bin_peaks(list(c(100.00, 100.15)), list(c(5, 9)))
  expect_length(b$bins, 1L)
  expect_equal(b$bins[[1]]$mz, 100.00)
})

test_that("a supplied mz_range drops out-of-range peaks", {
  b <- bin_peaks(list(c(50, 100.02, 500)), list(c(1, 2, 3)), 0.05, c(100.00, 100.10))
  expect_length(b$bins, 1L)
  expect_equal(b$bins[[1]]$mz, 100.02)
})

test_that("the list structure and values are validated", {
  expect_error(bin_peaks(1, list(1)), "`mz_list` must be a list")
  expect_error(bin_peaks(list(1), 1), "`int_list` must be a list")
  expect_error(bin_peaks(list(1, 2), list(1)), "same length")
  expect_error(bin_peaks(list(c(1, 2)), list(1)), "matching length")
  expect_error(bin_peaks(list("a"), list(1)), "finite numeric")
  expect_error(bin_peaks(list(c(1, NA)), list(c(1, 2))), "finite numeric")
})

test_that("bin_width and mz_range are validated", {
  expect_error(bin_peaks(list(1), list(1), bin_width = 0), "positive")
  expect_error(bin_peaks(list(1), list(1), bin_width = -1), "positive")
  expect_error(bin_peaks(list(1), list(1), mz_range = 100), "length-2")
  expect_error(bin_peaks(list(1), list(1), mz_range = c(100, 100)), "mz_range\\[1\\]")
  expect_error(bin_peaks(list(1), list(1), mz_range = c(200, 100)), "mz_range\\[1\\]")
  expect_error(bin_peaks(list(1), list(1), mz_range = c(100, NA)), "length-2")
})
