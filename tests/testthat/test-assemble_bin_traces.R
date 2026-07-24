# test-assemble_bin_traces.R

test_that("a bin expands to per-scan traces with a mean EIC", {
  tr <- assemble_bin_traces(
    list(scan = c(1L, 2L), mz = c(100.01, 100.02), intensity = c(10, 40)),
    3L
  )
  expect_identical(tr$mz_list, list(100.01, 100.02, numeric(0)))
  expect_identical(tr$int_list, list(10, 40, numeric(0)))
  expect_equal(tr$eic, c(10, 40, 0))
})

test_that("several peaks in one scan are grouped and averaged", {
  tr <- assemble_bin_traces(
    list(scan = c(1L, 1L, 2L), mz = c(100.01, 100.03, 100.02), intensity = c(10, 20, 40)),
    3L
  )
  expect_identical(tr$mz_list[[1]], c(100.01, 100.03))
  expect_equal(tr$eic, c(15, 40, 0))
})

test_that("an empty bin gives empty traces and a zero EIC", {
  tr <- assemble_bin_traces(list(scan = integer(0), mz = numeric(0), intensity = numeric(0)), 2L)
  expect_identical(tr$mz_list, list(numeric(0), numeric(0)))
  expect_equal(tr$eic, c(0, 0))
})

test_that("the bin structure is validated", {
  expect_error(assemble_bin_traces(42, 3), "list with")
  expect_error(assemble_bin_traces(list(scan = 1, mz = 1), 3), "list with")
  expect_error(assemble_bin_traces(list(scan = "a", mz = 1, intensity = 1), 3), "must be numeric")
  expect_error(assemble_bin_traces(list(scan = c(1, 2), mz = 1, intensity = 1), 3), "same length")
})

test_that("scan indices and n_scans are validated", {
  expect_error(assemble_bin_traces(list(scan = 5, mz = 1, intensity = 1), 3), "within \\[1, n_scans\\]")
  expect_error(assemble_bin_traces(list(scan = 0, mz = 1, intensity = 1), 3), "within \\[1, n_scans\\]")
  expect_error(assemble_bin_traces(list(scan = 2.5, mz = 1, intensity = 1), 3), "within \\[1, n_scans\\]")
  expect_error(assemble_bin_traces(list(scan = 1, mz = NA_real_, intensity = 1), 3), "finite")
  expect_error(assemble_bin_traces(list(scan = 1, mz = 1, intensity = 1), 0), "positive integer")
})
