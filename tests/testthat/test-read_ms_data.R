# test-read_ms_data.R

test_that("the file path is validated", {
  expect_error(read_ms_data(42), "single file path")
  expect_error(read_ms_data(c("a", "b")), "single file path")
  expect_error(read_ms_data(NA_character_), "single file path")
  expect_error(read_ms_data(character(0)), "single file path")
  expect_error(read_ms_data("/no/such/file.mzML"), "File not found")
})

test_that("ms_level is validated before any file is read", {
  tmp <- tempfile(fileext = ".mzML")
  file.create(tmp)
  expect_error(read_ms_data(tmp, ms_level = 0), "positive integer")
  expect_error(read_ms_data(tmp, ms_level = 2.5), "positive integer")
  expect_error(read_ms_data(tmp, ms_level = c(1, 2)), "single")
  expect_error(read_ms_data(tmp, ms_level = NA_integer_), "positive integer")
})

test_that("a real file is read into parallel per-scan traces", {
  skip_if_not_installed("Spectra")
  skip_if_not_installed("msdata")
  f <- system.file("microtofq/MM14.mzML", package = "msdata")
  skip_if(f == "", "msdata example file not found")
  d <- read_ms_data(f)
  expect_named(d, c("mz", "intensity", "rtime"))
  expect_type(d$mz, "list")
  expect_type(d$intensity, "list")
  expect_type(d$rtime, "double")
  expect_equal(length(d$mz), length(d$rtime))
  expect_equal(length(d$intensity), length(d$rtime))
  expect_true(all(lengths(d$mz) == lengths(d$intensity)))
})
