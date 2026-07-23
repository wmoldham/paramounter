# test-validate.R

test_that("check_numeric_vector accepts valid vectors", {
  expect_silent(check_numeric_vector(c(1, 2, 3), "x"))
  expect_silent(check_numeric_vector(numeric(0), "x"))
  expect_silent(check_numeric_vector(1L, "x"))
})

test_that("check_numeric_vector rejects non-numeric and non-finite input", {
  expect_error(check_numeric_vector("a", "x"), "numeric")
  expect_error(check_numeric_vector(list(1), "x"), "numeric")
  expect_error(check_numeric_vector(factor(1), "x"), "numeric")
  expect_error(check_numeric_vector(TRUE, "x"), "numeric")
  expect_error(check_numeric_vector(c(1, NA), "x"), "NA")
  expect_error(check_numeric_vector(c(1, NaN), "x"), "NA")
  expect_error(check_numeric_vector(c(1, Inf), "x"), "NA")
})

test_that("check_numeric_vector honours allow_empty", {
  expect_error(check_numeric_vector(numeric(0), "x", allow_empty = FALSE), "non-empty")
})

test_that("the error message names the offending argument", {
  expect_error(check_numeric_vector("a", "rtime"), "rtime")
  expect_error(check_nonneg_number(-1, "cutoff"), "cutoff")
})

test_that("check_nonneg_number bounds at zero", {
  expect_silent(check_nonneg_number(0, "x"))
  expect_silent(check_nonneg_number(3.5, "x"))
  expect_error(check_nonneg_number(-1, "x"), "non-negative")
  expect_error(check_nonneg_number(c(1, 2), "x"), "single")
  expect_error(check_nonneg_number(NA_real_, "x"), "non-negative")
  expect_error(check_nonneg_number("1", "x"), "single")
  expect_error(check_nonneg_number(Inf, "x"), "non-negative")
  expect_error(check_nonneg_number(TRUE, "x"), "single")
})

test_that("check_positive_number excludes zero", {
  expect_silent(check_positive_number(0.001, "x"))
  expect_error(check_positive_number(0, "x"), "positive")
  expect_error(check_positive_number(-1, "x"), "positive")
  expect_error(check_positive_number(Inf, "x"), "positive")
})

test_that("check_positive_number honours allow_infinite", {
  expect_silent(check_positive_number(Inf, "x", allow_infinite = TRUE))
  expect_error(check_positive_number(NA_real_, "x", allow_infinite = TRUE), "positive")
  expect_error(check_positive_number(0, "x", allow_infinite = TRUE), "positive")
  expect_error(check_positive_number(-Inf, "x", allow_infinite = TRUE), "positive")
})

test_that("check_nonneg_integer requires whole numbers from zero up", {
  expect_silent(check_nonneg_integer(0, "x"))
  expect_silent(check_nonneg_integer(3, "x"))
  expect_silent(check_nonneg_integer(3L, "x"))
  expect_error(check_nonneg_integer(-1, "x"), "non-negative")
  expect_error(check_nonneg_integer(2.5, "x"), "integer")
  expect_error(check_nonneg_integer(NA_integer_, "x"), "integer")
  expect_error(check_nonneg_integer(Inf, "x"), "integer")
  expect_error(check_nonneg_integer(c(1, 2), "x"), "single")
})

test_that("check_count requires whole numbers from one up", {
  expect_silent(check_count(1, "x"))
  expect_silent(check_count(10L, "x"))
  expect_error(check_count(0, "x"), "positive")
  expect_error(check_count(-1, "x"), "positive")
  expect_error(check_count(2.5, "x"), "integer")
  expect_error(check_count(Inf, "x"), "integer")
  expect_error(check_count("10", "x"), "single")
})

test_that("integer checks reject values that would overflow, without warning", {
  expect_error(check_count(3e9, "x"), "integer")
  expect_warning(try(check_count(3e9, "x"), silent = TRUE), NA)
})

test_that("check_index bounds against the container length", {
  expect_silent(check_index(1, "i", 5, "rtime"))
  expect_silent(check_index(5, "i", 5, "rtime"))
  expect_error(check_index(0, "i", 5, "rtime"), "within")
  expect_error(check_index(6, "i", 5, "rtime"), "within")
  expect_error(check_index(2.5, "i", 5, "rtime"), "integer")
  expect_error(check_index(6, "i", 5, "rtime"), "rtime")
})

test_that("check_flag requires a single non-missing logical", {
  expect_silent(check_flag(TRUE, "x"))
  expect_silent(check_flag(FALSE, "x"))
  expect_error(check_flag(NA, "x"), "TRUE or FALSE")
  expect_error(check_flag(1, "x"), "TRUE or FALSE")
  expect_error(check_flag(c(TRUE, FALSE), "x"), "TRUE or FALSE")
  expect_error(check_flag("TRUE", "x"), "TRUE or FALSE")
})
