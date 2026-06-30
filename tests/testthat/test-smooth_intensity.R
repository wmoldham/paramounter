# test-smooth_intensity.R

test_that("half_window = 0 returns the input unchanged", {
  x <- c(10, 100, 40, 80, 20)
  expect_identical(smooth_intensity(x, half_window = 0L), x)
  expect_identical(smooth_intensity(x), x) # default
})

test_that("interior points use symmetric triangular weights", {
  # half_window = 1, kernel (1, 2, 1):
  #   y[2] = (1*10 + 2*100 + 1*40) / 4 = 62.5
  #   y[3] = (1*100 + 2*40 + 1*80) / 4 = 65
  out <- smooth_intensity(c(10, 100, 40, 80, 20), half_window = 1L)
  expect_equal(out[2], 62.5)
  expect_equal(out[3], 65)
})

test_that("edges are renormalised over available points", {
  # left edge kernel (2, 1); right edge kernel (1, 2)
  out <- smooth_intensity(c(10, 100, 40, 80, 20), half_window = 1L)
  expect_equal(out[1], (2 * 10 + 1 * 100) / 3)
  expect_equal(out[5], (1 * 80 + 2 * 20) / 3)
})

test_that("output length always equals input length", {
  for (n in 0:5) {
    expect_length(smooth_intensity(runif(30), half_window = n), 30L)
  }
})

test_that("a constant trace is returned unchanged", {
  expect_equal(smooth_intensity(rep(7, 10), half_window = 3L), rep(7, 10))
})

test_that("length-0 and length-1 inputs are returned as-is", {
  expect_identical(smooth_intensity(numeric(0), half_window = 2L), numeric(0))
  expect_identical(smooth_intensity(42, half_window = 2L), 42)
})

test_that("smoothing reduces the variance of a noisy trace", {
  set.seed(42)
  x <- runif(200, 0, 1000)
  expect_lt(var(smooth_intensity(x, half_window = 3L)), var(x))
})

test_that("integer-valued doubles are accepted for half_window", {
  x <- c(10, 100, 40, 80, 20)
  expect_equal(
    smooth_intensity(x, half_window = 1),
    smooth_intensity(x, half_window = 1L)
  )
})

test_that("non-numeric intensity errors", {
  expect_error(smooth_intensity("a", half_window = 1L), "numeric")
  expect_error(smooth_intensity(list(1, 2), half_window = 1L), "numeric")
  expect_error(smooth_intensity(factor(c(1, 2)), half_window = 1L), "numeric")
})

test_that("non-finite intensity errors", {
  expect_error(smooth_intensity(c(1, NA, 3), half_window = 1L), "NA")
  expect_error(smooth_intensity(c(1, NaN, 3), half_window = 1L))
  expect_error(smooth_intensity(c(1, Inf, 3), half_window = 1L))
})

test_that("invalid half_window errors", {
  x <- c(1, 2, 3, 4, 5)
  expect_error(smooth_intensity(x, half_window = -1L), "non-negative")
  expect_error(smooth_intensity(x, half_window = 2.5), "integer")
  expect_error(smooth_intensity(x, half_window = c(1L, 2L)), "single")
  expect_error(smooth_intensity(x, half_window = NA_integer_))
  expect_error(smooth_intensity(x, half_window = "2"))
  expect_error(smooth_intensity(x, half_window = Inf))
})

test_that("matches an independent reference across random inputs", {
  ref_smooth <- function(x, n) {
    if (n == 0L || length(x) <= 1L) return(x)
    kernel <- c(seq_len(n), n + 1L, rev(seq_len(n)))
    out <- numeric(length(x))
    for (i in seq_along(x)) {
      d <- seq(max(-n, 1L - i), min(n, length(x) - i))
      w <- kernel[d + n + 1L]
      out[i] <- sum(w * x[i + d]) / sum(w)
    }
    out
  }
  set.seed(1)
  for (n in 1:4) {
    for (k in 1:50) {
      x <- round(runif(sample(5:80, 1), 0, 1e5), 2)
      expect_equal(smooth_intensity(x, half_window = n), ref_smooth(x, n))
    }
  }
})
