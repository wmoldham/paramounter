# test-estimate_noise.R

test_that("a trace with no positive values returns zeros", {
  expect_identical(
    estimate_noise(rep(0, 8)),
    list(cutoff = 0, noise_mean = 0, noise_sd = 0)
  )
})

test_that("fallback path uses the maximum when too few non-zero values", {
  # 4 non-zero values, default block_size 10
  r <- estimate_noise(c(0, 50, 30, 70, 40))
  expect_equal(r$cutoff, 70)
  expect_equal(r$noise_mean, 47.5)
  expect_equal(r$noise_sd, sd(c(30, 40, 50, 70)))
})

test_that("a single positive value gives a zero standard deviation", {
  expect_identical(
    estimate_noise(c(0, 0, 42)),
    list(cutoff = 42, noise_mean = 42, noise_sd = 0)
  )
})

test_that("non-positive values are ignored", {
  expect_equal(
    estimate_noise(c(-5, 0, 50, 30, 70, 40)),
    estimate_noise(c(50, 30, 70, 40))
  )
})

test_that("iterative path stops at the noise floor and excludes a spike", {
  # tight floor at 100 (11 values) then a spike: cutoff sits at the floor
  r <- estimate_noise(c(rep(100, 11), 99999))
  expect_equal(r$cutoff, 100)
  expect_equal(r$noise_mean, 100)
  expect_equal(r$noise_sd, 0)
})

test_that("integer-valued doubles are accepted for block_size", {
  set.seed(8)
  x <- rlnorm(50, 6, 0.5)
  expect_equal(
    estimate_noise(x, block_size = 10),
    estimate_noise(x, block_size = 10L)
  )
})

test_that("matches an independent oracle across random inputs and block sizes", {
  ref_noise <- function(x, block = 10L, sdf = 3) {
    s <- sort(x[x > 0])
    n <- length(s)
    if (n == 0L) return(list(cutoff = 0, noise_mean = 0, noise_sd = 0))
    if (n <= block) {
      return(list(
        cutoff = s[n],
        noise_mean = mean(s),
        noise_sd = if (n >= 2L) sd(s) else 0
      ))
    }
    cps <- seq.int(block, n, by = block)
    pick <- cps[length(cps)]
    for (N in cps) {
      if (N >= 2 && N < n && s[N + 1] >= mean(s[1:N]) + sdf * sd(s[1:N])) {
        pick <- N
        break
      }
    }
    list(cutoff = s[pick], noise_mean = mean(s[1:pick]), noise_sd = sd(s[1:pick]))
  }
  make_eic <- function(L) {
    x <- numeric(L)
    nz <- sample(L, max(1, rbinom(1, L, 0.6)))
    x[nz] <- rlnorm(length(nz), 6, 0.8)
    sp <- sample(L, sample(0:6, 1))
    x[sp] <- x[sp] + rlnorm(length(sp), 11, 0.5)
    round(x, 2)
  }
  set.seed(3)
  for (b in c(1, 5, 10, 25)) {
    for (k in 1:100) {
      x <- make_eic(sample(5:1000, 1))
      expect_equal(estimate_noise(x, block_size = b), ref_noise(x, block = b))
    }
  }
})

test_that("non-numeric intensity errors", {
  expect_error(estimate_noise("a"), "numeric")
  expect_error(estimate_noise(list(1, 2)), "numeric")
  expect_error(estimate_noise(factor(c(1, 2))), "numeric")
})

test_that("non-finite intensity errors", {
  expect_error(estimate_noise(c(1, NA, 3)), "NA")
  expect_error(estimate_noise(c(1, NaN, 3)))
  expect_error(estimate_noise(c(1, Inf, 3)))
})

test_that("invalid block_size errors", {
  x <- c(1, 2, 3)
  expect_error(estimate_noise(x, block_size = 0), "positive")
  expect_error(estimate_noise(x, block_size = -1), "positive")
  expect_error(estimate_noise(x, block_size = 2.5), "integer")
  expect_error(estimate_noise(x, block_size = c(1L, 2L)), "single")
  expect_error(estimate_noise(x, block_size = NA_integer_))
  expect_error(estimate_noise(x, block_size = "10"))
  expect_error(estimate_noise(x, block_size = Inf))
})

test_that("invalid sd_factor errors", {
  x <- c(1, 2, 3)
  expect_error(estimate_noise(x, sd_factor = -1), "non-negative")
  expect_error(estimate_noise(x, sd_factor = c(1, 2)), "single")
  expect_error(estimate_noise(x, sd_factor = NA_real_))
  expect_error(estimate_noise(x, sd_factor = "3"))
  expect_error(estimate_noise(x, sd_factor = Inf))
})
