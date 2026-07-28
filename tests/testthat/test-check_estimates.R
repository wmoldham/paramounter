# test-check_estimates.R

full_dists <- function(...) utils::modifyList(
  stats::setNames(rep(list(numeric(0)), 9),
                  c("ppm", "mz_diff", "noise", "width_seconds", "width_scans", "sn", "height", "mass_shift", "rt_shift")),
  list(...))

make_params <- function(d, ...) {
  universal_parameters(
    distributions = d, files = c("a", "b"), config = pm_config(...)
  )
}

test_that("the reported estimate is the one the translators will use", {
  set.seed(1)
  d <- full_dists(
    ppm = runif(200, 1, 8), mz_diff = runif(200, 0.001, 0.01),
    noise = runif(200, 100, 5000), width_seconds = runif(200, 3, 25),
    width_scans = sample(3:40, 200, TRUE), sn = runif(200, 3, 50),
    height = runif(200, 1e4, 1e6), mass_shift = runif(30, 0.001, 0.01),
    rt_shift = runif(30, 1, 30)
  )
  p <- make_params(d)
  chk <- check_estimates(p)
  v <- xcms_values(p)

  expect_s3_class(chk, "data.frame")
  expect_named(chk, c("quantity", "kind", "estimate", "median", "ratio"))

  pick <- function(q) chk$estimate[chk$quantity == q]
  expect_equal(pick("ppm"), v$ppm)
  expect_equal(pick("noise"), v$noise)
  expect_equal(pick("sn"), v$snthresh)
  expect_equal(pick("width_scans"), v$prefilter[1])
  # width_seconds reports the upper peak-width bound, not the lower
  expect_equal(pick("width_seconds"), v$peakwidth[2])
})

test_that("ratio measures distance from the bulk in the right direction", {
  set.seed(2)
  d <- full_dists(
    ppm = runif(200, 1, 8), noise = runif(200, 100, 5000),
    width_seconds = runif(200, 3, 25), width_scans = sample(3:40, 200, TRUE),
    sn = runif(200, 3, 50), height = runif(200, 1e4, 1e6)
  )
  chk <- check_estimates(make_params(d))
  # upper-end quantities: estimate above the median, so ratio > 1
  expect_gt(chk$ratio[chk$quantity == "ppm"], 1)
  # lower-end quantities: estimate below the median, and ratio is inverted so it
  # still reads as "how many times away from the bulk"
  expect_gt(chk$ratio[chk$quantity == "noise"], 1)
  expect_lt(chk$estimate[chk$quantity == "noise"], chk$median[chk$quantity == "noise"])

  for (q in c("ppm", "noise", "sn", "height", "width_scans", "width_seconds")) {
    expect_equal(
      chk$median[chk$quantity == q], stats::median(d[[q]]),
      info = q
    )
  }
})

test_that("a heavy tail shows up as a large ratio, and legacy makes it worse", {
  set.seed(3)
  rest <- list(
    noise = runif(500, 100, 5000), width_seconds = runif(500, 3, 25),
    width_scans = sample(3:40, 500, TRUE), sn = runif(500, 3, 50),
    height = runif(500, 1e4, 1e6)
  )
  tight <- do.call(full_dists, c(list(ppm = runif(500, 1.8, 2.2)), rest))
  # the Orbitrap shape: a tight bulk around 2 ppm with a long upper tail
  heavy <- do.call(full_dists, c(list(ppm = stats::rlnorm(500, log(2), 1.2)), rest))

  r <- function(d, legacy) {
    chk <- check_estimates(make_params(d, legacy = legacy))
    chk$ratio[chk$quantity == "ppm"]
  }

  # a tight distribution puts the estimate right next to the bulk
  expect_lt(r(tight, FALSE), 2)

  # a heavy tail pushes it away, which is the signal the user needs to see
  expect_gt(r(heavy, FALSE), 3)

  # and the legacy estimator, which takes the maximum, is pushed much further --
  # this is the whole reason tolerance_quantile exists
  expect_gt(r(heavy, TRUE), r(heavy, FALSE))
})

test_that("tolerances and thresholds are labelled", {
  set.seed(4)
  d <- full_dists(
    ppm = runif(50, 1, 8), mz_diff = runif(50, 0.001, 0.01),
    noise = runif(50, 100, 5000), width_seconds = runif(50, 3, 25),
    width_scans = sample(3:40, 50, TRUE), sn = runif(50, 3, 50),
    height = runif(50, 1e4, 1e6)
  )
  chk <- check_estimates(make_params(d))
  kind <- function(q) chk$kind[chk$quantity == q]
  expect_identical(kind("ppm"), "tolerance")
  expect_identical(kind("mz_diff"), "tolerance")
  expect_identical(kind("width_seconds"), "tolerance")
  expect_identical(kind("noise"), "threshold")
  expect_identical(kind("sn"), "threshold")
  expect_identical(kind("height"), "threshold")
  expect_identical(kind("width_scans"), "threshold")
})

test_that("unmeasured quantities are omitted and inputs are validated", {
  set.seed(5)
  # no mass_shift or rt_shift, as with a single file
  d <- full_dists(
    ppm = runif(50, 1, 8), noise = runif(50, 100, 5000),
    width_seconds = runif(50, 3, 25), width_scans = sample(3:40, 50, TRUE),
    sn = runif(50, 3, 50), height = runif(50, 1e4, 1e6)
  )
  chk <- check_estimates(make_params(d))
  expect_false(any(c("mass_shift", "rt_shift", "mz_diff") %in% chk$quantity))
  expect_equal(nrow(chk), 6L)

  expect_equal(nrow(check_estimates(make_params(full_dists()))), 0L)
  expect_error(check_estimates("x"), "universal_parameters object")
})
