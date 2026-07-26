# test-to_msdial.R

full_dists <- function(...) utils::modifyList(
  stats::setNames(rep(list(numeric(0)), 9),
                  c("ppm", "mz_diff", "noise", "width_seconds", "width_scans", "sn", "height", "mass_shift", "rt_shift")),
  list(...))

test_that("to_msdial reproduces the original MS-DIAL values", {
  orig <- function(d, multi) {
    mzd <- ceiling(max(d$mz_diff) * 100) / 100; ph <- floor(min(d$height)); ps <- floor(min(d$width_scans))
    if (multi) round(c(mzd, ph, mzd, ps, max(d$mass_shift), max(d$rt_shift) / 60), 3) else round(c(mzd, ph, mzd, ps), 3)
  }
  set.seed(11)
  for (rep in 1:200) {
    multi <- sample(c(TRUE, FALSE), 1); n <- sample(40:120, 1)
    d <- full_dists(mz_diff = runif(n, 0.001, 0.01), height = runif(n, 1e4, 1e6), width_scans = sample(3:40, n, TRUE),
                    mass_shift = if (multi) runif(sample(5:30, 1), 0.001, 0.01) else numeric(0),
                    rt_shift = if (multi) runif(sample(5:30, 1), 1, 30) else numeric(0))
    up <- universal_parameters(distributions = d, files = if (multi) c("a", "b") else "a")
    tbl <- to_msdial(up)
    expect_equal(tbl$value, orig(d, multi))
    expect_equal(nrow(tbl), if (multi) 6L else 4L)
  }
})

test_that("to_msdial validates and can write a file", {
  expect_error(to_msdial("x"), "universal_parameters object")
  d <- full_dists(height = 1:5 * 1000, width_scans = sample(3:9, 5, TRUE))
  expect_error(to_msdial(universal_parameters(distributions = d, files = "a")), "no measurements for mz_diff")
  d2 <- full_dists(mz_diff = runif(20, 0.001, 0.01), height = runif(20, 1e4, 1e6), width_scans = sample(3:20, 20, TRUE))
  up <- universal_parameters(distributions = d2, files = "a")
  f <- tempfile(fileext = ".csv")
  to_msdial(up, file = f)
  expect_true(file.exists(f))
  expect_equal(nrow(utils::read.csv(f)), 4L)
})
