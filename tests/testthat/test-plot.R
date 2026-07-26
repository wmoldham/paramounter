# test-plot.R

full_dists <- function(...) utils::modifyList(
  stats::setNames(rep(list(numeric(0)), 9),
    c("ppm", "mz_diff", "noise", "width_seconds", "width_scans", "sn", "height", "mass_shift", "rt_shift")),
  list(...))

test_that("plot runs on measured distributions", {
  set.seed(9)
  d <- full_dists(
    ppm = runif(100, 1, 8), mz_diff = runif(100, 0.001, 0.01), noise = runif(100, 100, 5000),
    width_seconds = runif(100, 4, 30), width_scans = sample(3:20, 100, TRUE), sn = runif(100, 3, 40),
    height = 10^runif(100, 3, 6), mass_shift = runif(10, 0.001, 0.01), rt_shift = runif(10, 1, 30)
  )
  up <- universal_parameters(distributions = d, files = c("a", "b"))
  grDevices::pdf(tempfile(fileext = ".pdf"))
  on.exit(grDevices::dev.off())
  expect_no_error(plot(up))
})

test_that("plot omits empty quantities and handles constant values", {
  set.seed(1)
  up_one <- universal_parameters(distributions = full_dists(ppm = runif(50, 1, 8)), files = "a")
  up_const <- universal_parameters(
    distributions = full_dists(ppm = rep(4, 30), noise = rep(100, 30), width_seconds = rep(10, 30),
                               width_scans = rep(5L, 30), sn = rep(6, 30), height = rep(1e4, 30)),
    files = "a")
  grDevices::pdf(tempfile(fileext = ".pdf"))
  on.exit(grDevices::dev.off())
  expect_no_error(plot(up_one))
  expect_no_error(plot(up_const))
})

test_that("plot errors when there is nothing to show", {
  grDevices::pdf(tempfile(fileext = ".pdf"))
  on.exit(grDevices::dev.off())
  expect_error(plot(universal_parameters()), "No non-empty")
})
