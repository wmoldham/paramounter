# test-to_mzmine.R

full_dists <- function(...) utils::modifyList(
  stats::setNames(rep(list(numeric(0)), 9),
                  c("ppm", "mz_diff", "noise", "width_seconds", "width_scans", "sn", "height", "mass_shift", "rt_shift")),
  list(...))

test_that("to_mzmine reproduces the original MZmine values", {
  orig <- function(d, multi) {
    mn <- floor(min(d$noise)); ps <- floor(min(d$width_scans)); ph <- floor(min(d$height)); mp <- ceiling(max(d$ppm))
    W <- mean(d$width_seconds, trim = 0.05); H <- mean(d$height, trim = 0.05); ratio <- H / W
    lo <- min(d$width_seconds); hi <- max(d$width_seconds)
    # these tests build universal_parameters with the default config, so legacy is
    # FALSE and the lower bound stays data-driven in both branches
    wide <- hi > 35 & ratio > 515
    lo <- ceiling(lo) + 4
    hi <- if (wide) (ceiling(hi) + 7) / 2 else ceiling(hi) + 5
    if (multi) round(c(mn, ps, mn, ph, mp, lo / 60, hi / 60, max(d$mass_shift), max(d$rt_shift) / 60), 3)
    else round(c(mn, ps, mn, ph, mp, lo / 60, hi / 60), 3)
  }
  gen <- function(branch, multi) {
    n <- sample(40:120, 1)
    ws <- if (branch == "wide") runif(n, 5, 60) else runif(n, 3, 28)
    ht <- if (branch == "lowratio") runif(n, 5, 40) else runif(n, 1e4, 1e6)
    full_dists(ppm = runif(n, 1, 8), noise = runif(n, 100, 5000), width_seconds = ws, width_scans = sample(3:40, n, TRUE),
               height = ht, mass_shift = if (multi) runif(sample(5:30, 1), 0.001, 0.01) else numeric(0),
               rt_shift = if (multi) runif(sample(5:30, 1), 1, 30) else numeric(0))
  }
  set.seed(12)
  for (rep in 1:200) {
    multi <- sample(c(TRUE, FALSE), 1); d <- gen(sample(c("wide", "narrow", "lowratio"), 1), multi)
    up <- universal_parameters(distributions = d, files = if (multi) c("a", "b") else "a")
    tbl <- to_mzmine(up)
    expect_equal(tbl$value, orig(d, multi))
    expect_equal(nrow(tbl), if (multi) 9L else 7L)
  }
})

test_that("to_mzmine validates and can write a file", {
  expect_error(to_mzmine("x"), "universal_parameters object")
  d <- full_dists(ppm = runif(20, 1, 8), noise = runif(20, 100, 5000), width_seconds = runif(20, 5, 30),
                  width_scans = sample(3:20, 20, TRUE), height = runif(20, 1e4, 1e6))
  up <- universal_parameters(distributions = d, files = "a")
  f <- tempfile(fileext = ".csv")
  to_mzmine(up, file = f)
  expect_true(file.exists(f))
  expect_equal(nrow(utils::read.csv(f)), 7L)
})
