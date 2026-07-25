# test-universal_parameters.R

full_dists <- function(...) {
  base <- stats::setNames(
    rep(list(numeric(0)), 9),
    c("ppm", "mz_diff", "noise", "width_seconds", "width_scans",
      "sn", "height", "mass_shift", "rt_shift")
  )
  utils::modifyList(base, list(...))
}

test_that("the summary is computed from the distributions", {
  up <- universal_parameters(
    distributions = full_dists(ppm = c(3, 5, 4, 6, 5), noise = c(100, 120, 110)),
    files = c("a.mzML", "b.mzML")
  )
  s <- up@summary
  expect_s3_class(s, "data.frame")
  expect_equal(nrow(s), 9L)
  ppm <- s[s$quantity == "ppm", ]
  expect_equal(ppm$n, 5L)
  expect_equal(ppm$min, 3)
  expect_equal(ppm$max, 6)
  expect_equal(ppm$mean, 4.6)
  expect_equal(ppm$median, 5)
})

test_that("empty quantities give n zero and NA statistics", {
  up <- universal_parameters(distributions = full_dists(ppm = c(1, 2)))
  rt <- up@summary[up@summary$quantity == "rt_shift", ]
  expect_equal(rt$n, 0L)
  expect_true(is.na(rt$min))
  expect_true(is.na(rt$mean))
})

test_that("the summary tracks changes to the distributions", {
  up <- universal_parameters(distributions = full_dists(ppm = c(3, 5, 4)))
  expect_equal(up@summary[up@summary$quantity == "ppm", ]$max, 5)
  up@distributions$ppm <- c(1, 2, 3)
  expect_equal(up@summary[up@summary$quantity == "ppm", ]$max, 3)
})

test_that("provenance is stored and the config type is enforced", {
  up <- universal_parameters(
    distributions = full_dists(),
    files = c("a.mzML", "b.mzML"),
    config = pm_config(legacy = FALSE)
  )
  expect_identical(up@files, c("a.mzML", "b.mzML"))
  expect_false(up@config@legacy)
  expect_true(S7::S7_inherits(up@config, pm_config))
})

test_that("an empty object is valid", {
  e <- universal_parameters()
  expect_equal(nrow(e@summary), 9L)
  expect_true(all(e@summary$n == 0L))
  expect_length(e@files, 0L)
})

test_that("distributions are validated", {
  expect_error(universal_parameters(distributions = list(1, 2)), "named list")
  expect_error(universal_parameters(distributions = list(ppm = 1)), "exactly these names")
  expect_error(
    universal_parameters(distributions = full_dists(ppm = c(1, NA))),
    "distributions\\$ppm.*finite"
  )
  expect_error(
    universal_parameters(distributions = full_dists(noise = "a")),
    "distributions\\$noise.*finite"
  )
})

test_that("files and config are validated", {
  expect_error(universal_parameters(files = 42), "character vector")
  expect_error(universal_parameters(config = "nope"))
})

test_that("summary is read-only", {
  up <- universal_parameters()
  expect_error({up@summary <- 1})
})

test_that("printing shows the file count and the parameter table", {
  up <- universal_parameters(distributions = full_dists(ppm = c(3, 5)), files = "a.mzML")
  expect_output(print(up), "universal_parameters")
  expect_output(print(up), "ppm")
})
