# test-pm_config.R

test_that("defaults match the published method", {
  d <- pm_config()
  expect_equal(d@bin_width, 0.05)
  expect_equal(d@smooth_half_window, 0L)
  expect_equal(d@noise_block_size, 10L)
  expect_equal(d@noise_sd_factor, 3)
  expect_equal(d@mass_sd_range, 2)
  expect_null(d@ppm_cutoff)
  expect_equal(d@ppm_quantile, 0.95)
  expect_equal(d@min_masses, 2L)
  expect_equal(d@max_masses, 199L)
  expect_equal(d@isolation_gap, 300)
  expect_equal(d@match_mz_tol, 0.015)
  expect_equal(d@match_rt_tol, 30)
  expect_equal(d@trim, 0.97)
  expect_false(d@legacy)
})

test_that("custom values are stored", {
  c1 <- pm_config(bin_width = 0.02, legacy = TRUE, ppm_cutoff = 20, smooth_half_window = 4L)
  expect_equal(c1@bin_width, 0.02)
  expect_true(c1@legacy)
  expect_equal(c1@ppm_cutoff, 20)
  expect_equal(c1@smooth_half_window, 4L)
})

test_that("ppm_cutoff accepts NULL (auto) or a positive number", {
  expect_null(pm_config(ppm_cutoff = NULL)@ppm_cutoff)
  expect_equal(pm_config(ppm_cutoff = 15)@ppm_cutoff, 15)
  expect_error(pm_config(ppm_cutoff = -1), "ppm_cutoff.*positive")
  expect_error(pm_config(ppm_cutoff = 0), "positive")
})

test_that("numeric settings are validated", {
  expect_error(pm_config(bin_width = -1), "bin_width.*positive")
  expect_error(pm_config(bin_width = 0), "positive")
  expect_error(pm_config(bin_width = c(1, 2)), "single")
  expect_error(pm_config(noise_sd_factor = -1), "non-negative")
  expect_error(pm_config(mass_sd_range = 0), "positive")
})

test_that("integer settings are validated", {
  expect_error(pm_config(smooth_half_window = 2.5), "integer")
  expect_error(pm_config(smooth_half_window = -1), "non-negative")
  expect_error(pm_config(noise_block_size = 0), "positive integer")
})

test_that("proportion settings must lie in (0, 1]", {
  expect_error(pm_config(ppm_quantile = 0), "0, 1")
  expect_error(pm_config(ppm_quantile = 1.5), "0, 1")
  expect_error(pm_config(trim = 0), "0, 1")
  expect_error(pm_config(trim = 1.5), "0, 1")
  expect_silent(pm_config(trim = 1))
  expect_silent(pm_config(ppm_quantile = 1))
})

test_that("legacy must be a flag", {
  expect_error(pm_config(legacy = NA), "TRUE or FALSE")
  expect_error(pm_config(legacy = 1), "TRUE or FALSE")
})

test_that("min_masses must not exceed max_masses", {
  expect_error(pm_config(min_masses = 10, max_masses = 5), "exceed")
})

test_that("settings are re-validated on modification", {
  cc <- pm_config()
  expect_error({cc@bin_width <- -5}, "positive")
  cc@bin_width <- 0.1
  expect_equal(cc@bin_width, 0.1)
})

test_that("printing shows the settings and the auto ppm cutoff", {
  expect_output(print(pm_config()), "pm_config")
  expect_output(print(pm_config()), "auto")
  expect_output(print(pm_config(ppm_cutoff = 20)), "20")
})

test_that("tolerance zone settings default to every zone and are validated", {
  d <- pm_config()
  expect_false(d@tolerance_isolated)
  expect_equal(d@tolerance_min_scans, 1L)
  expect_true(pm_config(tolerance_isolated = TRUE)@tolerance_isolated)
  expect_equal(pm_config(tolerance_min_scans = 5L)@tolerance_min_scans, 5L)
  expect_error(pm_config(tolerance_isolated = NA), "TRUE or FALSE")
  expect_error(pm_config(tolerance_min_scans = 0), "positive integer")
  expect_error(pm_config(tolerance_min_scans = 2.5), "integer")
})

test_that("printing shows the tolerance settings", {
  out <- capture.output(print(pm_config(tolerance_isolated = TRUE, tolerance_min_scans = 5L)))
  expect_true(any(grepl("tolerance_quantile +0.95", out)))
  expect_true(any(grepl("tolerance_isolated +TRUE", out)))
  expect_true(any(grepl("tolerance_min_scans +5", out)))
})
