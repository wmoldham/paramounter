# test-zoi_features.R

test_that("a valid zone of interest produces the expected measurements", {
  rt <- seq(10, 28, by = 2)
  masses <- c(220.0098, 220.0100, 220.0102, 220.0101, 220.0099)
  r <- zoi_features(masses, 220.0100, 2, 6, 4, rt, 5000, 100, 30, 120)
  mass_sd <- sd(masses)
  expect_equal(r$reference_mz, 220.0100)
  expect_equal(r$apex_rt, 16)
  expect_equal(r$ppm, 2 * mass_sd / 220.0100 * 1e6)
  expect_equal(r$mz_diff, 2 * mass_sd)
  expect_equal(r$width_seconds, 8)
  expect_identical(r$width_scans, 6L)
  expect_equal(r$sn, (5000 - 100) / 30)
  expect_equal(r$height, 5000)
})

test_that("too few or too many masses are filtered to NULL", {
  rt <- seq_len(400)
  expect_null(zoi_features(220.01, 220.01, 2, 6, 4, seq(10, 28, by = 2), 5000, 100, 30, 120))
  expect_false(is.null(
    zoi_features(c(220.010, 220.011), 220.01, 2, 6, 4, seq(10, 28, by = 2), 5000, 100, 30, 120)
  ))
  set.seed(1)
  expect_false(is.null(
    zoi_features(220 + rnorm(199, 0, 1e-4), 220, 2, 250, 100, rt, 5000, 100, 30, 120)
  ))
  expect_null(
    zoi_features(220 + rnorm(200, 0, 1e-4), 220, 2, 250, 100, rt, 5000, 100, 30, 120)
  )
})

test_that("a relative tolerance at or above the cutoff is filtered to NULL", {
  rt <- seq(10, 28, by = 2)
  wide <- c(220.00, 220.02, 220.04, 220.01, 220.03)
  expect_null(zoi_features(wide, 220.00, 2, 6, 4, rt, 5000, 100, 30, 120, ppm_cutoff = 5))
  expect_false(is.null(
    zoi_features(wide, 220.00, 2, 6, 4, rt, 5000, 100, 30, 120, ppm_cutoff = Inf)
  ))
})

test_that("signal-to-noise falls back to the apex intensity when cutoff is zero", {
  rt <- seq(10, 28, by = 2)
  masses <- c(220.0098, 220.0100, 220.0102, 220.0101, 220.0099)
  expect_equal(zoi_features(masses, 220.0100, 2, 6, 4, rt, 5000, 100, 30, 0)$sn, 5000)
})

test_that("signal-to-noise is infinite when the noise standard deviation is zero", {
  rt <- seq(10, 28, by = 2)
  masses <- c(220.0098, 220.0100, 220.0102, 220.0101, 220.0099)
  expect_true(is.infinite(zoi_features(masses, 220.0100, 2, 6, 4, rt, 5000, 100, 0, 120)$sn))
})

test_that("matches the original feature arithmetic across random zones", {
  orig_features <- function(masses, refMz, leftInd, rightInd, apex, eicApex, rtime,
                            blk, sdv, cutOFF, massSDrange, ppmCut) {
    if (!(length(masses) > 1 && length(masses) < 200)) return(NULL)
    ppmCheck <- (massSDrange * sd(masses)) / refMz * 1e6
    if (!(ppmCheck < ppmCut)) return(NULL)
    list(
      ppm = ppmCheck,
      da = massSDrange * sd(masses),
      width = rtime[rightInd - 1] - rtime[leftInd + 1],
      scans = rightInd - leftInd,
      sn = if (cutOFF > 0) (eicApex - blk) / sdv else eicApex,
      height = eicApex,
      apex_rt = rtime[apex]
    )
  }
  set.seed(31)
  for (rep in 1:800) {
    n_rt <- sample(30:200, 1)
    rtime <- cumsum(c(runif(1, 0, 5), runif(n_rt - 1, 0.2, 3)))
    center <- sample(c(120, 220, 450, 800), 1)
    nm <- sample(1:60, 1)
    masses <- center + rnorm(nm, 0, center * runif(1, 1e-6, 3e-5))
    left_edge <- sample(1:(n_rt - 2), 1)
    right_edge <- min(left_edge + sample(1:20, 1), n_rt)
    apex <- sample(left_edge:right_edge, 1)
    apex_int <- runif(1, 500, 1e5)
    nmean <- runif(1, 0, 500)
    nsd <- if (runif(1) < 0.1) 0 else runif(1, 1, 200)
    cutoff <- runif(1, 1, 400)
    ppm_cut <- sample(c(Inf, runif(1, 1, 30)), 1)
    z <- zoi_features(
      masses, center, left_edge, right_edge, apex, rtime,
      apex_int, nmean, nsd, cutoff, sd_range = 2, ppm_cutoff = ppm_cut
    )
    o <- orig_features(
      masses, center, left_edge - 1, right_edge + 1, apex,
      apex_int, rtime, nmean, nsd, cutoff, 2, ppm_cut
    )
    if (is.null(o)) {
      expect_null(z)
    } else {
      expect_equal(z$ppm, o$ppm)
      expect_equal(z$mz_diff, o$da)
      expect_equal(z$width_seconds, o$width)
      expect_equal(z$width_scans, o$scans)
      expect_equal(z$sn, o$sn)
      expect_equal(z$height, o$height)
      expect_equal(z$apex_rt, o$apex_rt)
    }
  }
})

test_that("vector and index inputs are validated", {
  rt <- seq(10, 28, by = 2)
  m <- c(220.0098, 220.0100, 220.0102)
  expect_error(zoi_features("x", 220, 2, 6, 4, rt, 5000, 100, 30, 120), "masses")
  expect_error(zoi_features(c(220, NA), 220, 2, 6, 4, rt, 5000, 100, 30, 120), "NA")
  expect_error(zoi_features(m, 220, 2, 6, 4, "x", 5000, 100, 30, 120), "rtime")
  expect_error(zoi_features(m, 220, 2, 6, 4, c(1, NA, 3), 5000, 100, 30, 120), "NA")
  expect_error(zoi_features(m, 220, 0, 6, 4, rt, 5000, 100, 30, 120), "within")
  expect_error(zoi_features(m, 220, 2, 99, 4, rt, 5000, 100, 30, 120), "within")
  expect_error(zoi_features(m, 220, 6, 2, 4, rt, 5000, 100, 30, 120), "exceed")
  expect_error(zoi_features(m, 220, 2, 6, 8, rt, 5000, 100, 30, 120), "between")
  expect_error(zoi_features(m, 220, 2.5, 6, 4, rt, 5000, 100, 30, 120), "integer")
})

test_that("scalar inputs are validated", {
  rt <- seq(10, 28, by = 2)
  m <- c(220.0098, 220.0100, 220.0102)
  expect_error(zoi_features(m, 0, 2, 6, 4, rt, 5000, 100, 30, 120), "positive")
  expect_error(zoi_features(m, c(220, 221), 2, 6, 4, rt, 5000, 100, 30, 120), "single")
  expect_error(zoi_features(m, 220, 2, 6, 4, rt, -5, 100, 30, 120), "non-negative")
  expect_error(zoi_features(m, 220, 2, 6, 4, rt, 5000, -1, 30, 120), "non-negative")
  expect_error(zoi_features(m, 220, 2, 6, 4, rt, 5000, 100, -1, 120), "non-negative")
  expect_error(zoi_features(m, 220, 2, 6, 4, rt, 5000, 100, 30, -1), "non-negative")
  expect_error(zoi_features(m, 220, 2, 6, 4, rt, 5000, 100, 30, 120, sd_range = 0), "positive")
  expect_error(zoi_features(m, 220, 2, 6, 4, rt, 5000, 100, 30, 120, ppm_cutoff = 0), "positive")
  expect_error(zoi_features(m, 220, 2, 6, 4, rt, 5000, 100, 30, 120, ppm_cutoff = NA_real_), "positive")
  expect_error(zoi_features(m, 220, 2, 6, 4, rt, 5000, 100, 30, 120, min_masses = 0), "positive")
  expect_error(zoi_features(m, 220, 2, 6, 4, rt, 5000, 100, 30, 120, min_masses = 10, max_masses = 5), "exceed")
})

test_that("the scan-count toggle switches between legacy and corrected spans", {
  rt <- seq(10, 28, by = 2)
  masses <- c(220.0098, 220.0100, 220.0102, 220.0101, 220.0099)
  legacy <- zoi_features(masses, 220.0100, 2, 6, 4, rt, 5000, 100, 30, 120)
  corrected <- zoi_features(
    masses, 220.0100, 2, 6, 4, rt, 5000, 100, 30, 120,
    legacy_scan_count = FALSE
  )
  expect_identical(legacy$width_scans, 6L) # right - left + 2, as published
  expect_identical(corrected$width_scans, 5L) # inclusive span
  expect_identical(
    legacy[names(legacy) != "width_scans"],
    corrected[names(corrected) != "width_scans"]
  )
})

test_that("legacy_scan_count must be a flag", {
  rt <- seq(10, 28, by = 2)
  m <- c(220.0098, 220.0100, 220.0102)
  expect_error(
    zoi_features(m, 220, 2, 6, 4, rt, 5000, 100, 30, 120, legacy_scan_count = NA),
    "TRUE or FALSE"
  )
  expect_error(
    zoi_features(m, 220, 2, 6, 4, rt, 5000, 100, 30, 120, legacy_scan_count = 1),
    "TRUE or FALSE"
  )
})
