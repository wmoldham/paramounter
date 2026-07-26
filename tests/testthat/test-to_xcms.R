# test-to_xcms.R

full_dists <- function(...) utils::modifyList(
  stats::setNames(rep(list(numeric(0)), 9),
                  c("ppm", "mz_diff", "noise", "width_seconds", "width_scans", "sn", "height", "mass_shift", "rt_shift")),
  list(...)
)

test_that("xcms_values reproduces the original point estimates", {
  orig_xcms_values <- function(d, legacy = TRUE) {
    maxppm <- ceiling(max(d$ppm)); minnoise <- floor(min(d$noise))
    W <- mean(d$width_seconds, trim = 0.05, na.rm = TRUE); H <- mean(d$height, trim = 0.05, na.rm = TRUE); ratio <- H / W
    minpeakwidth <- min(d$width_seconds); maxpeakwidth <- max(d$width_seconds)
    if (maxpeakwidth > 35 & ratio > 515) { minpeakwidth <- 0; maxpeakwidth <- (ceiling(maxpeakwidth) + 7) / 2 }
    else { minpeakwidth <- ceiling(minpeakwidth) + 4; maxpeakwidth <- ceiling(maxpeakwidth) + 5 }
    minpeakscan <- floor(min(d$width_scans)); minSN <- min(d$sn); if (minSN < 3) minSN <- 3
    bw <- 5; if (!legacy && length(d$rt_shift) > 0) bw <- max(d$rt_shift)
    list(ppm = maxppm, min_pw = minpeakwidth, max_pw = maxpeakwidth, snthresh = minSN,
         pf_k = minpeakscan, pf_I = minnoise, noise = minnoise, bw = bw, binSize = max(d$mass_shift))
  }
  gen <- function(branch) {
    n <- sample(50:150, 1)
    ws <- if (branch == "wide") runif(n, 5, 60) else runif(n, 3, 28)
    ht <- if (branch == "lowratio") runif(n, 5, 40) else runif(n, 1e4, 1e6)
    full_dists(ppm = runif(n, 1, 8), noise = runif(n, 100, 5000), width_seconds = ws,
               width_scans = sample(3:40, n, TRUE), sn = runif(n, 3, 50), height = ht,
               mass_shift = runif(sample(5:30, 1), 0.001, 0.01), rt_shift = runif(sample(5:30, 1), 1, 30))
  }
  set.seed(505)
  for (rep in 1:300) {
    d <- gen(sample(c("wide", "narrow", "lowratio"), 1)); legacy <- sample(c(TRUE, FALSE), 1)
    up <- universal_parameters(distributions = d, files = c("a", "b"), config = pm_config(legacy = legacy))
    v <- xcms_values(up); a <- orig_xcms_values(d, legacy)
    expect_equal(v$ppm, a$ppm)
    expect_equal(v$peakwidth, c(a$min_pw, a$max_pw))
    expect_equal(v$snthresh, a$snthresh)
    expect_equal(v$prefilter, c(a$pf_k, a$pf_I))
    expect_equal(v$noise, a$noise)
    expect_equal(v$bw, a$bw)
    expect_equal(v$binSize, a$binSize)
  }
})

test_that("xcms_values validates, floors bw under legacy, and honours sample_groups", {
  d <- full_dists(ppm = runif(50, 1, 8), noise = runif(50, 100, 5000), width_seconds = runif(50, 5, 45),
                  width_scans = sample(3:40, 50, TRUE), sn = runif(50, 3, 50), height = runif(50, 1e4, 1e6),
                  mass_shift = runif(10, 0.001, 0.01), rt_shift = runif(10, 1, 30))
  expect_equal(xcms_values(universal_parameters(distributions = d, files = c("a", "b"), config = pm_config(legacy = TRUE)))$bw, 5)
  expect_gt(xcms_values(universal_parameters(distributions = d, files = c("a", "b"), config = pm_config(legacy = FALSE)))$bw, 0)
  up <- universal_parameters(distributions = d, files = c("a", "b"))
  expect_identical(xcms_values(up, sample_groups = c("A", "B"))$sample_groups, c("A", "B"))

  d_empty <- full_dists(noise = runif(5, 1, 10), width_seconds = runif(5, 1, 10),
                        width_scans = sample(3:9, 5, TRUE), sn = runif(5, 3, 10), height = runif(5, 1e3, 1e4))
  up_empty <- universal_parameters(distributions = d_empty, files = "a")
  expect_error(xcms_values(up_empty), "no measurements for ppm")
})

test_that("to_xcms validates its inputs", {
  d <- full_dists(ppm = runif(20, 1, 8), noise = runif(20, 100, 5000), width_seconds = runif(20, 5, 30),
                  width_scans = sample(3:20, 20, TRUE), sn = runif(20, 3, 50), height = runif(20, 1e4, 1e6),
                  mass_shift = runif(5, 0.001, 0.01), rt_shift = runif(5, 1, 30))
  up <- universal_parameters(distributions = d, files = c("a", "b"))
  expect_error(to_xcms("x"), "universal_parameters object")
  expect_error(to_xcms(up, sample_groups = "A"), "one entry per file")
})

test_that("to_xcms builds the xcms parameter objects", {
  skip_if_not_installed("xcms")
  d <- full_dists(ppm = runif(40, 1, 8), noise = runif(40, 100, 5000), width_seconds = runif(40, 5, 30),
                  width_scans = sample(3:20, 40, TRUE), sn = runif(40, 3, 50), height = runif(40, 1e4, 1e6),
                  mass_shift = runif(10, 0.001, 0.01), rt_shift = runif(10, 1, 30))
  up <- universal_parameters(distributions = d, files = c("a", "b"), config = pm_config(legacy = FALSE))
  xp <- to_xcms(up)
  expect_named(xp, c("chrom_peaks", "group", "retention"))
  expect_s4_class(xp$chrom_peaks, "CentWaveParam")
  expect_s4_class(xp$group, "PeakDensityParam")
  expect_s4_class(xp$retention, "ObiwarpParam")
  expect_equal(xp$chrom_peaks@ppm, ceiling(max(d$ppm)))
})

test_that("to_xcms warns and falls back when shifts are unavailable", {
  skip_if_not_installed("xcms")
  d <- full_dists(ppm = runif(20, 1, 8), noise = runif(20, 100, 5000), width_seconds = runif(20, 5, 30),
                  width_scans = sample(3:20, 20, TRUE), sn = runif(20, 3, 50), height = runif(20, 1e4, 1e6))
  up <- universal_parameters(distributions = d, files = "a")
  expect_warning(to_xcms(up), "shifts unavailable")
})
