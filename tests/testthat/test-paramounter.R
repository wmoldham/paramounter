# test-paramounter.R

make_multi_traces <- function(n_files, n_shared, n_scans, mz_lo = 100, mz_hi = 110) {
  smz <- sort(runif(n_shared, mz_lo + 1, mz_hi - 1)); srt <- runif(n_shared, 10, n_scans - 10); sw <- runif(n_shared, 2.5, 6); shh <- runif(n_shared, 1e4, 1e6)
  lapply(seq_len(n_files), function(fi) {
    dmz <- rnorm(1, 0, 0.001); drt <- rnorm(1, 0, 0.5); mzD <- vector("list", n_scans); intD <- vector("list", n_scans)
    for (s in seq_len(n_scans)) { pm <- numeric(0); pv <- numeric(0)
    for (k in seq_len(n_shared)) { g <- shh[k] * exp(-0.5 * ((s - (srt[k] + drt)) / sw[k])^2)
    if (g > shh[k] * 0.02) { pm <- c(pm, smz[k] + dmz + rnorm(1, 0, smz[k] * 2e-6)); pv <- c(pv, g) } }
    nn <- rpois(1, 2); pm <- c(pm, runif(nn, mz_lo, mz_hi)); pv <- c(pv, runif(nn, 1e2, 1.5e3)); mzD[[s]] <- pm; intD[[s]] <- pv }
    mzD[[1]] <- c(mzD[[1]], mz_lo + 0.007); intD[[1]] <- c(intD[[1]], 300); mzD[[2]] <- c(mzD[[2]], mz_hi - 0.007); intD[[2]] <- c(intD[[2]], 300)
    for (s in seq_len(n_scans)) { o <- order(mzD[[s]]); mzD[[s]] <- mzD[[s]][o]; intD[[s]] <- intD[[s]][o] }
    list(mz = mzD, intensity = intD, rtime = seq(0, by = 6, length.out = n_scans)) })
}

# temp files whose contents are ignored by the fake reader (they only need to exist)
make_temp_files <- function(n) vapply(seq_len(n), function(i) { f <- tempfile(fileext = ".mzML"); file.create(f); f }, character(1))

test_that("paramounter matches a manual read + measure + aggregate", {
  set.seed(404)
  traces <- make_multi_traces(3, 16, 80)
  files <- make_temp_files(3)
  by_file <- stats::setNames(traces, files)
  reader <- function(f) by_file[[f]]
  config <- pm_config()

  result <- paramounter(files, config, reader = reader)
  per_file <- lapply(traces, function(t) measure_file(t$mz, t$intensity, t$rtime, config))
  manual <- aggregate_files(per_file, files, config)

  expect_true(S7::S7_inherits(result, universal_parameters))
  expect_identical(result@files, files)
  for (q in names(result@distributions)) {
    expect_equal(result@distributions[[q]], manual@distributions[[q]])
  }
})

test_that("the reader is called once per file", {
  set.seed(1)
  traces <- make_multi_traces(3, 12, 70)
  files <- make_temp_files(3)
  by_file <- stats::setNames(traces, files)
  calls <- 0L
  reader <- function(f) { calls <<- calls + 1L; by_file[[f]] }
  paramounter(files, reader = reader)
  expect_equal(calls, 3L)
})

test_that("a single file yields empty shift distributions", {
  set.seed(2)
  traces <- make_multi_traces(1, 12, 70)
  files <- make_temp_files(1)
  by_file <- stats::setNames(traces, files)
  result <- paramounter(files, reader = function(f) by_file[[f]])
  expect_length(result@distributions$mass_shift, 0L)
  expect_length(result@distributions$rt_shift, 0L)
  expect_gt(length(result@distributions$ppm), 0L)
})

test_that("configuration threads through to the result", {
  set.seed(3)
  traces <- make_multi_traces(2, 14, 75)
  files <- make_temp_files(2)
  by_file <- stats::setNames(traces, files)
  result <- paramounter(files, pm_config(ppm_cutoff = 3), reader = function(f) by_file[[f]])
  expect_true(all(result@distributions$ppm < 3))
  expect_equal(result@config@ppm_cutoff, 3)
})

test_that("inputs are validated", {
  files <- make_temp_files(1)
  reader <- function(f) make_multi_traces(1, 5, 40)[[1]]
  expect_error(paramounter(42), "non-empty character")
  expect_error(paramounter(character(0)), "non-empty character")
  expect_error(paramounter(c(files, NA)), "non-empty character")
  expect_error(paramounter("/no/such/file.mzML"), "not found")
  expect_error(paramounter(files, config = "x", reader = reader), "pm_config object")
  expect_error(paramounter(files, reader = 42), "must be a function")
})

test_that("tolerance zone settings thread through to the translated parameters", {
  set.seed(505)
  traces <- make_multi_traces(3, 16, 80)
  files <- make_temp_files(3)
  by_file <- stats::setNames(traces, files)
  reader <- function(f) by_file[[f]]

  all <- paramounter(files, pm_config(), reader = reader)
  off <- paramounter(files, pm_config(tolerance_isolated = FALSE, tolerance_min_scans = 1L), reader = reader)
  on <- paramounter(files, pm_config(tolerance_isolated = TRUE, tolerance_min_scans = 5L), reader = reader)

  # explicit defaults are the defaults
  for (q in names(all@distributions)) {
    expect_equal(off@distributions[[q]], all@distributions[[q]])
  }
  # restricted tolerance distributions are a non-empty subset of the full ones;
  # threshold and shift distributions are untouched
  for (q in c("ppm", "mz_diff", "width_seconds")) {
    expect_gt(length(on@distributions[[q]]), 0L)
    expect_lt(length(on@distributions[[q]]), length(all@distributions[[q]]))
    expect_true(all(on@distributions[[q]] %in% all@distributions[[q]]))
  }
  for (q in c("noise", "width_scans", "sn", "height", "mass_shift", "rt_shift")) {
    expect_equal(on@distributions[[q]], all@distributions[[q]])
  }
  # and the translator reads the restricted distributions
  est <- point_estimates(on@distributions, FALSE, on@config@tolerance_quantile)
  v <- xcms_values(on)
  expect_equal(v$ppm, est$max_ppm)
  expect_equal(v$peakwidth, est$peakwidth)
})
