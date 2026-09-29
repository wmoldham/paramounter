# test-aggregate_files.R

test_that("trim_distribution keeps the requested tail", {
  x <- 1:100
  expect_equal(trim_distribution(x, 0.97, "high"), as.numeric(1:97))
  expect_equal(trim_distribution(x, 0.97, "low"), as.numeric(100:4))
  expect_length(trim_distribution(numeric(0)), 0L)
  expect_equal(trim_distribution(c(5, 1, 3, NA, 2, 4), 1, "high"), c(1, 2, 3, 4, 5))
})

test_that("compute_ppm_cutoff takes the quantile of the trimmed distribution", {
  expect_equal(compute_ppm_cutoff(1:1000, 0.97, 0.95), 922)
  expect_equal(compute_ppm_cutoff(numeric(0), 0.97, 0.95), Inf)
})

test_that("aggregate_files reproduces part 1 + part 2 across multi-file data", {
  orig_match <- function(zl, mt = 0.015, rt_ = 30) {
    nf <- length(zl); f1 <- zl[[1]]; f2 <- zl[[2]]; al <- list()
    for (i in seq_len(nrow(f1))) { tp <- f2[f2$mz >= f1$mz[i] - mt & f2$mz <= f1$mz[i] + mt & f2$rt >= f1$rt[i] - rt_ & f2$rt <= f1$rt[i] + rt_, , drop = FALSE]
    if (nrow(tp) > 0) al[[length(al) + 1]] <- list(mz = c(f1$mz[i], tp$mz[1]), rt = c(f1$rt[i], tp$rt[1])) }
    if (nf >= 3) for (k in 3:nf) { fk <- zl[[k]]; na <- list()
    for (a in al) { tp <- fk[fk$mz >= a$mz[1] - mt & fk$mz <= a$mz[1] + mt & fk$rt >= a$rt[1] - rt_ & fk$rt <= a$rt[1] + rt_, , drop = FALSE]
    if (nrow(tp) > 0) na[[length(na) + 1]] <- list(mz = c(a$mz, tp$mz[1]), rt = c(a$rt, tp$rt[1])) }
    al <- na }
    if (length(al) == 0) return(list(mz = matrix(numeric(0), 0, nf), rt = matrix(numeric(0), 0, nf)))
    list(mz = do.call(rbind, lapply(al, `[[`, "mz")), rt = do.call(rbind, lapply(al, `[[`, "rt")))
  }
  orig_aggregate <- function(per_file, trim = 0.97, quantile = 0.95, mz_tol = 0.015, rt_tol = 30) {
    all_ppm <- unlist(lapply(per_file, function(x) x$zoi$ppm))
    s <- sort(all_ppm); nt <- round(length(s) * trim); tr <- s[1:nt]; cutoff <- tr[round(length(tr) * quantile)]
    ppm <- c(); mzd <- c(); pw <- c(); ps <- c(); sn <- c(); ht <- c(); noise <- c(); clean <- list()
    for (f in seq_along(per_file)) { z <- per_file[[f]]$zoi; noise <- c(noise, per_file[[f]]$noise); keep <- z$ppm < cutoff
    ppm <- c(ppm, z$ppm[keep]); mzd <- c(mzd, z$mz_diff[keep]); pw <- c(pw, z$width_seconds[keep]); ps <- c(ps, z$width_scans[keep])
    sn <- c(sn, z$sn[keep]); ht <- c(ht, z$height[keep]); ck <- keep & z$isolated
    clean[[f]] <- data.frame(mz = z$reference_mz[ck], rt = z$apex_rt[ck]) }
    noise <- sort(noise)[1:round(length(noise) * trim)]; mzd <- sort(mzd)[1:round(length(mzd) * trim)]
    pw <- sort(pw)[1:round(length(pw) * trim)]; ps <- sort(ps)[1:round(length(ps) * trim)]
    sn <- sort(sn, decreasing = TRUE)[1:round(length(sn) * trim)]; ht <- sort(ht, decreasing = TRUE)[1:round(length(ht) * trim)]
    ms <- numeric(0); rs <- numeric(0)
    if (length(per_file) > 1) { matched <- orig_match(clean, mz_tol, rt_tol)
    if (nrow(matched$mz) > 0) { ms <- apply(matched$mz, 1, function(r) max(r) - min(r)); rs <- apply(matched$rt, 1, function(r) max(r) - min(r))
    ms <- sort(ms)[1:round(length(ms) * trim)]; rs <- sort(rs)[1:round(length(rs) * trim)] } }
    list(ppm = ppm, mz_diff = mzd, noise = noise, width_seconds = pw, width_scans = ps, sn = sn, height = ht, mass_shift = ms, rt_shift = rs)
  }
  make_multi <- function(n_files, n_shared, n_scans, mz_lo = 100, mz_hi = 110) {
    smz <- sort(runif(n_shared, mz_lo + 1, mz_hi - 1)); srt <- runif(n_shared, 10, n_scans - 10); sw <- runif(n_shared, 2.5, 6); shh <- runif(n_shared, 1e4, 1e6)
    lapply(seq_len(n_files), function(fi) {
      dmz <- rnorm(1, 0, 0.001); drt <- rnorm(1, 0, 0.5); mzD <- vector("list", n_scans); intD <- vector("list", n_scans)
      for (s in seq_len(n_scans)) { pm <- numeric(0); pv <- numeric(0)
      for (k in seq_len(n_shared)) { g <- shh[k] * exp(-0.5 * ((s - (srt[k] + drt)) / sw[k])^2)
      if (g > shh[k] * 0.02) { pm <- c(pm, smz[k] + dmz + rnorm(1, 0, smz[k] * 2e-6)); pv <- c(pv, g) } }
      nn <- rpois(1, 2); pm <- c(pm, runif(nn, mz_lo, mz_hi)); pv <- c(pv, runif(nn, 1e2, 1.5e3)); mzD[[s]] <- pm; intD[[s]] <- pv }
      mzD[[1]] <- c(mzD[[1]], mz_lo + 0.007); intD[[1]] <- c(intD[[1]], 300); mzD[[2]] <- c(mzD[[2]], mz_hi - 0.007); intD[[2]] <- c(intD[[2]], 300)
      for (s in seq_len(n_scans)) { o <- order(mzD[[s]]); mzD[[s]] <- mzD[[s]][o]; intD[[s]] <- intD[[s]][o] }
      list(mz = mzD, int = intD, rtime = seq(0, by = 6, length.out = n_scans)) })
  }
  # orig_match takes the first candidate in the window (temp$mz[1]) and
  # orig_aggregate consumes per_file verbatim, so both the measurement and the
  # aggregation side must be pinned to legacy rather than relying on defaults.
  set.seed(202)
  legacy <- pm_config(legacy = TRUE)
  for (rep in 1:12) {
    files <- make_multi(sample(2:4, 1), sample(10:22, 1), sample(60:100, 1))
    per_file <- lapply(files, function(f) measure_file(f$mz, f$int, f$rtime, legacy))
    fn <- paste0("f", seq_along(files), ".mzML")
    up <- aggregate_files(per_file, fn, legacy)
    a <- orig_aggregate(per_file)
    for (q in names(a)) expect_equal(sort(up@distributions[[q]]), sort(as.numeric(a[[q]])))
  }
})

test_that("aggregate_files honours a manual cutoff and handles a single file", {
  make_one <- function(n_scans, mz_lo = 100, mz_hi = 106) {
    smz <- sort(runif(10, mz_lo + 1, mz_hi - 1)); srt <- runif(10, 10, n_scans - 10); sw <- runif(10, 3, 6); shh <- runif(10, 1e4, 1e6)
    mzD <- vector("list", n_scans); intD <- vector("list", n_scans)
    for (s in seq_len(n_scans)) { pm <- numeric(0); pv <- numeric(0)
    for (k in 1:10) { g <- shh[k] * exp(-0.5 * ((s - srt[k]) / sw[k])^2); if (g > shh[k] * 0.02) { pm <- c(pm, smz[k] + rnorm(1, 0, smz[k] * 2e-6)); pv <- c(pv, g) } }
    mzD[[s]] <- pm; intD[[s]] <- pv }
    mzD[[1]] <- c(mzD[[1]], mz_lo + 0.007); intD[[1]] <- c(intD[[1]], 300); mzD[[2]] <- c(mzD[[2]], mz_hi - 0.007); intD[[2]] <- c(intD[[2]], 300)
    for (s in seq_len(n_scans)) { o <- order(mzD[[s]]); mzD[[s]] <- mzD[[s]][o]; intD[[s]] <- intD[[s]][o] }
    list(mz = mzD, int = intD, rtime = seq(0, by = 6, length.out = n_scans))
  }
  set.seed(7)
  f <- make_one(80)
  pf <- measure_file(f$mz, f$int, f$rtime)
  up_manual <- aggregate_files(list(pf), "a.mzML", pm_config(ppm_cutoff = 3))
  expect_true(all(up_manual@distributions$ppm < 3))
  up1 <- aggregate_files(list(pf), "a.mzML")
  expect_length(up1@distributions$mass_shift, 0L)
  expect_length(up1@distributions$rt_shift, 0L)
  expect_true(S7::S7_inherits(up1, universal_parameters))
})

test_that("aggregate_files validates its inputs", {
  pf <- list(noise = c(1, 2), zoi = data.frame(ppm = 1, mz_diff = 1, width_seconds = 1, width_scans = 1L, sn = 1, height = 1, reference_mz = 200, apex_rt = 1, isolated = TRUE))
  expect_error(aggregate_files(list(), "a"), "non-empty list")
  expect_error(aggregate_files(list(list(noise = 1)), "a"), "numeric `noise` and a data frame")
  expect_error(aggregate_files(list(pf), c("a", "b")), "one per file")
  expect_error(aggregate_files(list(pf), "a", config = "x"), "pm_config object")
})

test_that("tolerance zone settings restrict only the tolerance distributions", {
  zoi <- data.frame(
    ppm = c(1, 2, 40, 60, 3), mz_diff = c(0.001, 0.002, 0.02, 0.03, 0.003),
    width_seconds = c(10, 20, 1, 1.2, 30), width_scans = c(8, 15, 2, 2, 25),
    sn = c(50, 80, 4, 5, 90), height = c(1e5, 2e5, 1e4, 2e4, 3e5),
    reference_mz = 100 + 1:5, apex_rt = c(60, 400, 120, 700, 900),
    isolated = c(TRUE, TRUE, TRUE, FALSE, FALSE)
  )
  pf <- list(list(noise = c(100, 200, 300), zoi = zoi))
  all <- aggregate_files(pf, "f.mzML", pm_config(ppm_cutoff = 100))
  sel <- aggregate_files(
    pf, "f.mzML",
    pm_config(ppm_cutoff = 100, tolerance_isolated = TRUE, tolerance_min_scans = 5L)
  )
  # rows 1 and 2 are isolated and wide enough; row 3 is isolated but a 2-scan
  # blip, which is the case isolation alone does not catch
  expect_equal(sort(sel@distributions$ppm), c(1, 2))
  expect_equal(sort(sel@distributions$mz_diff), c(0.001, 0.002))
  expect_equal(sort(sel@distributions$width_seconds), c(10, 20))
  expect_equal(sort(all@distributions$ppm), sort(zoi$ppm))
  for (q in c("noise", "width_scans", "sn", "height")) {
    expect_equal(sel@distributions[[q]], all@distributions[[q]])
  }
  only_wide <- aggregate_files(pf, "f.mzML", pm_config(ppm_cutoff = 100, tolerance_min_scans = 5L))
  expect_equal(sort(only_wide@distributions$ppm), c(1, 2, 3))
  # isolation alone keeps the isolated 2-scan blip (row 3)
  only_iso <- aggregate_files(pf, "f.mzML", pm_config(ppm_cutoff = 100, tolerance_isolated = TRUE))
  expect_equal(sort(only_iso@distributions$ppm), c(1, 2, 40))
  # the restriction applies after the ppm cutoff, not instead of it
  cut <- aggregate_files(
    pf, "f.mzML",
    pm_config(ppm_cutoff = 1.5, tolerance_isolated = TRUE, tolerance_min_scans = 5L)
  )
  expect_equal(cut@distributions$ppm, 1)
})

test_that("tolerance zone settings are ignored under legacy and fail loudly when empty", {
  zoi <- data.frame(
    ppm = c(1, 2), mz_diff = c(0.001, 0.002), width_seconds = c(1, 2),
    width_scans = c(2, 2), sn = c(5, 6), height = c(1e4, 2e4),
    reference_mz = c(101, 102), apex_rt = c(10, 20), isolated = c(FALSE, FALSE)
  )
  pf <- list(list(noise = c(100, 200), zoi = zoi))
  expect_error(
    aggregate_files(pf, "f.mzML", pm_config(ppm_cutoff = 100, tolerance_isolated = TRUE)),
    "No zones left"
  )
  leg <- aggregate_files(
    pf, "f.mzML",
    pm_config(ppm_cutoff = 100, legacy = TRUE, tolerance_isolated = TRUE, tolerance_min_scans = 5L)
  )
  expect_length(leg@distributions$ppm, 2L)
})
