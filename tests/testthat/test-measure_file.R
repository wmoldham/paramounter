# test-measure_file.R

test_that("measure_file reproduces the original inner loop across synthetic files", {
  peak_smooth <- function(x, level) {
    n <- level
    if (length(x) < 2 * n) {
      return(x)
    } else if (length(unique(x)) == 1) {
      return(x)
    } else {
      y <- vector(length = length(x))
      for (i in 1:n) y[i] <- sum(c((n - i + 2):(n + 1), n:1) * x[1:(i + n)]) / sum(c((n - i + 2):(n + 1), n:1))
      for (i in (n + 1):(length(y) - n)) y[i] <- sum(c(1:(n + 1), n:1) * x[(i - n):(i + n)]) / sum(c(1:(n + 1), n:1))
      for (i in (length(y) - n + 1):length(y)) y[i] <- sum(c(1:n, (n + 1):(n + i - length(x) + 1)) * x[(i - n):length(x)]) / sum(c(1:n, (n + 1):(n + i - length(x) + 1)))
      return(y)
    }
  }
  orig_measure_file <- function(mzData, intData, rtime, binwidth = 0.05, smooth = 0, massSDrange = 2) {
    allmz <- unlist(mzData)
    ROI <- seq(min(allmz), max(allmz), binwidth)
    noiseALL <- c(); refMzV <- c(); apexRtV <- c(); ppmV <- c(); daV <- c()
    widthV <- c(); scansV <- c(); snV <- c(); htV <- c(); isoV <- c()
    tk <- function(v) { Q1 <- as.numeric(summary(v)[2]); Q3 <- as.numeric(summary(v)[5]); l <- v[length(v)]; l < Q1 - 1.5 * (Q3 - Q1) || l > Q3 + 1.5 * (Q3 - Q1) }
    for (i in seq_len(length(ROI) - 1)) {
      lo <- ROI[i]; hi <- ROI[i + 1]
      tmpMZ <- vector("list", length(mzData)); tmpINT <- vector("list", length(mzData))
      for (j in seq_along(mzData)) { idx <- which(mzData[[j]] >= lo & mzData[[j]] < hi); tmpMZ[[j]] <- mzData[[j]][idx]; tmpINT[[j]] <- intData[[j]][idx] }
      eicINTraw <- numeric(length(mzData))
      for (k in seq_along(mzData)) eicINTraw[k] <- if (length(tmpINT[[k]]) > 0) mean(tmpINT[[k]]) else 0
      if (sum(eicINTraw != 0) == 0) next
      eicINT <- peak_smooth(eicINTraw, smooth); eicNon0 <- sort(eicINT[eicINT > 0]); blk <- 0; sdv <- 0
      if (length(eicNon0) > 10) {
        for (x in seq(10, length(eicNon0), 10)) { sdv <- sd(eicNon0[1:x]); blk <- sum(eicNon0[1:x]) / x; thres <- blk + 3 * sdv
        if (x + 1 <= length(eicNon0)) { if (eicNon0[x + 1] >= thres) break } }
        cutOFF <- eicNon0[x]
      } else { cutOFF <- max(eicNon0); sdv <- 0; blk <- 0 }
      noiseALL <- c(noiseALL, cutOFF)
      aboveTHindex <- which(eicINT > cutOFF); if (length(aboveTHindex) == 0) next
      seg <- split(aboveTHindex, cumsum(c(1, diff(aboveTHindex) != 1))); peakInd <- c()
      for (x in seq_along(seg)) peakInd[x] <- which(eicINT[seg[[x]]] == max(eicINT[seg[[x]]]))[1] + min(seg[[x]]) - 1
      refMZvec <- c()
      for (y in seq_along(peakInd)) { h2 <- which(tmpINT[[peakInd[y]]] == max(tmpINT[[peakInd[y]]]))[1]; refMZvec[y] <- tmpMZ[[peakInd[y]]][h2] }
      for (z in seq_along(peakInd)) {
        currPeakInd <- peakInd[z]; currRefMz <- refMZvec[z]; csm <- currRefMz; leftInd <- currPeakInd - 1; rightInd <- currPeakInd + 1
        if (leftInd > 0) while (length(tmpMZ[[leftInd]]) > 0 & mean(tmpINT[[leftInd]]) >= cutOFF) {
          if (length(tmpMZ[[leftInd]]) == 1) csm <- c(csm, tmpMZ[[leftInd]]) else { ab <- abs(tmpMZ[[leftInd]] - currRefMz); csm <- c(csm, tmpMZ[[leftInd]][which(ab == min(ab))[1]]) }
          if (eicINT[leftInd] > eicINT[leftInd + 1] & length(csm) > 5 && tk(csm)) break
          leftInd <- leftInd - 1; if (leftInd <= 0) break
        }
        if (rightInd <= length(tmpMZ)) while (length(tmpMZ[[rightInd]]) > 0 & mean(tmpINT[[rightInd]]) >= cutOFF) {
          if (length(tmpMZ[[rightInd]]) == 1) csm <- c(csm, tmpMZ[[rightInd]]) else { ab <- abs(tmpMZ[[rightInd]] - currRefMz); csm <- c(csm, tmpMZ[[rightInd]][which(ab == min(ab))[1]]) }
          if (eicINT[rightInd] > eicINT[rightInd - 1] & length(csm) > 5 && tk(csm)) break
          rightInd <- rightInd + 1; if (rightInd > length(tmpMZ)) break
        }
        if (sum(is.na(csm)) > 0) next
        if (length(csm) > 1 && length(csm) < 200) {
          cpw <- rtime[[rightInd - 1]] - rtime[[leftInd + 1]]; cps <- rightInd - leftInd
          in_shift <- FALSE
          if (z == 1) { if (is.na(peakInd[z + 1])) in_shift <- TRUE
          if (!is.na(peakInd[z + 1])) { if (rtime[[currPeakInd]] < (rtime[[peakInd[z + 1]]] - 300)) in_shift <- TRUE } }
          if (z > 1) { if (is.na(peakInd[z + 1]) & rtime[[currPeakInd]] > (rtime[[peakInd[z - 1]]] + 300)) in_shift <- TRUE
          if (!is.na(peakInd[z + 1])) { if (rtime[[currPeakInd]] < (rtime[[peakInd[z + 1]]] - 300)) in_shift <- TRUE } }
          refMzV <- c(refMzV, currRefMz); apexRtV <- c(apexRtV, rtime[[currPeakInd]]); ppmV <- c(ppmV, massSDrange * sd(csm) / currRefMz * 1e6)
          daV <- c(daV, massSDrange * sd(csm)); widthV <- c(widthV, cpw); scansV <- c(scansV, cps)
          snV <- c(snV, if (cutOFF > 0) (eicINT[currPeakInd] - blk) / sdv else eicINT[currPeakInd]); htV <- c(htV, eicINT[currPeakInd]); isoV <- c(isoV, in_shift)
        }
      }
    }
    list(noise = noiseALL, zoi = data.frame(ppm = ppmV, mz_diff = daV, width_seconds = widthV, width_scans = scansV, sn = snV, height = htV, reference_mz = refMzV, apex_rt = apexRtV, isolated = isoV))
  }
  make_lcms <- function(n_scans, n_feat, mz_lo = 100, mz_hi = 105) {
    fmz <- runif(n_feat, mz_lo, mz_hi); frt <- runif(n_feat, 1, n_scans); fw <- runif(n_feat, 2, 7); fh <- runif(n_feat, 5e3, 1e6)
    nc <- max(1, n_feat %/% 5)
    fmz <- c(fmz, fmz[seq_len(nc)] + 0.001); frt <- c(frt, (frt[seq_len(nc)] + n_scans / 2) %% n_scans + 1)
    fw <- c(fw, runif(nc, 2, 7)); fh <- c(fh, runif(nc, 5e3, 1e6)); K <- length(fmz)
    mzD <- vector("list", n_scans); intD <- vector("list", n_scans)
    for (s in seq_len(n_scans)) {
      pm <- numeric(0); pv <- numeric(0)
      for (f in seq_len(K)) { g <- fh[f] * exp(-0.5 * ((s - frt[f]) / fw[f])^2)
      if (g > fh[f] * 0.02) { pm <- c(pm, fmz[f] + rnorm(1, 0, fmz[f] * 3e-6)); pv <- c(pv, g) } }
      nn <- rpois(1, 3); pm <- c(pm, runif(nn, mz_lo, mz_hi)); pv <- c(pv, runif(nn, 1e2, 2e3))
      o <- order(pm); mzD[[s]] <- pm[o]; intD[[s]] <- pv[o]
    }
    list(mz = mzD, int = intD, rtime = seq(0, by = 6, length.out = n_scans))
  }
  set.seed(101)
  for (rep in 1:15) {
    sp <- make_lcms(sample(50:90, 1), sample(8:20, 1))
    a <- orig_measure_file(sp$mz, sp$int, sp$rtime)
    b <- measure_file(sp$mz, sp$int, sp$rtime)
    expect_equal(b$noise, a$noise)
    expect_equal(nrow(b$zoi), nrow(a$zoi))
    if (nrow(a$zoi) > 0) {
      for (col in names(a$zoi)) expect_equal(b$zoi[[col]], a$zoi[[col]])
    }
  }
})

test_that("measure_file returns the documented structure", {
  n <- 20; mz <- vector("list", n); intn <- vector("list", n)
  set.seed(5)
  for (s in seq_len(n)) {
    pm <- numeric(0); pv <- numeric(0)
    for (fc in c(200.123, 200.377)) { g <- 1e5 * exp(-0.5 * ((s - 10) / 3)^2)
    if (g > 1e3) { pm <- c(pm, fc + rnorm(1, 0, fc * 2e-6)); pv <- c(pv, g) } }
    mz[[s]] <- pm; intn[[s]] <- pv
  }
  mz[[1]] <- c(mz[[1]], 200.001); intn[[1]] <- c(intn[[1]], 500)
  mz[[2]] <- c(mz[[2]], 200.548); intn[[2]] <- c(intn[[2]], 500)
  for (s in seq_len(n)) { o <- order(mz[[s]]); mz[[s]] <- mz[[s]][o]; intn[[s]] <- intn[[s]][o] }
  res <- measure_file(mz, intn, seq(0, by = 6, length.out = n))
  expect_named(res, c("noise", "zoi"))
  expect_type(res$noise, "double")
  expect_s3_class(res$zoi, "data.frame")
  expect_named(res$zoi, c("ppm", "mz_diff", "width_seconds", "width_scans", "sn", "height", "reference_mz", "apex_rt", "isolated"))
  expect_equal(nrow(res$zoi), 2L)
  expect_true(all(res$zoi$isolated))
})

test_that("empty input gives empty results", {
  res <- measure_file(list(numeric(0), numeric(0)), list(numeric(0), numeric(0)), c(0, 6))
  expect_length(res$noise, 0L)
  expect_equal(nrow(res$zoi), 0L)
})

test_that("inputs are validated", {
  expect_error(measure_file(42, list(), 1), "must be lists")
  expect_error(measure_file(list(1), list(1, 2), 1), "same length")
  expect_error(measure_file(list(1), list(1), c(1, NA)), "NA")
  expect_error(measure_file(list(1), list(1), "x"), "numeric")
  expect_error(measure_file(list(1, 2), list(1, 2), 1), "one value per scan")
  expect_error(measure_file(list(1), list(1), 1, config = "nope"), "pm_config object")
})
