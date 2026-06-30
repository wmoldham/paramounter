# test-collect_zoi_masses.R

empty_lists <- function(n) replicate(n, numeric(0), simplify = FALSE)

test_that("an apex with no valid neighbours returns only the reference", {
  mz <- empty_lists(5)
  int <- empty_lists(5)
  mz[[3]] <- 220.01
  int[[3]] <- 5000
  r <- collect_zoi_masses(mz, int, c(0, 0, 5000, 0, 0), cutoff = 100, apex = 3, reference_mz = 220.01)
  expect_identical(r$masses, 220.01)
  expect_equal(r$left_edge, 3L)
  expect_equal(r$right_edge, 3L)
})

test_that("a symmetric single peak collects masses apex-out and spans the peak", {
  mz <- empty_lists(7)
  int <- empty_lists(7)
  mzv <- c(NA, 220.0102, 220.0098, 220.0100, 220.0101, 220.0103, NA)
  intv <- c(0, 1000, 3000, 5000, 3000, 1000, 0)
  for (s in 2:6) {
    mz[[s]] <- mzv[s]
    int[[s]] <- intv[s]
  }
  r <- collect_zoi_masses(mz, int, intv, cutoff = 500, apex = 4, reference_mz = 220.0100)
  # order is reference, then left-going, then right-going
  expect_identical(r$masses, c(220.0100, 220.0098, 220.0102, 220.0101, 220.0103))
  expect_equal(r$left_edge, 2L)
  expect_equal(r$right_edge, 6L)
})

test_that("the nearest m/z to the reference is taken when a scan has several", {
  mz <- empty_lists(5)
  int <- empty_lists(5)
  mz[[2]] <- c(220.0500, 220.0099)
  int[[2]] <- c(800, 900)
  mz[[3]] <- 220.0100
  int[[3]] <- 5000
  r <- collect_zoi_masses(mz, int, c(0, 850, 5000, 0, 0), cutoff = 100, apex = 3, reference_mz = 220.0100)
  expect_identical(r$masses, c(220.0100, 220.0099))
})

test_that("the walk stops where mean intensity drops below the cutoff", {
  mz <- empty_lists(5)
  int <- empty_lists(5)
  mz[[2]] <- 220.0102
  int[[2]] <- 300 # below cutoff 500
  mz[[3]] <- 220.0099
  int[[3]] <- 4000
  mz[[4]] <- 220.0100
  int[[4]] <- 5000
  r <- collect_zoi_masses(mz, int, c(0, 300, 4000, 5000, 0), cutoff = 500, apex = 4, reference_mz = 220.0100)
  expect_identical(r$masses, c(220.0100, 220.0099))
  expect_equal(r$left_edge, 3L)
  expect_equal(r$right_edge, 4L)
})

test_that("an apex at the first scan does not read out of bounds", {
  mz <- empty_lists(4)
  int <- empty_lists(4)
  mz[[1]] <- 220.0100
  int[[1]] <- 5000
  mz[[2]] <- 220.0101
  int[[2]] <- 2000
  r <- collect_zoi_masses(mz, int, c(5000, 2000, 0, 0), cutoff = 500, apex = 1, reference_mz = 220.0100)
  expect_identical(r$masses, c(220.0100, 220.0101))
  expect_equal(r$left_edge, 1L)
  expect_equal(r$right_edge, 2L)
})

test_that("a mass-shifted adjacent peak is cut off by the outlier test", {
  n <- 12
  mz <- empty_lists(n)
  int <- empty_lists(n)
  mz_first <- c(220.0098, 220.0099, 220.0100, 220.0101, 220.0102, 220.0103, 220.0104)
  int_first <- c(1000, 2000, 4000, 8000, 4000, 2000, 1000) # apex at scan 5
  for (i in 1:7) {
    mz[[i + 1]] <- mz_first[i]
    int[[i + 1]] <- int_first[i]
  }
  mz_second <- c(220.0500, 220.0501, 220.0502, 220.0503)
  int_second <- c(1500, 3000, 5000, 3000) # rising second peak
  for (i in 1:4) {
    mz[[i + 8]] <- mz_second[i]
    int[[i + 8]] <- int_second[i]
  }
  eic <- vapply(int, function(v) if (length(v)) mean(v) else 0, numeric(1))
  r <- collect_zoi_masses(mz, int, eic, cutoff = 500, apex = 5, reference_mz = 220.0102)
  expect_length(r$masses, 8L)
  expect_equal(r$masses[8], 220.0500) # outlier appended, then the walk breaks
  expect_equal(r$left_edge, 2L)
  expect_equal(r$right_edge, 8L)
  expect_false(any(r$masses > 220.0501)) # later scans never reached
})

test_that("the rewrite matches the original walk across random bins", {
  orig_walk <- function(mzL, intL, eic, cut, apex, ref) {
    m <- ref
    L <- apex - 1
    R <- apex + 1
    n <- length(mzL)
    tk <- function(v) {
      q1 <- as.numeric(summary(v)[2])
      q3 <- as.numeric(summary(v)[5])
      last <- v[length(v)]
      last < q1 - 1.5 * (q3 - q1) || last > q3 + 1.5 * (q3 - q1)
    }
    if (L > 0) {
      while (length(mzL[[L]]) > 0 & mean(intL[[L]]) >= cut) {
        if (length(mzL[[L]]) == 1) {
          m <- c(m, mzL[[L]])
        } else {
          ab <- abs(mzL[[L]] - ref)
          m <- c(m, mzL[[L]][which(ab == min(ab))[1]])
        }
        if (eic[L] > eic[L + 1] & length(m) > 5 && tk(m)) break
        L <- L - 1
        if (L <= 0) break
      }
    }
    if (R <= n) {
      while (length(mzL[[R]]) > 0 & mean(intL[[R]]) >= cut) {
        if (length(mzL[[R]]) == 1) {
          m <- c(m, mzL[[R]])
        } else {
          ab <- abs(mzL[[R]] - ref)
          m <- c(m, mzL[[R]][which(ab == min(ab))[1]])
        }
        if (eic[R] > eic[R - 1] & length(m) > 5 && tk(m)) break
        R <- R + 1
        if (R > n) break
      }
    }
    list(masses = m, left_edge = L + 1, right_edge = R - 1)
  }
  make_bin <- function(n, peaks, center) {
    mzL <- replicate(n, numeric(0), simplify = FALSE)
    intL <- replicate(n, numeric(0), simplify = FALSE)
    for (pk in peaks) {
      for (s in 1:n) {
        g <- pk$h * exp(-0.5 * ((s - pk$a) / pk$s)^2)
        if (g > pk$h * 0.02) {
          mzL[[s]] <- c(mzL[[s]], center + pk$off + rnorm(1, 0, center * 2e-6))
          intL[[s]] <- c(intL[[s]], g)
        }
      }
    }
    for (s in 1:n) {
      if (length(intL[[s]]) > 0 && runif(1) < 0.12) {
        mzL[[s]] <- c(mzL[[s]], center + rnorm(1, 0, center * 8e-6))
        intL[[s]] <- c(intL[[s]], runif(1, 0, max(intL[[s]])))
      }
    }
    eic <- vapply(intL, function(v) if (length(v)) mean(v) else 0, numeric(1))
    list(mzL = mzL, intL = intL, eic = eic)
  }
  set.seed(21)
  for (rep in 1:400) {
    n <- sample(20:120, 1)
    center <- sample(c(120, 220, 450, 800), 1)
    two <- runif(1) < 0.5
    peaks <- if (two) {
      a1 <- sample(8:(n %/% 2), 1)
      a2 <- a1 + sample(4:12, 1)
      if (a2 > n - 8) next
      list(
        list(a = a1, s = runif(1, 2, 5), h = runif(1, 1e4, 1e5), off = 0),
        list(a = a2, s = runif(1, 2, 5), h = runif(1, 1e4, 1e5), off = runif(1, 5e-3, 2e-2))
      )
    } else {
      list(list(a = sample(8:(n - 8), 1), s = runif(1, 2, 6), h = runif(1, 1e3, 1e5), off = 0))
    }
    b <- make_bin(n, peaks, center)
    if (max(b$eic) == 0) next
    apex <- which.max(b$eic)
    ref <- b$mzL[[apex]][which.max(b$intL[[apex]])]
    cut <- max(b$eic) * runif(1, 0.02, 0.4)
    a <- orig_walk(b$mzL, b$intL, b$eic, cut, apex, ref)
    z <- collect_zoi_masses(b$mzL, b$intL, b$eic, cut, apex, ref)
    expect_equal(z$masses, a$masses)
    expect_equal(z$left_edge, a$left_edge)
    expect_equal(z$right_edge, a$right_edge)
  }
})

test_that("structural inputs are validated", {
  mz <- empty_lists(3)
  int <- empty_lists(3)
  eic <- c(1, 2, 3)
  expect_error(collect_zoi_masses("x", int, eic, 1, 2, 220), "mz_list")
  expect_error(collect_zoi_masses(mz, "x", eic, 1, 2, 220), "int_list")
  expect_error(collect_zoi_masses(mz, int, "x", 1, 2, 220), "numeric")
  expect_error(collect_zoi_masses(mz, int, c(1, NA, 3), 1, 2, 220), "NA")
  expect_error(collect_zoi_masses(mz, int, c(1, 2), 1, 2, 220), "same length")
  expect_error(collect_zoi_masses(empty_lists(2), int, eic, 1, 2, 220), "same length")
  expect_error(
    collect_zoi_masses(list(1, c(1, 2), 3), list(1, 1, 1), eic, 1, 2, 220),
    "matching length"
  )
})

test_that("scalar inputs are validated", {
  mz <- empty_lists(3)
  int <- empty_lists(3)
  eic <- c(1, 2, 3)
  expect_error(collect_zoi_masses(mz, int, eic, -1, 2, 220), "non-negative")
  expect_error(collect_zoi_masses(mz, int, eic, c(1, 2), 2, 220), "single")
  expect_error(collect_zoi_masses(mz, int, eic, 1, 0, 220), "within")
  expect_error(collect_zoi_masses(mz, int, eic, 1, 4, 220), "within")
  expect_error(collect_zoi_masses(mz, int, eic, 1, 2.5, 220), "integer")
  expect_error(collect_zoi_masses(mz, int, eic, 1, 2, NA_real_), "finite")
  expect_error(collect_zoi_masses(mz, int, eic, 1, 2, c(220, 221)), "single")
  expect_error(collect_zoi_masses(mz, int, eic, 1, 2, 220, min_points = -1), "non-negative")
  expect_error(collect_zoi_masses(mz, int, eic, 1, 2, 220, min_points = 2.5), "integer")
  expect_error(collect_zoi_masses(mz, int, eic, 1, 2, 220, fence = -1), "non-negative")
  expect_error(collect_zoi_masses(mz, int, eic, 1, 2, 220, fence = Inf), "non-negative")
})
