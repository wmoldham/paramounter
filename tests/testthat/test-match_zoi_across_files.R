# test-match_zoi_across_files.R

df <- function(mz, rt) data.frame(mz = mz, rt = rt)

test_that("features shared by two files are matched", {
  r <- match_zoi_across_files(list(df(c(200, 300), c(500, 800)), df(c(200.003, 300.002), c(505, 803))))
  expect_equal(r$mz, rbind(c(200, 200.003), c(300, 300.002)))
  expect_equal(r$rt, rbind(c(500, 505), c(800, 803)))
})

test_that("a feature absent from any file is dropped", {
  r <- match_zoi_across_files(list(df(c(200, 300), c(500, 800)), df(200.003, 505)))
  expect_equal(nrow(r$mz), 1L)
  expect_equal(r$mz, matrix(c(200, 200.003), 1))
  expect_equal(r$rt, matrix(c(500, 505), 1))
})

test_that("a crowded window takes the first candidate in legacy mode and the nearest otherwise", {
  f1 <- df(200.000, 500)
  f2 <- df(c(200.010, 200.002), c(520, 505)) # far candidate listed before near candidate
  legacy <- match_zoi_across_files(list(f1, f2), legacy_first_match = TRUE)
  nearest <- match_zoi_across_files(list(f1, f2), legacy_first_match = FALSE)
  expect_equal(legacy$mz, matrix(c(200.000, 200.010), 1))
  expect_equal(nearest$mz, matrix(c(200.000, 200.002), 1))
})

test_that("matching extends across three files", {
  r <- match_zoi_across_files(list(df(200, 500), df(200.003, 505), df(200.005, 510)))
  expect_equal(r$mz, matrix(c(200, 200.003, 200.005), 1))
  expect_equal(r$rt, matrix(c(500, 505, 510), 1))
})

test_that("a feature must be present in every file", {
  r <- match_zoi_across_files(list(df(200, 500), df(200.003, 505), df(400, 900)))
  expect_equal(nrow(r$mz), 0L)
  expect_equal(ncol(r$mz), 3L)
})

test_that("an empty anchor yields no features", {
  r <- match_zoi_across_files(list(df(numeric(0), numeric(0)), df(200, 500)))
  expect_equal(nrow(r$mz), 0L)
  expect_equal(ncol(r$mz), 2L)
})

test_that("the tolerance windows are inclusive at the edge", {
  expect_equal(nrow(match_zoi_across_files(list(df(200.000, 500), df(200.015, 500)))$mz), 1L)
  expect_equal(nrow(match_zoi_across_files(list(df(200.000, 500), df(200.016, 500)))$mz), 0L)
})

test_that("one file's zone may match several anchor zones", {
  r <- match_zoi_across_files(list(df(c(200.000, 200.005), c(500, 500)), df(200.003, 500)))
  expect_equal(nrow(r$mz), 2L)
})

test_that("column names are taken from the list names", {
  r <- match_zoi_across_files(list(a = df(200, 500), b = df(200.003, 505)))
  expect_identical(colnames(r$mz), c("a", "b"))
  expect_identical(colnames(r$rt), c("a", "b"))
})

test_that("legacy mode matches the original algorithm across random datasets", {
  orig_match <- function(zoi_list, mz_tol = 0.015, rt_tol = 30) {
    nf <- length(zoi_list)
    f1 <- zoi_list[[1]]
    f2 <- zoi_list[[2]]
    align <- list()
    for (i in seq_len(nrow(f1))) {
      temp <- f2[
        f2$mz >= f1$mz[i] - mz_tol & f2$mz <= f1$mz[i] + mz_tol &
          f2$rt >= f1$rt[i] - rt_tol & f2$rt <= f1$rt[i] + rt_tol, ,
        drop = FALSE
      ]
      if (nrow(temp) > 0) {
        align[[length(align) + 1]] <- list(mz = c(f1$mz[i], temp$mz[1]), rt = c(f1$rt[i], temp$rt[1]))
      }
    }
    if (nf >= 3) {
      for (k in 3:nf) {
        fk <- zoi_list[[k]]
        na <- list()
        for (a in align) {
          temp <- fk[
            fk$mz >= a$mz[1] - mz_tol & fk$mz <= a$mz[1] + mz_tol &
              fk$rt >= a$rt[1] - rt_tol & fk$rt <= a$rt[1] + rt_tol, ,
            drop = FALSE
          ]
          if (nrow(temp) > 0) {
            na[[length(na) + 1]] <- list(mz = c(a$mz, temp$mz[1]), rt = c(a$rt, temp$rt[1]))
          }
        }
        align <- na
      }
    }
    if (length(align) == 0) {
      return(list(mz = matrix(numeric(0), 0, nf), rt = matrix(numeric(0), 0, nf)))
    }
    list(
      mz = do.call(rbind, lapply(align, `[[`, "mz")),
      rt = do.call(rbind, lapply(align, `[[`, "rt"))
    )
  }
  make_files <- function(n_files, n_true, n_noise) {
    true_mz <- sort(runif(n_true, 100, 900))
    true_rt <- runif(n_true, 30, 1500)
    present <- matrix(runif(n_true * n_files) < 0.8, n_true, n_files)
    present[, 1] <- TRUE
    lapply(seq_len(n_files), function(f) {
      d <- data.frame(
        mz = c(true_mz[present[, f]] + rnorm(sum(present[, f]), 0, 0.002), runif(n_noise, 100, 900)),
        rt = c(true_rt[present[, f]] + rnorm(sum(present[, f]), 0, 8), runif(n_noise, 30, 1500))
      )
      d[sample(nrow(d)), , drop = FALSE]
    })
  }
  set.seed(41)
  for (rep in 1:300) {
    fl <- make_files(sample(2:5, 1), sample(5:40, 1), sample(0:15, 1))
    # orig_match takes temp$mz[1], the first candidate in the window
    a <- orig_match(fl)
    b <- match_zoi_across_files(fl, legacy_first_match = TRUE)
    expect_equal(unname(b$mz), unname(a$mz))
    expect_equal(unname(b$rt), unname(a$rt))
  }
})

test_that("the list structure is validated", {
  expect_error(match_zoi_across_files(df(1, 1)), "list of data frames")
  expect_error(match_zoi_across_files(list(df(1, 1))), "at least two")
  expect_error(match_zoi_across_files(list(df(1, 1), 42)), "data frame with")
  expect_error(match_zoi_across_files(list(df(1, 1), data.frame(mz = 1))), "data frame with")
  expect_error(match_zoi_across_files(list(data.frame(mz = "a", rt = 1), df(1, 1))), "numeric")
  expect_error(
    match_zoi_across_files(list(data.frame(mz = c(1, 2), rt = c(1, NA_real_)), df(1, 1))),
    "NA"
  )
})

test_that("scalar arguments are validated", {
  ok <- list(df(1, 1), df(1, 1))
  expect_error(match_zoi_across_files(ok, mz_tol = 0), "positive")
  expect_error(match_zoi_across_files(ok, mz_tol = -1), "positive")
  expect_error(match_zoi_across_files(ok, rt_tol = c(1, 2)), "single")
  expect_error(match_zoi_across_files(ok, legacy_first_match = NA), "TRUE or FALSE")
  expect_error(match_zoi_across_files(ok, legacy_first_match = 1), "TRUE or FALSE")
})
