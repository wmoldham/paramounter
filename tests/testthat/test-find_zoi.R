# test-find_zoi.R

test_that("no points above the cutoff returns a zero-row frame with integer columns", {
  z <- find_zoi(c(1, 2, 3), cutoff = 100)
  expect_equal(nrow(z), 0L)
  expect_type(z$start, "integer")
  expect_type(z$end, "integer")
  expect_type(z$apex, "integer")
})

test_that("an empty intensity vector returns zero rows", {
  expect_equal(nrow(find_zoi(numeric(0), cutoff = 5)), 0L)
})

test_that("contiguous runs and their apexes are identified", {
  z <- find_zoi(c(2, 1, 50, 80, 40, 1, 1, 90, 30, 1, 100), cutoff = 10)
  expect_identical(z$start, c(3L, 8L, 11L))
  expect_identical(z$end, c(5L, 9L, 11L))
  expect_identical(z$apex, c(4L, 8L, 11L))
})

test_that("a trace entirely above the cutoff is a single ZOI", {
  z <- find_zoi(c(10, 20, 99, 20, 10), cutoff = 5)
  expect_identical(z$start, 1L)
  expect_identical(z$end, 5L)
  expect_identical(z$apex, 3L)
})

test_that("apex ties resolve to the first maximum", {
  expect_identical(find_zoi(c(1, 77, 77, 1), cutoff = 10)$apex, 2L)
})

test_that("the cutoff is exclusive", {
  # points exactly equal to the cutoff are not part of any ZOI
  expect_equal(nrow(find_zoi(c(5, 5, 5), cutoff = 5)), 0L)
})

test_that("a zero cutoff keeps positive points and drops zeros", {
  z <- find_zoi(c(0, 3, 4, 0, 5), cutoff = 0)
  expect_identical(z$start, c(2L, 5L))
  expect_identical(z$end, c(3L, 5L))
  expect_identical(z$apex, c(3L, 5L))
})

test_that("result columns are integer", {
  z <- find_zoi(c(0, 3, 4, 0, 5), cutoff = 0)
  expect_type(z$start, "integer")
  expect_type(z$end, "integer")
  expect_type(z$apex, "integer")
})

test_that("matches the original segmentation across random inputs", {
  orig_zoi <- function(eic, cut) {
    a <- which(eic > cut)
    if (length(a) == 0) return(list(s = integer(0), e = integer(0), p = integer(0)))
    g <- split(a, cumsum(c(1, diff(a) != 1)))
    p <- integer(length(g))
    s <- integer(length(g))
    e <- integer(length(g))
    for (x in seq_along(g)) {
      seg <- g[[x]]
      p[x] <- which(eic[seg] == max(eic[seg]))[1] + min(seg) - 1
      s[x] <- min(seg)
      e[x] <- max(seg)
    }
    list(s = s, e = e, p = p)
  }
  set.seed(11)
  for (rep in 1:500) {
    L <- sample(3:400, 1)
    x <- round(runif(L, 0, 1000), 2)
    sp <- sample(L, min(L, sample(0:8, 1)))
    x[sp] <- x[sp] + runif(length(sp), 1000, 5000)
    cut <- as.numeric(quantile(x, runif(1, 0.3, 0.95)))
    a <- orig_zoi(x, cut)
    b <- find_zoi(x, cut)
    expect_equal(b$start, a$s)
    expect_equal(b$end, a$e)
    expect_equal(b$apex, a$p)
  }
})

test_that("non-numeric intensity errors", {
  expect_error(find_zoi("a", 1), "numeric")
  expect_error(find_zoi(list(1, 2), 1), "numeric")
  expect_error(find_zoi(factor(c(1, 2)), 1), "numeric")
})

test_that("non-finite intensity errors", {
  expect_error(find_zoi(c(1, NA, 3), 1), "NA")
  expect_error(find_zoi(c(1, NaN, 3), 1))
  expect_error(find_zoi(c(1, Inf, 3), 1))
})

test_that("invalid cutoff errors", {
  x <- c(1, 2, 3)
  expect_error(find_zoi(x, -1), "non-negative")
  expect_error(find_zoi(x, c(1, 2)), "single")
  expect_error(find_zoi(x, NA_real_))
  expect_error(find_zoi(x, "10"))
  expect_error(find_zoi(x, Inf))
})
