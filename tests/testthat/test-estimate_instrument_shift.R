# test-estimate_instrument_shift.R

test_that("two-file shifts are the absolute per-feature differences", {
  r <- estimate_instrument_shift(list(
    mz = rbind(c(200.000, 200.003), c(300.000, 300.002)),
    rt = rbind(c(500, 505), c(800, 803))
  ))
  expect_equal(r$mass_shift, c(0.003, 0.002))
  expect_equal(r$rt_shift, c(5, 3))
})

test_that("shifts span the full range across three files", {
  r <- estimate_instrument_shift(list(
    mz = matrix(c(200, 200.003, 200.005), 1),
    rt = matrix(c(500, 505, 510), 1)
  ))
  expect_equal(r$mass_shift, 0.005)
  expect_equal(r$rt_shift, 10)
})

test_that("no matched features gives empty distributions", {
  r <- estimate_instrument_shift(list(
    mz = matrix(numeric(0), 0, 3),
    rt = matrix(numeric(0), 0, 3)
  ))
  expect_identical(r$mass_shift, numeric(0))
  expect_identical(r$rt_shift, numeric(0))
})

test_that("identical values across files give zero shift", {
  r <- estimate_instrument_shift(list(
    mz = matrix(c(200, 200, 200), 1),
    rt = matrix(c(9, 9, 9), 1)
  ))
  expect_equal(r$mass_shift, 0)
  expect_equal(r$rt_shift, 0)
})

test_that("shifts match a naive per-row max-minus-min across random matrices", {
  set.seed(7)
  for (rep in 1:500) {
    nr <- sample(0:30, 1)
    nc <- sample(2:6, 1)
    mz <- matrix(runif(nr * nc, 100, 900), nr, nc)
    rt <- matrix(runif(nr * nc, 0, 1500), nr, nc)
    r <- estimate_instrument_shift(list(mz = mz, rt = rt))
    expected_mass <- if (nr == 0) {
      numeric(0)
    } else {
      vapply(seq_len(nr), function(i) max(mz[i, ]) - min(mz[i, ]), numeric(1))
    }
    expected_rt <- if (nr == 0) {
      numeric(0)
    } else {
      vapply(seq_len(nr), function(i) max(rt[i, ]) - min(rt[i, ]), numeric(1))
    }
    expect_equal(r$mass_shift, expected_mass)
    expect_equal(r$rt_shift, expected_rt)
  }
})

test_that("the matched structure is validated", {
  expect_error(estimate_instrument_shift(42), "list with")
  expect_error(estimate_instrument_shift(list(mz = matrix(1, 1, 2))), "list with")
  expect_error(estimate_instrument_shift(list(mz = c(1, 2), rt = matrix(1, 1, 2))), "numeric matrix")
  expect_error(estimate_instrument_shift(list(mz = matrix(1, 1, 2), rt = c(1, 2))), "numeric matrix")
  expect_error(estimate_instrument_shift(list(mz = matrix("a", 1, 2), rt = matrix(1, 1, 2))), "numeric matrix")
  expect_error(estimate_instrument_shift(list(mz = matrix(1, 2, 2), rt = matrix(1, 1, 2))), "same dimensions")
  expect_error(estimate_instrument_shift(list(mz = matrix(1, 1, 1), rt = matrix(1, 1, 1))), "at least two")
  expect_error(estimate_instrument_shift(list(mz = matrix(c(1, NA), 1, 2), rt = matrix(1, 1, 2))), "NA")
  expect_error(estimate_instrument_shift(list(mz = matrix(1, 1, 2), rt = matrix(c(1, Inf), 1, 2))), "NA")
})
