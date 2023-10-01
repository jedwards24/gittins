test_that("bmab_gi gives expected output", {
  x1 <- bmab_gi_ab(3, 2, gamma = 0.8, N = 20, tol = 1e-4)
  x2 <- bmab_gi(3, 5, gamma = 0.8, N = 20, tol = 1e-4)
  expect_identical(x1, x2)
  expect_equal(x1, 0.67143994)
})

test_that("bmab_kgi gives expected output", {
  expect_equal(bmab_kgi(2, 4, 0.9), 0.58181818)
})

test_that("bmab_giplus gives expected output", {
  x3 <- bmab_giplus(3, 5, 0.8, 1e-4)
  x4 <- bmab_giplus(3, 5, 0.8, 1e-4, TRUE)
  expect_gte(x4, x3)
  expect_equal(x3, 0.7168457)
})

test_that("bmab_gi() arg checks work", {
  expect_error(bmab_gi(3, 5, gamma = 0.8, N = 1),"`N` must be at least 2")
  expect_error(bmab_gi(3, 5, gamma = 0.8, N = 9.5),"`N` must be a whole number")
  expect_error(bmab_gi(3, 5, gamma = 0, N = 1),"`gamma` must be greater than 0")
  expect_error(bmab_gi(Sigma = 1, n = 1, gamma = 0.5, N = 20), "`n` must be greater than `Sigma`")
})
