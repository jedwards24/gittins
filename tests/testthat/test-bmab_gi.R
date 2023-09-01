test_that("bmab_gi gives expected output", {
  x1 <- bmab_gi_ab(3, 2, 0.8, 1e-4, 20)
  x2 <- bmab_gi(3, 5, 0.8, 1e-4, 20)
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
