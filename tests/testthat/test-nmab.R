test_that("nmab_gi gives expected output", {
  x1 <- nmab_gi(0, 1, 0.8, 1, 30, 3, 0.1, 1e-4)
  x2 <- nmab_gi(1, 1, 0.8, 1, 30, 3, 0.1, 1e-4)
  expect_identical(x1 + 1, x2)
  expect_equal(x1, 0.50552724)
})

test_that("nmab_kgi and nmab_giplus give expected output", {
  x3 <- nmab_kgi(0, 1, 0.8, 1, 1e-4)
  x4 <- nmab_giplus(0, 1, 0.8, 1e-4)
  x5 <- nmab_giplus(0, 1, 0.8, 1e-4, upper = TRUE)
  expect_gte(x5, x4)
  expect_equal(x3, 0.44973755)
  expect_equal(x4, 0.63601685)
})

