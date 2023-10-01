test_that("`check_numeric()` works.", {
  expect_null(check_numeric(5, "arg"))
  expect_null(check_numeric(5, "arg", 1))
  expect_null(check_numeric(5, "arg", 1, 10))
  expect_error(check_numeric(0, "arg", 1, 10), "`arg` must be greater than 1")
  expect_error(check_numeric(1, "arg", 1, 10), "`arg` must be greater than 1")
  expect_error(check_numeric(5, "arg", 1, 5), "`arg` must be less than 5")
  expect_null(check_numeric(5.5, "arg", 1, 10))
})

test_that("`check_integerish()` works.", {
  expect_null(check_integerish(5, "arg"))
  expect_error(check_integerish(5.5, "arg"), "`arg` must be a whole number")
  expect_null(check_integerish(5, "arg", 5, 9))
  expect_null(check_integerish(5, "arg", 0, 5))
  expect_error(check_integerish(0, "arg", 1, 10), "`arg` must be at least 1")
  expect_error(check_integerish(6, "arg", 1, 5), "`arg` must be no greater than 5")
})
