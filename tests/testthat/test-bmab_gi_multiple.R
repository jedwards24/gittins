test_that("bmab_gi_multiple() arg checks work", {
  expect_true(check_start_args(alpha_start = 1, beta_start = 1))
  expect_false(check_start_args(alpha_start = 1))
  expect_false(check_start_args(alpha_start = 1, beta_start = 1, n_start = 1))
  expect_false(check_start_args(alpha_start = 1, beta_start = 1, n_start = 1, Sigma_start = 1))
  expect_true(check_start_args(n_start = 1, Sigma_start = 1))
  expect_false(check_start_args(alpha_start = 1, Sigma_start = 1))

  expect_error(bmab_gi_multiple(alpha_start = 1, gamma = 0.5, N = 20, num_actions = 2), "Incorrect number")
  expect_error(bmab_gi_multiple(1, 1, gamma = 0.5, N = 20, num_actions = 1.9),
               "`num_actions` must be at least 2")

})

test_that("bmab_gi_multiple() gives correct output", {
  b1 <- bmab_gi_multiple(alpha_start = 2, beta_start = 1, num_actions = 2, gamma = 0.5, N = 20)
  b2 <- bmab_gi_multiple(Sigma_start = 2, n_start = 3, num_actions = 2, gamma = 0.5, N = 20)
  expect_identical(b1, b2)
  expect_s3_class(b1, "data.frame")
  expect_true(all(c("params", "gi_matrix", "gi_matrix_ns") %in% names(attributes(b1))))
  expect_snapshot_value(b1, style = "json2")
})
