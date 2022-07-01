test_that("check_prev_out passes if all check_names are present", {
  prev_out <- list(a = 1, b = 2, c = 3)
  check_prev_out(prev_out, c("a", "b", "c"))
  # avoid skipping empty test
  expect_true(TRUE)
})

test_that("check_prev_out fails if not all check_names are present", {
  prev_out <- list(a = 1, b = 2, c = 3)
  expect_error(check_prev_out(prev_out, c("a", "b", "c", "d")))
})
