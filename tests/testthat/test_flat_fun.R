
test_that("get_mu_beta returns values < 1", {
  expect_true(get_mu_beta(0.3) < 1)
})

test_that("get_omega returns error when beta =  1", {
  expect_error(get_omega(1))
})
