test_that("vec_to_array works", {
  expect_equal(vec_to_array(c(1:2)), matrix(1:2))
  foo <- c(1:2)
  names(foo) <- c("a", "b")
  bar <- vec_to_array(foo)
  expect_equal(rownames(bar), names(foo))
})


test_that("array_pre/post gives errors as expected", {
  expect_error(array_pre(c(1,2,3), 0))
  expect_error(array_pre(c(1,2,3), 5))
  expect_error(array_pre(c(1,2,3), NULL))
  expect_error(array_pre(c(1,2,3), 1:5))
  
  expect_error(array_post(c(1,2,3), 0, c(1)))
  expect_error(array_post(c(1,2,3), 5, c(1)))
  expect_error(array_post(c(1,2,3), 2:3, c(1)))
  expect_error(array_post(c(1,2,3), 2:4, c(1)))
})
