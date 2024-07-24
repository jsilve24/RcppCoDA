test_that("glr and glrInv are inverses and correct", {
  X <- abs(matrix(rnorm(10), 5, 2))
  X <- clo(X)
  V <- ilrContrast(5)
  Y <- glr(X, V)
  
  # They are inverses
  expect_equal(X, glrInv(Y, V))
  
  # glr Correct and therepy glrInv correct given they are inverses
  expect_equal(Y, V%*%log(X))
})
