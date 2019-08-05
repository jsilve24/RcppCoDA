# Simple Cases
test_that("glr and glrInv are inverses and correct", {
  X <- abs(matrix(rnorm(10), 5, 2))
  X <- clo(X)
  V <- ilrContrast(5)
  Y <- glr(X, V)
  
  # They are inverses
  expect_equal(X, glrInv(Y, V))
  
  # glr Correct and therepy glrInv correct given they are inverses
  expect_equal(Y, V%*%log(X))
  
  X <- array(X, c(3,3,3))
  X <- clo(X, b=2)
  V <- ilrContrast(5)
  expect_error(glr(X, V))
  V <- ilrContrast(3)
  Y <- glr(X, V, b=2)
  expect_equal(Y[,,2], t(glr(t(X[,,2]), V)))
  expect_equal(X, glrInv(Y, V, b=2))
})

# More Difficult Cases
test_that("glr Works", {
  X <- abs(matrix(rnorm(10), 5, 2))
  X[1:3,] <- clo(X[1:3,])
  X[4:5,] <- clo(X[4:5,])
  V1 <- ilrContrast(3)
  V2 <- ilrContrast(2)
  
  Y <- matrix(0, 4, 2)
  Y[1:2,] <- glr(X[1:3,], V1)
  Y[3:4,] <- X[4:5,]
  expect_equal(Y, glr(X, V1))
  expect_equal(X, glrInv(Y, V1))
  
  Y <- matrix(0, 3, 2)
  Y[1:2,] <- glr(X[1:3,], V1)
  Y[3,] <- glr(X[4:5,], V2)
  expect_equal(Y, glr(X, V1, V2))
  expect_equal(X, glrInv(Y, V1, V2))
})