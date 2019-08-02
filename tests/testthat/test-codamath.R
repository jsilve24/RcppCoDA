x <- matrix(runif(300), 100, 3)
x <- clo(x)

test_that("clo sums to 1", {
  expect_equal(colSums(x), rep(1, ncol(x)))
})

test_that("clo correct on arrays", {
  x <- array(x, dim=c(2,4,3))
  y <- clo(x)
  expect_true(all(abs(colSums(y)-1)<1e-10))
  y <- clo(x, b=2)
  expect_false(all(abs(colSums(y)-1)<1e-10))
  expect_true(  
    all(abs(rowSums(y[,,1])-1) < 1e-10) &
      all(abs(rowSums(y[,,2])-1) < 1e-10) &
      all(abs(rowSums(y[,,3])-1) < 1e-10)
  )
})

test_that("center sums to 0", {
  x <- matrix(runif(300), 100, 3)
  expect_equal(colSums(center(x)), rep(0, ncol(x)))
  
  x <- array(x, dim=c(2,4,3))
  y <- center(x)
  expect_true(all(abs(colSums(y))<1e-10))
  y <- center(x, b=2)
  expect_false(all(abs(colSums(y))<1e-10))
  expect_true(  
    all(abs(rowSums(y[,,1])) < 1e-10) &
      all(abs(rowSums(y[,,2])) < 1e-10) &
      all(abs(rowSums(y[,,3])) < 1e-10)
  )
})


test_that("linForm Correct", {
  X <- matrix(runif(300), 100, 3)
  X <- array(x, dim=c(2,4,3))
  
  # Test Simple case
  V <- matrix(rnorm(4*5), 5,4)
  Y = array(0, dim=c(c(2,5,3)))
  for (i in 1:2){
    Y[i,,] <- V%*%X[i,,]
  }
  expect_error(linForm(X, V))
  expect_equal(Y, linForm(X, V,b=2))
  
  # Test case where V2=I
  V <- matrix(rnorm(5*3), 5,3)
  Y <- array(0, dim=c(2, 6, 3))
  for (i in 1:2){
    Y[i,1:5,] <- V%*%X[i,1:3,]
    Y[i,6,] <- X[i,4,]
  }
  expect_equal(Y, linForm(X, V, b=2))
  
  # test simple case where V1 and V2 are given
  V1 <- matrix(rnorm(5*3), 5,3)
  V2 <- matrix(rnorm(1))
  X <- matrix(rnorm(4*2), 4, 2)
  Y <- matrix(0, 6, 2)
  Y[1:5,] <- V1 %*% X[1:3,]
  Y[6,] <- V2 %*% X[4,]
  expect_equal(Y, linForm(X, V1, V2))
  
  # Test case where V2 is given
  X <- matrix(runif(300), 100, 3)
  X <- array(x, dim=c(2,4,3))
  V1 <- matrix(rnorm(5*3), 5,3)
  V2 <- matrix(rnorm(1))
  Y <- array(0, dim=c(2, 6, 3))
  for (i in 1:2){
    Y[i,1:5,] <- V1%*%X[i,1:3,]
    Y[i,6,] <- V2%*%X[i,4,]
  }
  expect_equal(Y, linForm(X, V1,V2, b=2))
})


test_that("quadForm Correct", {
  X <- runif(5*5*3)
  X <- array(X, dim=c(5,5,3))
  for (i in 1:3){
    X[,,i] <- tcrossprod(X[,,i])
  }
  
  # Test Simple case
  V <- matrix(rnorm(6*5), 6,5)
  Y = array(0, dim=c(c(6,6,3)))
  for (i in 1:3){
    Y[,,i] <- V%*%X[,,i]%*%t(V)
  }
  expect_equal(Y, quadForm(X, V))
   
  # Test case where V2=I
  V <- matrix(rnorm(6*4), 6,4)
  Y = array(0, dim=c(c(7,7,3)))
  for (i in 1:3){
    Y[1:6,1:6,i] <- V%*%X[1:4,1:4,i]%*%t(V)
    Y[1:6,7, i] <- V %*% X[1:4,5,i]
    Y[7,1:6, i] <- t(Y[1:6,7, i])
    Y[7,7,i] <- X[5,5,i]
  }
  expect_equal(Y, quadForm(X, V))
  
  # Test case where V1 and V2 are specified
  V1 <- matrix(rnorm(6*4), 6,4)
  V2 <- matrix(rnorm(1))
  Y = array(0, dim=c(c(7,7,3)))
  for (i in 1:3){
    Y[1:6,1:6,i] <- V1%*%X[1:4,1:4,i]%*%t(V1)
    Y[1:6,7, i] <- V1 %*% X[1:4,5,i] %*% t(V2)
    Y[7,1:6, i] <- t(Y[1:6,7, i])
    Y[7,7,i] <- V2 %*%X[5,5,i] %*% t(V2)
  }
  expect_equal(Y, quadForm(X, V1, V2))
})



