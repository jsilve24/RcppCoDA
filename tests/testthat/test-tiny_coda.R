context("Tiny CoDA Functions")

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

test_that("alrContrast correct", {
  B.inv.true <- cbind(diag(3), 0)
  B.true <- cbind(diag(3), -1)
  expect_equal(B.inv.true, alrContrast(4,4, TRUE))
  expect_equal(B.true, alrContrast(4,4,FALSE))
  expect_error(alrContrast(4, 5, inv=TRUE))
  
  B.true <- matrix(0, 4, 5)
  foo <- cbind(1:4, c(1,2,4,5))
  B.true[foo] <- 1
  B.true[,3] <- -1
  expect_equal(B.true, alrContrast(5, 3, FALSE))
}) 

test_that("ilrContrast correct", {
  B <- ilrContrast(4)
  B.test <- t(qr.Q(qr(t(alrContrast(4,4,FALSE)))))
  expect_equal(B, B.test)
  expect_equal(B%*%t(B), diag(3))
  expect_equal(t(B)%*%B, diag(4)-1/4)
})


test_that("clrContrast correct", {
  B <- matrix(-1, 4, 4) + 4*diag(4)
  B <- B/4
  expect_equal(clrContrast(4, FALSE), B)
  expect_equal(clrContrast(4, TRUE), diag(4))
})


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
  Y <- glr(X, V, 2)
  expect_equal(Y[,,2], t(glr(t(X[,,2]), V)))
  expect_equal(X, glrInv(Y, V, 2))
})

test_that("alr and alrInv are inverses and correct", {
  X <- abs(matrix(rnorm(10), 5, 2))
  X <- clo(X)
  Y <- alr(X, 3)
  Y.manual <- X[-3,]
  Y.manual <- sweep(Y.manual, 2, X[3,], FUN="/")
  Y.manual <- log(Y.manual)
  # Correctness
  expect_equal(Y, Y.manual)

  # Inverse - and alrInv by association
  expect_equal(X, alrInv(Y, 3))
  
  X <- array(X, c(3,3,3))
  X <- clo(X, b=2)
  expect_error(alr(X, 10))
  Y <- alr(X, 2, 2)
  expect_equal(Y[,,2], t(alr(t(X[,,2]), 2)))
  expect_equal(X, alrInv(Y, 2, 2))
})


test_that("clr and clrInv are inverses and correct", {
  X <- abs(matrix(rnorm(10), 5, 2))
  X <- clo(X)
  Y <- clr(X)
  Y.manual <- glr(X, clrContrast(5, FALSE))
  # Correctness
  expect_equal(Y, Y.manual)
  
  # Inverse - and alrInv by association
  expect_equal(X, clrInv(Y))
  
  X <- array(X, c(3,3,3))
  X <- clo(X, b=2)
  Y <- clr(X, 2)
  expect_equal(Y[,,2], t(clr(t(X[,,2]))))
  expect_equal(X, clrInv(Y, 2))
})

test_that("ilr and ilrInv are inverses and correct", {
  X <- abs(matrix(rnorm(10), 5, 2))
  X <- clo(X)
  Y <- ilr(X)
  Y.manual <- ilr(X, ilrContrast(5))
  # Correctness
  expect_equal(Y, Y.manual)
  
  # Inverse - and alrInv by association
  expect_equal(X, ilrInv(Y))
  
  X <- array(X, c(3,3,3))
  X <- clo(X, b=2)
  expect_error(ilr(X, matrix(1:10, 5,2)))
  Y <- ilr(X, b=2)
  expect_equal(Y[,,2], t(ilr(t(X[,,2]))))
  expect_equal(X, ilrInv(Y, b=2))
})


test_that("array_pre and array_post correctness", {
  x  <- array(1:24, 2:4)
  y <- array_pre(x, 2)
  y.manual <- rbind(
    c(1,2,7,8,13,14,19,20), 
    c(3,4,9,10,15,16,21,22),
    c(5,6,11,12,17,18,23,24)
  )
  expect_equal(y, y.manual)
  expect_equal(array_post(y, 2, 2:4), x)
})

