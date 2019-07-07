context("Tiny CoDA Functions")

x <- matrix(runif(300), 100, 3)
x <- clo(x)

test_that("clo sums to 1", {
  expect_equal(colSums(x), rep(1, ncol(x)))
})

test_that("center sums to 0", {
  x <- matrix(runif(300), 100, 3)
  expect_equal(colSums(center(x)), rep(0, ncol(x)))
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
})


