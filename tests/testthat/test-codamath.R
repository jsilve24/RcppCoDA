x <- matrix(runif(300), 100, 3)
x <- clo(x)

test_that("clo sums to 1", {
  expect_equal(colSums(x), rep(1, ncol(x)))
})


test_that("clo correct on vectors", {
  x1 <- x[,1,drop=T]
  y <- clo(x1)
  expect_equal(matrix(x1/sum(x1)), y)
  expect_equal(sum(y), 1)
})


test_that("clo correct on matricies", {
  y <- clo(x)
  y.true <- sweep(x, 2, colSums(x), FUN=`/`)
  expect_equal(y.true, y)
  expect_equal(colSums(y), rep(1, 3))
})


test_that("center correct on vectors", {
  x <- runif(100)
  expect_equal(matrix(x-mean(x)), center(x))
  expect_equal(sum(center(x)), 0)
})


test_that("center correct on matricies", {
  x <- matrix(runif(300), 100, 3)
  y <- center(x)
  y.true <- sweep(x, 2, mean(x), FUN=`-`)
  expect_equal(center(x), y.true)
  expect_equal(colSums(y), rep(0, 3))
})
