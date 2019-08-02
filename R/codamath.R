#' Linear form Y=VX
#' @param X vector, matrix, or array of data (b is the dimension on which the linear form acts 
#'    e.g., the dimension on which V1 and V2 are multiplied). 
#' @param V1 (P1 X D1, e.g., transfer Contrast matrix)
#' @param V2 (P2 x D2, e.g., transfer Contrast matrix)
#' @param b (index of dimension to operate on; default: 1)
#' @details Calculates
#'
#'   Y = bdiag[V1, V2]  X
#' where bdiag a block diagonal matrix with blocks V1 and V2.
#' If V2 is not given then assumes that second block is identity.
#' In that case determines the size of second block based on size of 1.
#'
#' @return Y
#' @export
#' @name linForm
linForm <- function(X, V1, V2=NULL, b=1){
  X <- vec_to_array(X)
  s <- dim(X)
  if (is.null(V2)){
    p <- nrow(V1) + s[b] - ncol(V1)  
  } else {
    p <- nrow(V1) + nrow(V2)
  }
  X <- array_pre(X, b)
  X <- linForm_internal(X, V1, V2)
  s[b] <- p
  X <- array_post(X, b, s)
  return(X)
}


#' Quadratic Form Y = V X V^T
#' @param X matrix or array (If array, assumes first two dimensions 
#'   are the ones to operate on, e.g., X_i = X[,,i] if 3D array)
#' @param V1 (e.g., transfer Contrast matrix)
#' @param V2 (e.g., transfer Contrast matrix)
#' @return Y
#' @details Calculates:
#'
#' Y_i = bdiag[V1, V2]  X_i  bdiag[V1, V2]^T
#'
#' where Y_i and X_i refer to the i-th set of D-columns of Y and X respectively
#' and where bdiag a block diagonal matrix with blocks V1 and V2.
#' If V2 is not given then assumes that second block is identity.
#' In that case determines the size of second block based on size of V1.
#'
#' @name quadForm
#' @export
quadForm <- function(X, V1, V2=NULL){
  b <- 1:2
  X <- vec_to_array(X)
  s <- dim(X)
  if (is.null(V2)){
    p <- nrow(V1) + s[1] - ncol(V1)  
  } else {
    p <- nrow(V1) + nrow(V2)
  }
  X <- array_pre(X, b)
  X <- quadForm_internal(X, V1, V2)
  s[b] <- p
  X <- array_post(X, b, s)
  return(X)
}