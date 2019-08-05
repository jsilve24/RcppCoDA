
#' Generalized LR Transform
#' @param X vector, matrix, or array of data (b is the dimension which is compositionally relevant)
#' @param V1 (P1 X D1, e.g., transfer Contrast matrix)
#' @param V2 (P2 x D2, e.g., transfer Contrast matrix)
#' @param b index of dimension to operate on (e.g., index of dimension of parts or coords in X;
#' default is 1 meaning that compositions/log-ratios are rows)
#' @return Eigen::MatrixXd
#' @details calculates
#' 
#'   let Y = cbind[Y1, Y2] 
#'   let X = cbind[X1, X2] 
#'   Calculates Y = GLR(X) or Y=GLRINV(X) defined by 
#'   Y1 = GLR(X1, V1) or Y1 = GLRINV(X1, V1)
#'   Y2 = GLR(X2, V2) or Y2 = GLRINV(X2, V2)
#'   
#' If V2 is not given then assumes that Y2 = X2. 
#' In that case determines the size of X2 and Y2 based on on size of V1. 
#' @examples 
#' # CUSTOM - Be careful if your custom matrix is not
#' # orthogonal the inverse transform may not be given by just the transpose!
#' # For example, this is the case for the ALR
#' x <- matrix(runif(30), 10, 3)
#' x <- clo(x)
#' V <- matrix(c(1, 1, -1), 1, 3)
#' x.custom <- glr(x, V)
#' @name glr_transforms
NULL

#' @rdname glr_transforms
#' @export
glr <- function(X, V1, V2=NULL, b=1){
  X <- vec_to_array(X)
  s <- dim(X)
  if (is.null(V2)){
    p <- nrow(V1) + s[b] - ncol(V1)  
  } else {
    p <- nrow(V1) + nrow(V2)
  }
  X <- array_pre(X, b)
  X <- glr_internal(X, V1, V2)
  s[b] <- p
  X <- array_post(X, b, s)
  return(X)
}

#' @rdname glr_transforms
#' @export
glrInv <- function(X, V1, V2=NULL, b=1){
  X <- vec_to_array(X)
  s <- dim(X)
  if (is.null(V2)){
    p <- ncol(V1) + s[b] - nrow(V1)  
  } else {
    p <- ncol(V1) + ncol(V2)
  }
  X <- array_pre(X, b)
  X <- glrInv_internal(X, V1, V2)
  s[b] <- p
  X <- array_post(X, b, s)
  return(X)
}