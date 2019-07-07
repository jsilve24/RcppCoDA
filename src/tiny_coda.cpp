#include "tiny_coda.h"

using namespace Rcpp;

using namespace Rcpp;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::Map;

//' Closure divide each column by column sum
//' @param X matrix (D x N)
//' @return matrix of same dimension as X
//' @export
// [[Rcpp::export]]
Eigen::MatrixXd clo(Eigen::Map<Eigen::MatrixXd> X){
  return coda::clo(X);
}

//' Center - subtract each column by column mean
//' @param X matrix (D x N)
//' @return matrix of same dimension as X
//' @export
// [[Rcpp::export]]
Eigen::MatrixXd center(Eigen::Map<Eigen::MatrixXd> X){
  return coda::center(X);
}

//' Log-Ratio Transformation
//'
//' \code{glr} is generic log-ratio transform, code used by other
//' transforms, can be called directly. *Contrast functions produce contrast 
//' matricies (e.g., "V") that define the various transforms. 
//'
//' @param X vector or matrix (rows are parts/coordinates, columns are samples) of data
//' @param V transformation matrix (defines transform; P x D where D is number of parts)
//'   if NULL then uses ilr default basis (see details)
//' @param d for ALR, which component (integer position) to take as reference
//' (default is nrow(x)) for alrInv corresponds to column position in untransformed
//' matrix.
//' @param inv for ALR and CLR, transformation matrix is different forward and inverse
//' @param D the number of parts (e.g., number of columns in untransformed data)
//' @return matrix (converts vectors to column matricies)
//' @details The default ILR base formed by Gram-Schmidt orthogonalization of an ALR_D basis.
//' @name base_lr_transforms
//' @examples
//' #ALR Transform
//' x <- matrix(runif(30), 10, 3)
//' x <- clo(x)
//' x.alr <- alr(x, 2)
//' x <- alrInv(x.alr, 2)
//'
//' # ILR
//' x.ilr <- ilr(x)
//' x <- ilrInv(x.ilr)
//'
//' # CLR
//' x.clr <- clr(x)
//' x <- clrInv(x.clr)
//'
//' # CUSTOM - Be careful if your custom matrix is not
//' # orthogonal the inverse transform may not be given by just the transpose!
//' # For example, this is the case for the ALR
//' V <- matrix(c(1, 1, -1), 1, 3)
//' x.custom <- glr(x, V)
//' @export
// [[Rcpp::export]]
Eigen::MatrixXd alrContrast(int D, int d, bool inv){
  return coda::alrContrast(D, d, inv);
}

//' @rdname base_lr_transforms
//' @export
// [[Rcpp::export]]
Eigen::MatrixXd clrContrast(int D, bool inv){
  return coda::clrContrast(D, inv);
}

//' @rdname base_lr_transforms
//' @export
// [[Rcpp::export]]
Eigen::MatrixXd ilrContrast(int D){
  return coda::ilrContrast(D);
}

//' @rdname base_lr_transforms
//' @export
// [[Rcpp::export]]
Eigen::MatrixXd glr(Eigen::Map<Eigen::MatrixXd>& X, Eigen::Map<Eigen::MatrixXd>& V){
  return coda::glr(X, V);
}

//' @rdname base_lr_transforms
//' @export
// [[Rcpp::export]]
Eigen::MatrixXd glrInv(Eigen::Map<Eigen::MatrixXd>& X, Eigen::Map<Eigen::MatrixXd>& V){
  return coda::glrInv(X,V);
}

//' @rdname base_lr_transforms
//' @export
// [[Rcpp::export]]
Eigen::MatrixXd alr(Eigen::Map<Eigen::MatrixXd>& X, int d){
  return coda::alr(X,d);
}

//' @rdname base_lr_transforms
//' @export
// [[Rcpp::export]]
Eigen::MatrixXd alrInv(Eigen::Map<Eigen::MatrixXd>& X, int d){
  return coda::alrInv(X,d);
}

//' @rdname base_lr_transforms
//' @export
// [[Rcpp::export]]
Eigen::MatrixXd clr(Eigen::Map<Eigen::MatrixXd>& X){
  return coda::clr(X);
}

//' @rdname base_lr_transforms
//' @export
// [[Rcpp::export]]
Eigen::MatrixXd clrInv(Eigen::Map<Eigen::MatrixXd>& X){
  return coda::clrInv(X);
}

//' @rdname base_lr_transforms
//' @export
// [[Rcpp::export]]
Eigen::MatrixXd ilr(Eigen::Map<Eigen::MatrixXd>& X,
                    Rcpp::Nullable<Eigen::Map<Eigen::MatrixXd>&> V = R_NilValue){
  if (V.isNull()){
    return coda::ilr(X);
  } else {
    MatrixXd VV = Rcpp::as<Eigen::MatrixXd>(V);
    return coda::ilr(X, VV);
  }
}

//' @rdname base_lr_transforms
//' @export
// [[Rcpp::export]]
Eigen::MatrixXd ilrInv(Eigen::Map<Eigen::MatrixXd>& X,
                       Rcpp::Nullable<Eigen::Map<Eigen::MatrixXd>&> V = R_NilValue){
  if (V.isNull()){
    return coda::ilrInv(X);
  } else {
    MatrixXd VV = Rcpp::as<Eigen::MatrixXd>(V);
    return coda::ilrInv(X, VV);
  }
}
