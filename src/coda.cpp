#include "CoDA.h"

using namespace Rcpp;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::Map;

// ----- CODA MATH


//' Closure operation - divide elements by sum of elements in b dimension
//' @param X vector, matrix, or array of data 
//'   (b is the dimension which is relevant to be closed)
//' @param b index of dimension to operate on 
//'   (e.g., index of dimension of parts or coords in X;
//'   default is 1 meaning that compositions/log-ratios are rows)
//' @export
// [[Rcpp::export]]
Eigen::MatrixXd clo(Eigen::Map<Eigen::MatrixXd> X){
  return coda::clo(X);
}


//' Center operation - subtract from each element the mean of elements in b dimension
//' @param X vector, matrix, or array of data 
//'   (b is the dimension which is relevant to be centered)
//' @param b index of dimension to operate on 
//'   (e.g., index of dimension of parts or coords in X;
//'   default is 1 meaning that compositions/log-ratios are rows)
//' @export
// [[Rcpp::export]]
Eigen::MatrixXd center(Eigen::Map<Eigen::MatrixXd> X){
  return coda::center(X);
}

// ------- GLR TRANSFORMS



//' Log-Ratio Transformation
//'
//' \code{glr} is generic log-ratio transform, code used by other
//' transforms, can be called directly. *Contrast functions produce contrast 
//' matricies (e.g., "V") that define the various transforms. 
//'
//' @param X vector, matrix, or array of data (b is the dimension which is compositionally relevant)
//' @param V transformation matrix (defines transform; P x D where D is number of parts)
//'   if NULL then uses ilr default basis (see details)
//' @param d for ALR, which component (integer position) to take as reference
//' (default is nrow(x) for alr or nrow(x)+1 for alrInv) for alrInv corresponds 
//' to column position in untransformed matrix. ilrContrast can accept a D-1 x D sign matrix to produce a 
//' Sequential Binary Partition basis. 
//' @param b index of dimension to operate on (e.g., index of dimension of parts or coords in X;
//' default is 1 meaning that compositions/log-ratios are rows)
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
//' @export
// [[Rcpp::export]]
Eigen::MatrixXd glr(Eigen::Map<Eigen::MatrixXd>& X, Eigen::Map<Eigen::MatrixXd>& V){
  return coda::glr(X, V);
}

//' @rdname base_lr_transforms
//' @export
// [[Rcpp::export]]
Eigen::MatrixXd glrInv(Eigen::Map<Eigen::MatrixXd>& X, Eigen::Map<Eigen::MatrixXd>& V){
  return coda::glrInv(X, V);
}



// ------- CODA TRANSFORMS

//' @rdname base_lr_transforms
//' @export
// [[Rcpp::export]]
Eigen::MatrixXd alrContrast(int d, int D, bool inv){
  return coda::alrContrast(d, D, inv);
}

//' @rdname base_lr_transforms
//' @export
// [[Rcpp::export]]
Eigen::MatrixXd clrContrast(int D, bool inv){
  return coda::clrContrast(D, inv);
}

// [[Rcpp::export]]
Eigen::MatrixXd ilrContrastDefault_internal(int D){
  return coda::ilrContrast(D);
}


// [[Rcpp::export]]
Eigen::MatrixXd ilrContrastSign_internal(Eigen::Map<Eigen::MatrixXi>& S){
  return coda::ilrContrast(S);
}


//' @rdname base_lr_transforms
//' @export
// [[Rcpp::export]]
Eigen::MatrixXd alr(Eigen::Map<Eigen::MatrixXd>& X, Rcpp::Nullable<int> d = R_NilValue){
  int d_;
  if (d.isNull()){
    d_ = X.rows();
  } else {
    d_ = Rcpp::as<int>(d);
  }
  return coda::alr(X,d_);
}


//' @rdname base_lr_transforms
//' @export
// [[Rcpp::export]]
Eigen::MatrixXd alrInv(Eigen::Map<Eigen::MatrixXd>& X, Rcpp::Nullable<int> d = R_NilValue){
  int d_;
  if (d.isNull()){
    d_ = X.rows()+1;
  } else {
    d_ = Rcpp::as<int>(d);
  }
  return coda::alrInv(X,d_);
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
