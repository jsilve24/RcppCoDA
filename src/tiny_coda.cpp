#include "tiny_coda.h"

using namespace Rcpp;

using namespace Rcpp;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::Map;

// [[Rcpp::export]]
Eigen::MatrixXd clo_internal(Eigen::Map<Eigen::MatrixXd> X){
  return coda::clo(X);
}

// [[Rcpp::export]]
Eigen::MatrixXd center_internal(Eigen::Map<Eigen::MatrixXd> X){
  return coda::center(X);
}

//' @rdname base_lr_transforms
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

// [[Rcpp::export]]
Eigen::MatrixXd glr_internal(Eigen::Map<Eigen::MatrixXd>& X, Eigen::Map<Eigen::MatrixXd>& V){
  return coda::glr(X, V);
}

// [[Rcpp::export]]
Eigen::MatrixXd glrInv_internal(Eigen::Map<Eigen::MatrixXd>& X, Eigen::Map<Eigen::MatrixXd>& V){
  return coda::glrInv(X,V);
}

// [[Rcpp::export]]
Eigen::MatrixXd alr_internal(Eigen::Map<Eigen::MatrixXd>& X, int d){
  return coda::alr(X,d);
}

// [[Rcpp::export]]
Eigen::MatrixXd alrInv_internal(Eigen::Map<Eigen::MatrixXd>& X, int d){
  return coda::alrInv(X,d);
}

// [[Rcpp::export]]
Eigen::MatrixXd clr_internal(Eigen::Map<Eigen::MatrixXd>& X){
  return coda::clr(X);
}

// [[Rcpp::export]]
Eigen::MatrixXd clrInv_internal(Eigen::Map<Eigen::MatrixXd>& X){
  return coda::clrInv(X);
}

// [[Rcpp::export]]
Eigen::MatrixXd ilr_internal(Eigen::Map<Eigen::MatrixXd>& X,
                    Rcpp::Nullable<Eigen::Map<Eigen::MatrixXd>&> V = R_NilValue){
  if (V.isNull()){
    return coda::ilr(X);
  } else {
    MatrixXd VV = Rcpp::as<Eigen::MatrixXd>(V);
    return coda::ilr(X, VV);
  }
}

// [[Rcpp::export]]
Eigen::MatrixXd ilrInv_internal(Eigen::Map<Eigen::MatrixXd>& X,
                       Rcpp::Nullable<Eigen::Map<Eigen::MatrixXd>&> V = R_NilValue){
  if (V.isNull()){
    return coda::ilrInv(X);
  } else {
    MatrixXd VV = Rcpp::as<Eigen::MatrixXd>(V);
    return coda::ilrInv(X, VV);
  }
}
