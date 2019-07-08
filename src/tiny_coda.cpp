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
  
//' Transfer Contrasts for transfering from one coordinate system to another
//' @param Sigma covariance matrix in specified transformed space
//' @param V ILR contrast matrix (i.e., transformation matrix of ILR)
//' @param V1 ILR contrast matrix of basis Sigma is already in
//' @param V2 ILR contrast matrix of basis Sigma is desired in
//' @param d1 alr reference element Sigma is already expressed with respec to
//' @param d2 alr reference element Sigma is to be expressed with respect to
//' @return matrix
//' @name convert_coda
//' @export
// [[Rcpp::export]]
Eigen::MatrixXd iiContrast(Eigen::Map<Eigen::MatrixXd>& V1, 
                           Eigen::Map<Eigen::MatrixXd>& V2){
  return coda::iiContrast(V1, V2);
}

//' @rdname convert_coda
//' @export
// [[Rcpp::export]]
Eigen::MatrixXd icContrast(Eigen::Map<Eigen::MatrixXd>& V1){
  return coda::icContrast(V1);
}

//' @rdname convert_coda
//' @export
// [[Rcpp::export]]
Eigen::MatrixXd ciContrast(Eigen::Map<Eigen::MatrixXd>& V2){
  return coda::ciContrast(V2);
}

//' @rdname convert_coda
//' @export
// [[Rcpp::export]]
Eigen::MatrixXd iaContrast(Eigen::Map<Eigen::MatrixXd>& V1, int d2, int D){
  return coda::iaContrast(V1, d2, D);
}

//' @rdname convert_coda
//' @export
// [[Rcpp::export]]
Eigen::MatrixXd aiContrast(int d2, Eigen::Map<Eigen::MatrixXd>& V2, int D){
  return coda::aiContrast(d2, V2, D);
}

//' @rdname convert_coda
//' @export
// [[Rcpp::export]]
Eigen::MatrixXd caContrast(int d2, int D){
  return coda::caContrast(d2, D);
}

//' @rdname convert_coda
//' @export
// [[Rcpp::export]]
Eigen::MatrixXd acContrast(int d1, int D){
  return coda::acContrast(d1, D);
}

//' @rdname convert_coda
//' @export
// [[Rcpp::export]]
Eigen::MatrixXd aaContrast(int d1, int d2, int D){
  return coda::aaContrast(d1, d2, D);
}



 
