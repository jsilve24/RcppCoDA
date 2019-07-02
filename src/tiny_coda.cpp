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
Eigen::MatrixXd miniclo(Eigen::Map<Eigen::MatrixXd> X){
  return coda::miniclo(X);
}

//' Create ALR contrast matrix
//' @param D total number of parts
//' @param d reference part for ALR
//' @param inv (is this being used for the alr or alrInv?)
//' @return D-1 x D matrix
//' @export
// [[Rcpp::export]]
Eigen::MatrixXd alrContrast(int D, int d, bool inv){
  return coda::alrContrast(D, d, inv);
}

//' Create DEFAULT ILR contrast matrix
//' @param D total number of parts
//' @return D-1 x D matrix (a basis)
//' @export
// [[Rcpp::export]]
Eigen::MatrixXd ilrContrast(int D){
  return coda::ilrContrast(D);
}

//' Create CLR contrast matrix
//' @param D total number of parts
//' @return D-1 x D matrix (a basis)
//' @export
// [[Rcpp::export]]
Eigen::MatrixXd clrContrast(int D){
  return coda::clrContrast(D);
}

//' Generalized LR Transform
//' @param X data (parts x samples aka D x N)
//' @param V contrast matrix (P x D)
// [[Rcpp::export]]
Eigen::MatrixXd glr(Eigen::Map<Eigen::MatrixXd>& X, Eigen::Map<Eigen::MatrixXd>& V){
  return coda::glr(X, V);
}

//' Generalized LR Inverse Transform
//' @param X data (coords x samples aka P x N)
//' @param V contrast matrix (P x D)
// [[Rcpp::export]]
Eigen::MatrixXd glrInv(Eigen::Map<Eigen::MatrixXd>& X, Eigen::Map<Eigen::MatrixXd>& V){
  return coda::glrInv(X,V);
}

//' Additive LR Transform
//' @param X data (parts x samples aka P x N)
//' @param d index of part to take as reference for ALR 
//'   (e.g., number of the part that goes in the denominator)
// [[Rcpp::export]]
Eigen::MatrixXd alr(Eigen::Map<Eigen::MatrixXd>& X, int d){
  return coda::alr(X,d);
}

//' Inverse Additive LR Transform
//' @param X data (coords x samples aka P x N)
//' @param d index of part to take as reference for ALR 
//'   (e.g., number of the part that goes in the denominator)
// [[Rcpp::export]]
Eigen::MatrixXd alrInv(Eigen::Map<Eigen::MatrixXd>& X, int d){
  return coda::alrInv(X,d);
}