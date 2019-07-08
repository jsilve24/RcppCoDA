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
  
//' Transfer data and covariance matricies from one coordinate system to another
//' 
//' Naming Conventions:
//' xxContrast produces a contrast matrix for transforming between representations where 
//' i=ILR, c=CLR, a=ALR. x2x represent transformations of data between coordinate systems 
//' xvar2xvar represent transformations of covariance matricies between coordinate systems. 
//' By convention all functions take arguments as follows: all functions have the form function(a,b,c)
//' where a is the data/covariance, b is the argument for specifying the current coordinate system, and c 
//' represents the argument for specifying the desired coordinate system. 
//' @param X data to transform between representations (parts/coords x samples)
//' @param Sigma covariance matrix in specified transformed space
//' @param V ILR contrast matrix (i.e., transformation matrix of ILR)
//' @param V1 ILR contrast matrix of basis Sigma is already in
//' @param V2 ILR contrast matrix of basis Sigma is desired in
//' @param d1 alr reference element Sigma is already expressed with respec to
//' @param d2 alr reference element Sigma is to be expressed with respect to
//' @return matrix
//' @name convert_coda
//' @export
//' @examples
//' # For Data
//' X <- matrix(abs(rnorm(10)), 5, 2)
//' X <- clo(X)
//' X.clr <- clr(X)
//' X.alr <- clr2alr(X, 3)
//' 
//' # For a covariance matrix starting in the alr
//' Sigma.alr <- X.alr%*%t(X.alr) # Just a random covariance matrix 
//' Sigma.ilr <- alrvar2ilrvar(Sigma.alr, 3, ilrContrast(5)) # Covert to default ILR
//' 
//' # Another way of doing this more "manually" - use contrast matricies
//' V <- caContrast(3, nrow(X)) # create contrast from clr to alr
//' X.alr <- V %*% X.clr
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

//' @rdname convert_coda
// [[Rcpp::export]]
Eigen::MatrixXd ilr2ilr_internal(Eigen::Map<Eigen::MatrixXd>& X, 
                        Eigen::Map<Eigen::MatrixXd>& V1, 
                        Eigen::Map<Eigen::MatrixXd>& V2){
  return coda::ilr2ilr(X, V1, V2);
}

//' @rdname convert_coda
// [[Rcpp::export]]
Eigen::MatrixXd ilr2clr_internal(Eigen::Map<Eigen::MatrixXd>& X, 
                        Eigen::Map<Eigen::MatrixXd>& V1){
  return coda::ilr2clr(X, V1);
}

//' @rdname convert_coda
// [[Rcpp::export]]
Eigen::MatrixXd clr2ilr_internal(Eigen::Map<Eigen::MatrixXd>& X, 
                        Eigen::Map<Eigen::MatrixXd>& V2){
  return coda::clr2ilr(X, V2);
}

//' @rdname convert_coda
// [[Rcpp::export]]
Eigen::MatrixXd alr2clr_internal(Eigen::Map<Eigen::MatrixXd>& X, int d1){
  return coda::alr2clr(X, d1);
}

//' @rdname convert_coda
// [[Rcpp::export]]
Eigen::MatrixXd clr2alr_internal(Eigen::Map<Eigen::MatrixXd>& X, int d2){
  return coda::clr2alr(X, d2);
}

//' @rdname convert_coda
// [[Rcpp::export]]
Eigen::MatrixXd alr2alr_internal(Eigen::Map<Eigen::MatrixXd>& X, int d1, int d2){
  return coda::alr2alr(X, d1, d2);
}

//' @rdname convert_coda
// [[Rcpp::export]]
Eigen::MatrixXd ilr2alr_internal(Eigen::Map<Eigen::MatrixXd>& X, 
                        Eigen::Map<Eigen::MatrixXd>& V1, int d2){
  return coda::ilr2alr(X, V1, d2);
}

//' @rdname convert_coda
// [[Rcpp::export]]
Eigen::MatrixXd alr2ilr_internal(Eigen::Map<Eigen::MatrixXd>& X, 
                        int d1, Eigen::Map<Eigen::MatrixXd>& V2){
  return coda::alr2ilr(X, d1, V2);
}


//// - covariances

//' @rdname convert_coda
// [[Rcpp::export]]
Eigen::MatrixXd ilrvar2ilrvar_internal(Eigen::Map<Eigen::MatrixXd>& Sigma, 
                        Eigen::Map<Eigen::MatrixXd>& V1, 
                        Eigen::Map<Eigen::MatrixXd>& V2){
  return coda::ilrvar2ilrvar(Sigma, V1, V2);
}

//' @rdname convert_coda
// [[Rcpp::export]]
Eigen::MatrixXd ilrvar2clrvar_internal(Eigen::Map<Eigen::MatrixXd>& Sigma, 
                        Eigen::Map<Eigen::MatrixXd>& V1){
  return coda::ilrvar2clrvar(Sigma, V1);
}

//' @rdname convert_coda
// [[Rcpp::export]]
Eigen::MatrixXd clrvar2ilrvar_internal(Eigen::Map<Eigen::MatrixXd>& Sigma, 
                        Eigen::Map<Eigen::MatrixXd>& V2){
  return coda::clrvar2ilrvar(Sigma, V2);
}

//' @rdname convert_coda
// [[Rcpp::export]]
Eigen::MatrixXd alrvar2clrvar_internal(Eigen::Map<Eigen::MatrixXd>& Sigma, int d1){
  return coda::alrvar2clrvar(Sigma, d1);
}

//' @rdname convert_coda
// [[Rcpp::export]]
Eigen::MatrixXd clrvar2alrvar_internal(Eigen::Map<Eigen::MatrixXd>& Sigma, int d2){
  return coda::clrvar2alrvar(Sigma, d2);
}

//' @rdname convert_coda
// [[Rcpp::export]]
Eigen::MatrixXd alrvar2alrvar_internal(Eigen::Map<Eigen::MatrixXd>& Sigma, int d1, int d2){
  return coda::alrvar2alrvar(Sigma, d1, d2);
}

//' @rdname convert_coda
// [[Rcpp::export]]
Eigen::MatrixXd ilrvar2alrvar_internal(Eigen::Map<Eigen::MatrixXd>& Sigma, 
                        Eigen::Map<Eigen::MatrixXd>& V1, int d2){
  return coda::ilrvar2alrvar(Sigma, V1, d2);
}

//' @rdname convert_coda
// [[Rcpp::export]]
Eigen::MatrixXd alrvar2ilrvar_internal(Eigen::Map<Eigen::MatrixXd>& Sigma, 
                        int d1, Eigen::Map<Eigen::MatrixXd>& V2){
  return coda::alrvar2ilrvar(Sigma, d1, V2);
}
 
