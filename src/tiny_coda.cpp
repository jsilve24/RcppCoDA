#include "CoDA.h"

using namespace Rcpp;

using namespace Rcpp;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::Map;

// ----- CODA MATH

// [[Rcpp::export]]
Eigen::MatrixXd clo_internal(Eigen::Map<Eigen::MatrixXd> X){
  return coda::clo(X);
}

// [[Rcpp::export]]
Eigen::MatrixXd center_internal(Eigen::Map<Eigen::MatrixXd> X){
  return coda::center(X);
}


// [[Rcpp::export]]
Eigen::MatrixXd linForm_internal(Eigen::Map<Eigen::MatrixXd> X,
                        Eigen::Map<Eigen::MatrixXd> V1,
                        Rcpp::Nullable<Eigen::Map<Eigen::MatrixXd> > V2 = R_NilValue){
  if (V2.isNotNull()){
   Eigen::Map<Eigen::MatrixXd> V2tmp = Rcpp::as<Eigen::Map<Eigen::MatrixXd> >(V2);
   return coda::linForm(X, V1, V2tmp);
  }
  return coda::linForm(X, V1);
}


// [[Rcpp::export]]
Eigen::MatrixXd quadForm_internal(Eigen::Map<Eigen::MatrixXd> X,
                         Eigen::Map<Eigen::MatrixXd> V1,
                         Rcpp::Nullable<Eigen::Map<Eigen::MatrixXd> > V2 = R_NilValue){
  if (V2.isNotNull()){
    Eigen::Map<Eigen::MatrixXd> V2tmp = Rcpp::as<Eigen::Map<Eigen::MatrixXd> >(V2);
    return coda::quadForm(X, V1, V2tmp);
  }
  return coda::quadForm(X, V1);
}


// ------- GLR TRANSFORMS


// [[Rcpp::export]]
Eigen::MatrixXd glr_internal(Eigen::Map<Eigen::MatrixXd>& X, 
                             Eigen::Map<Eigen::MatrixXd>& V1, 
                             Rcpp::Nullable<Eigen::Map<Eigen::MatrixXd> > V2 = R_NilValue){
  if (V2.isNotNull()){
    Eigen::Map<Eigen::MatrixXd> V2tmp = Rcpp::as<Eigen::Map<Eigen::MatrixXd> >(V2);
    return coda::glr(X, V1, V2tmp);
  }
  return coda::glr(X, V1);
}

// [[Rcpp::export]]
Eigen::MatrixXd glrInv_internal(Eigen::Map<Eigen::MatrixXd>& X, 
                             Eigen::Map<Eigen::MatrixXd>& V1, 
                             Rcpp::Nullable<Eigen::Map<Eigen::MatrixXd> > V2 = R_NilValue){
  if (V2.isNotNull()){
    Eigen::Map<Eigen::MatrixXd> V2tmp = Rcpp::as<Eigen::Map<Eigen::MatrixXd> >(V2);
    return coda::glrInv(X, V1, V2tmp);
  }
  return coda::glrInv(X, V1);
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

//' @rdname base_lr_transforms
//' @export
// [[Rcpp::export]]
Eigen::MatrixXd ilrContrast(int D){
  return coda::ilrContrast(D);
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
  
//' @rdname convert_coda
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

// [[Rcpp::export]]
Eigen::MatrixXd ilr2ilr_internal(Eigen::Map<Eigen::MatrixXd>& X, 
                        Eigen::Map<Eigen::MatrixXd>& V1, 
                        Eigen::Map<Eigen::MatrixXd>& V2){
  return coda::ilr2ilr(X, V1, V2);
}

// [[Rcpp::export]]
Eigen::MatrixXd ilr2clr_internal(Eigen::Map<Eigen::MatrixXd>& X, 
                        Eigen::Map<Eigen::MatrixXd>& V1){
  return coda::ilr2clr(X, V1);
}

// [[Rcpp::export]]
Eigen::MatrixXd clr2ilr_internal(Eigen::Map<Eigen::MatrixXd>& X, 
                        Eigen::Map<Eigen::MatrixXd>& V2){
  return coda::clr2ilr(X, V2);
}

// [[Rcpp::export]]
Eigen::MatrixXd alr2clr_internal(Eigen::Map<Eigen::MatrixXd>& X, int d1){
  return coda::alr2clr(X, d1);
}

// [[Rcpp::export]]
Eigen::MatrixXd clr2alr_internal(Eigen::Map<Eigen::MatrixXd>& X, int d2){
  return coda::clr2alr(X, d2);
}

// [[Rcpp::export]]
Eigen::MatrixXd alr2alr_internal(Eigen::Map<Eigen::MatrixXd>& X, int d1, int d2){
  return coda::alr2alr(X, d1, d2);
}

// [[Rcpp::export]]
Eigen::MatrixXd ilr2alr_internal(Eigen::Map<Eigen::MatrixXd>& X, 
                        Eigen::Map<Eigen::MatrixXd>& V1, int d2){
  return coda::ilr2alr(X, V1, d2);
}

// [[Rcpp::export]]
Eigen::MatrixXd alr2ilr_internal(Eigen::Map<Eigen::MatrixXd>& X, 
                        int d1, Eigen::Map<Eigen::MatrixXd>& V2){
  return coda::alr2ilr(X, d1, V2);
}


//// - covariances

// [[Rcpp::export]]
Eigen::MatrixXd ilrvar2ilrvar_internal(Eigen::Map<Eigen::MatrixXd>& Sigma, 
                        Eigen::Map<Eigen::MatrixXd>& V1, 
                        Eigen::Map<Eigen::MatrixXd>& V2){
  return coda::ilrvar2ilrvar(Sigma, V1, V2);
}

// [[Rcpp::export]]
Eigen::MatrixXd ilrvar2clrvar_internal(Eigen::Map<Eigen::MatrixXd>& Sigma, 
                        Eigen::Map<Eigen::MatrixXd>& V1){
  return coda::ilrvar2clrvar(Sigma, V1);
}

// [[Rcpp::export]]
Eigen::MatrixXd clrvar2ilrvar_internal(Eigen::Map<Eigen::MatrixXd>& Sigma, 
                        Eigen::Map<Eigen::MatrixXd>& V2){
  return coda::clrvar2ilrvar(Sigma, V2);
}

// [[Rcpp::export]]
Eigen::MatrixXd alrvar2clrvar_internal(Eigen::Map<Eigen::MatrixXd>& Sigma, int d1){
  return coda::alrvar2clrvar(Sigma, d1);
}

// [[Rcpp::export]]
Eigen::MatrixXd clrvar2alrvar_internal(Eigen::Map<Eigen::MatrixXd>& Sigma, int d2){
  return coda::clrvar2alrvar(Sigma, d2);
}

// [[Rcpp::export]]
Eigen::MatrixXd alrvar2alrvar_internal(Eigen::Map<Eigen::MatrixXd>& Sigma, int d1, int d2){
  return coda::alrvar2alrvar(Sigma, d1, d2);
}

// [[Rcpp::export]]
Eigen::MatrixXd ilrvar2alrvar_internal(Eigen::Map<Eigen::MatrixXd>& Sigma, 
                        Eigen::Map<Eigen::MatrixXd>& V1, int d2){
  return coda::ilrvar2alrvar(Sigma, V1, d2);
}

// [[Rcpp::export]]
Eigen::MatrixXd alrvar2ilrvar_internal(Eigen::Map<Eigen::MatrixXd>& Sigma, 
                        int d1, Eigen::Map<Eigen::MatrixXd>& V2){
  return coda::alrvar2ilrvar(Sigma, d1, V2);
}

// [[Rcpp::export]]
Eigen::MatrixXd clrvar2vararray_internal(Eigen::Map<Eigen::MatrixXd>& Sigma){
  return coda::clrvar2vararray(Sigma);
}

//' Functions for implementing / calculating IQLR 
//' 
//' var2iqlrContrast calculates contrast matrix in iqlr based on part-wise variance 
//' and a specification for the lower and upper quantiles to use. 
//' 
//' clrvar2iqlrvar transforms covariance matrices into the iqlr. Note a different 
//' iqlr is calculated for each covariance matrix
//' 
//' These are advanced features and should not be used without prior knowledge of the
//' IQLR transform. 
//' 
//' @param S P-vector of "variances" for IQLR (e.g., diagonal of CLR covariance matrix)
//' @param Sigma (PN)xP covariance matrix in CLR
//' @param qLow lower quantile (between 0 and 1)
//' @param qHigh upper quantile (between 0 and 1; must be greater than qLow)
//' @return PxP contrast matrix (var2iqlrContrast) or 
//'    (PN)xP covariance matricies (clrvar2iqlrvar)
//' @name iqlr
//' @export
//' @references Jia R Wu, Jean M Macklaim, Briana L Genge, Gregory B Gloor. 
//'   Finding the centre: corrections for asymmetry in high-throughput sequencing datasets. 
//'   (2017) arXiv preprint arXiv:1704.01841
// [[Rcpp::export]]
Eigen::MatrixXd var2iqlrContrast(Eigen::Map<Eigen::VectorXd>& S, double qLow=0.25, double qHigh=0.75){
  return coda::var2iqlrContrast(S, qLow, qHigh);
}

// [[Rcpp::export]]
Eigen::MatrixXd clrvar2iqlrvar_internal(Eigen::Map<Eigen::MatrixXd>& Sigma, double qLow=0.25, double qHigh=0.75){
  return coda::clrvar2iqlrvar(Sigma, qLow, qHigh);
}


