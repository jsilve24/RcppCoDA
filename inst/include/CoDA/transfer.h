#ifndef CODA__TRANSFER__H
#define CODA__TRANSFER__H

#include "RcppEigen.h"
#include "contrast.h"

// Convention here is that samples are columns and parts are rows
// This is notably different than how Driver was originally implemented. Sorry
// for any confusion.
// 
// In contrast to the GLR Functions, none of the functions in this file
// accept "block" data / covariance matricies. That is, they assume 
// that there is a single composition and transform, not a blocked transform
// as the GLR functions allow for. 

namespace coda {

  // For all of these functions *1 refers to the parameters defining the current
  // coordinate system, *2 refers to parameters defining the desired coordinate
  // system (e.g., V1 is the contrast matrix of the current ILR; V2 is the 
  // contrast matrix of the desired ILR). 

  template <typename TX, typename TV1, typename TV2>
  Eigen::MatrixXd ilr2ilr(Eigen::MatrixBase<TX>& X, 
                          Eigen::MatrixBase<TV1>& V1, 
                          Eigen::MatrixBase<TV2>& V2){
    MatrixXd V = iiContrast(V1, V2);
    MatrixXd Y = linForm(X, V);
    return Y;
  }
  
  template <typename TX, typename TV1>
  Eigen::MatrixXd ilr2clr(Eigen::MatrixBase<TX>& X, 
                          Eigen::MatrixBase<TV1>& V1){
    MatrixXd V = icContrast(V1);
    MatrixXd Y = linForm(X, V);
    return Y;
  }
  
  template <typename TX, typename TV2>
  Eigen::MatrixXd clr2ilr(Eigen::MatrixBase<TX>& X, 
                          Eigen::MatrixBase<TV2>& V2){
    MatrixXd V = ciContrast(V2);
    MatrixXd Y = linForm(X, V);
    return Y;
  }
  
  template <typename TX>
  Eigen::MatrixXd alr2clr(Eigen::MatrixBase<TX>& X, int d1){
    MatrixXd V = acContrast(d1, X.rows()+1);
    MatrixXd Y = linForm(X, V);
    return Y;
  }
  
  template <typename TX>
  Eigen::MatrixXd clr2alr(Eigen::MatrixBase<TX>& X, int d2){
    MatrixXd V = caContrast(d2, X.rows());
    MatrixXd Y = linForm(X, V);
    return Y;
  }
  
  template <typename TX>
  Eigen::MatrixXd alr2alr(Eigen::MatrixBase<TX>& X, int d1, int d2){
    MatrixXd V = aaContrast(d1, d2, X.rows()+1);
    MatrixXd Y = linForm(X, V);
    return Y;
  }
  
  template <typename TX, typename TV>
  Eigen::MatrixXd ilr2alr(Eigen::MatrixBase<TX>& X, 
                          Eigen::MatrixBase<TV>& V1, int d2){
    MatrixXd V = iaContrast(V1, d2, X.rows()+1);
    MatrixXd Y = linForm(X, V);
    return Y;
  }
  
  template <typename TX, typename TV>
  Eigen::MatrixXd alr2ilr(Eigen::MatrixBase<TX>& X,int d1, Eigen::MatrixBase<TV>& V2){
    MatrixXd V = aiContrast(d1, V2, X.rows()+1);
    MatrixXd Y = linForm(X, V);
    return Y;
  }
  
  
  template <typename TS, typename TV1, typename TV2>
  Eigen::MatrixXd ilrvar2ilrvar(Eigen::MatrixBase<TS>& X, 
                                Eigen::MatrixBase<TV1>& V1, 
                                Eigen::MatrixBase<TV2>& V2){
    MatrixXd V = iiContrast(V1, V2);
    MatrixXd O = quadForm(X, V);
    return O;
  }
  
  template <typename TX, typename TV1>
  Eigen::MatrixXd ilrvar2clrvar(Eigen::MatrixBase<TX>& X, 
                                Eigen::MatrixBase<TV1>& V1){
    MatrixXd V = icContrast(V1);
    MatrixXd O = quadForm(X, V);
    return O;
  }
  
  template <typename TX, typename TV2>
  Eigen::MatrixXd clrvar2ilrvar(Eigen::MatrixBase<TX>& X, 
                                Eigen::MatrixBase<TV2>& V2){
    MatrixXd V = ciContrast(V2);
    MatrixXd O = quadForm(X, V);
    return O;
  }
  
  template <typename TX>
  Eigen::MatrixXd alrvar2clrvar(Eigen::MatrixBase<TX>& X, int d1){
    MatrixXd V = acContrast(d1, X.rows()+1);
    MatrixXd O = quadForm(X, V);
    return O;
  }
  
  template <typename TX>
  Eigen::MatrixXd clrvar2alrvar(Eigen::MatrixBase<TX>& X, int d2){
    MatrixXd V = caContrast(d2, X.rows());
    MatrixXd O = quadForm(X, V);
    return O;
  }
  
  template <typename TX>
  Eigen::MatrixXd alrvar2alrvar(Eigen::MatrixBase<TX>& X, int d1, int d2){
    MatrixXd V = aaContrast(d1, d2, X.rows()+1);
    MatrixXd O = quadForm(X, V);
    return O;
  }
  
  template <typename TX, typename TV>
  Eigen::MatrixXd ilrvar2alrvar(Eigen::MatrixBase<TX>& X, 
                                Eigen::MatrixBase<TV>& V1, int d2){
    MatrixXd V = iaContrast(V1, d2, X.rows()+1);
    MatrixXd O = quadForm(X, V);
    return O;
  }
  
  template <typename TX, typename TV>
  Eigen::MatrixXd alrvar2ilrvar(Eigen::MatrixBase<TX>& X,int d1, 
                                Eigen::MatrixBase<TV>& V2){
    MatrixXd V = aiContrast(d1, V2, X.rows()+1);
    MatrixXd O = quadForm(X, V);
    return O;
  }
}

#endif 