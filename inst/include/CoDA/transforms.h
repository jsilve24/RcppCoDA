#ifndef CODA__TRANSFORMS__H
#define CODA__TRANSFORMS__H

#include <RcppEigen.h>
#include "contrast.h"
#include "glr.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::Map;
using Eigen::MatrixBase;
using Eigen::SparseMatrix;


// Convention here is that samples are columns and parts are rows
// This is notably different than how Driver was originally implemented. Sorry
// for any confusion.

namespace coda {

  //' Additive LR Transform
  //' @param X data (parts x samples aka D x N)
  //' @param d index of part to take as reference for ALR
  //'   (e.g., number of the part that goes in the denominator)
  template <typename T>
  Eigen::MatrixXd alr(Eigen::MatrixBase<T>& X, int d){
    SparseMatrix<double> B = alrContrastSparse(d, X.rows(), false);
    return glr(X, B);
  }
  
  //' Inverse Additive LR Transform
  //' @param X data (coords x samples aka P x N)
  //' @param d index of part to take as reference for ALR
  //'   (e.g., number of the part that goes in the denominator)
  template <typename T>
  Eigen::MatrixXd alrInv(Eigen::MatrixBase<T>& X, int d){
    SparseMatrix<double> B = alrContrastSparse(d, X.rows()+1, true);
    return glrInv(X, B);
  }
  
  //' Centered LR Transform
  //' @param X data (parts x samples aka D x N)
  template <typename T>
  Eigen::MatrixXd clr(Eigen::MatrixBase<T>& X){
    MatrixXd Y = X.array().log().matrix();
    return center(Y);
  }
  
  //' Inverse Centered LR Transform
  //' @param X data (coords x samples aka P x N) here P = D because its the CLR
  template <typename T>
  Eigen::MatrixXd clrInv(Eigen::MatrixBase<T>& X){
    MatrixXd Y = X.array().exp().matrix();
    return clo(Y);
  }
  
  //' Isometric LR Transform
  //' @param X data (parts x samples aka D x N)
  //' @param V contrast matrix (P x D)
  template <typename TX, typename TV>
  Eigen::MatrixXd ilr(Eigen::MatrixBase<TX>& X, Eigen::MatrixBase<TV>& V){
    return glr(X, V);
  }
  
  //' Inverse Isometric LR Transform
  //' @param X data (coords x samples aka P x N) 
  //' @param V contrast matrix (P x D)
  template <typename TX, typename TV>
  Eigen::MatrixXd ilrInv(Eigen::MatrixBase<TX>& X, Eigen::MatrixBase<TV>& V){
    return glrInv(X, V);
  }
  
  //' Isometric LR Transform - default basis
  //' @param X data (parts x samples aka D x N)
  template <typename TX>
  Eigen::MatrixXd ilr(Eigen::MatrixBase<TX>& X){
    MatrixXd V = ilrContrast(X.rows());
    return glr(X, V);
  }
  
  //' Inverse Isometric LR Transform - default basis
  //' @param X data (coords x samples aka P x N) 
  template <typename TX>
  Eigen::MatrixXd ilrInv(Eigen::MatrixBase<TX>& X){
    MatrixXd V = ilrContrast(X.rows()+1);
    return glrInv(X, V);
  }
  
}

#endif 
