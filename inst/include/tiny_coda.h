#ifndef RCPPCODA__CODA__H
#define RCPPCODA__CODA__H

#include <RcppEigen.h>

using namespace Rcpp;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::Map;
using Eigen::MatrixBase;
using Eigen::SparseMatrix;


// Convention here is that samples are columns and parts are rows
// This is notably different than how Driver was originally implemented. Sorry
// for any confusion.

namespace coda {

  //' Closure divide each column by column sum
  //' @param X matrix
  //' @return matrix of same dimension as X
  template <typename T>
  Eigen::MatrixXd miniclo(Eigen::DenseBase<T>& X){
    VectorXd s = X.colwise().sum();
    int c = X.cols();
    MatrixXd Y(X.rows(), X.cols());
    for (int i=0; i<c; i++) Y.col(i) = X.col(i)/s(i);
    return Y;
  }
  
  
  //' Create alr contrast matrix
  //' @param D total number of parts
  //' @param d reference part for alr
  //' @param inv (is this being used for the alr or alrInv?)
  //' @return D-1 x D matrix
  Eigen::SparseMatrix<double> alrContrastSparse(int D, int d, bool inv){
    if ((d<1) || (d>D)) throw std::invalid_argument("d must be within [1,...,D]");
    SparseMatrix<double> B(D-1, D);// = MatrixXd::Zero(D-1, D);
    B.reserve(Eigen::VectorXi::Constant(D, 2)); // reserve 2 non-zero per column
    int pos=0;
    for (int i=0; i<D; i++){
      if (i == (d-1)) pos++;
      B.insert(i, pos) = 1;
      pos++;
    }
    if (!inv) {
      for (int i=0; i<(D-1); i++){
        B.insert(i,d-1) = -1;
      }
    }
    return B;
  }
  
  //' Create alr contrast matrix
  //' @param D total number of parts
  //' @param d reference part for alr
  //' @param inv (is this being used for the alr or alrInv?)
  //' @return D-1 x D matrix
  Eigen::MatrixXd alrContrast(int D, int d, bool inv){
    Eigen::MatrixXd B;
    SparseMatrix<double> Bsp = alrContrastSparse(D, d, inv);
    B = MatrixXd(Bsp);
    // if ((d<1) || (d>D)) throw std::invalid_argument("d must be within [1,...,D]");
    // MatrixXd B = MatrixXd::Zero(D-1, D);
    // int pos=0;
    // for (int i=0; i<D; i++){
    //   if (i == (d-1)) pos++;
    //   B(i, pos) = 1;
    //   pos++;
    // }
    // if (!inv) B.col(d-1) = -VectorXd::Ones(D);
    // //if (inv) B.col(d-1) = VectorXd::Zero(D);
    return B;
  }
  
  
  //' Create ilr default contrast matrix
  //' @param D total number of parts
  //' @return D-1xD contrast matrix (a basis)
  Eigen::MatrixXd ilrContrast(int D){
    if (D < 2) throw std::invalid_argument("D must be > 1");
    //Eigen::FullPivHouseholderQR<MatrixXd> qr(D-1, D);
    //qr.compute(alrContrast(D,D,false));
    MatrixXd B = alrContrast(D, D, false);
    MatrixXd Q = B.transpose().householderQr().householderQ().transpose();
    return Q.topRows(D-1);
  }
  
  //' Create CLR contrast matrix
  //' @param D total number of parts
  //' @return D x D contrast matrix 
  Eigen::MatrixXd clrContrast(int D){
    if (D < 2) throw std::invalid_argument("D must be >1");
    MatrixXd B = -MatrixXd::Ones(D,D);
    B.diagonal().array() += D;
    return B.array()/D;
  }
  
  
  //' Generalized LR Transform
  //' @param X data (parts x samples aka D x N)
  //' @param V contrast matrix (P x D)
  template <typename TX, typename TV>
  Eigen::MatrixXd glr(Eigen::MatrixBase<TX>& X, Eigen::MatrixBase<TV>& V){
    MatrixXd Y = X.array().log().matrix();
    return V*Y;
  }
  
  //' Generalized LR Inverse Transform
  //' @param X data (coords x samples aka P x N)
  //' @param V contrast matrix (P x D)
  template <typename TX, typename TV>
  Eigen::MatrixXd glrInv(Eigen::MatrixBase<TX>& X, Eigen::MatrixBase<TV>& V){
    MatrixXd Y = V.transpose()*X;
    Y = Y.array().exp().matrix();
    return miniclo(Y);
  }
  
  
  //' Generalized LR Transform
  //' @param X data (parts x samples aka D x N)
  //' @param V contrast matrix (P x D)
  template <typename TX, typename TV>
  Eigen::MatrixXd glr(Eigen::MatrixBase<TX>& X, Eigen::SparseMatrixBase<TV>& V){
    MatrixXd Y = X.array().log().matrix();
    return V*Y;
  }
  
  //' Generalized LR Inverse Transform
  //' @param X data (coords x samples aka P x N)
  //' @param V contrast matrix (P x D)
  template <typename TX, typename TV>
  Eigen::MatrixXd glrInv(Eigen::MatrixBase<TX>& X, Eigen::SparseMatrixBase<TV>& V){
    MatrixXd Y = V.transpose()*X;
    Y = Y.array().exp().matrix();
    return miniclo(Y);
  }

  //' Additive LR Transform
  //' @param X data (parts x samples aka D x N)
  //' @param d index of part to take as reference for ALR
  //'   (e.g., number of the part that goes in the denominator)
  template <typename T>
  Eigen::MatrixXd alr(Eigen::MatrixBase<T>& X, int d){
    SparseMatrix<double> B = alrContrastSparse(X.rows(), d, false);
    return glr(X, B);
  }
  
  //' Inverse Additive LR Transform
  //' @param X data (coords x samples aka P x N)
  //' @param d index of part to take as reference for ALR
  //'   (e.g., number of the part that goes in the denominator)
  template <typename T>
  Eigen::MatrixXd alrInv(Eigen::MatrixBase<T>& X, int d){
    SparseMatrix<double> B = alrContrastSparse(X.rows()+1, d, true);
    return glrInv(X, B);
  }
  
}



#endif
