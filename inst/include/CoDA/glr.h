#ifndef CODA__GLR__H
#define CODA__GLR__H

#include "codamath.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::Map;
using Eigen::MatrixBase;
using Eigen::SparseMatrix;

// Convention here is that samples are columns and parts are rows
// This is notably different than how Driver was originally implemented. Sorry
// for any confusion.

namespace coda {

  //' Generalized LR Transform
  //' @param X data (D x N) (e.g., parts by samples)
  //' @param V (P X D, e.g., transfer Contrast matrix)
  //' @return Eigen::MatrixXd
  //' @details calculates
  //' 
  //'   Calculates Y = GLR(X) or Y=GLRINV(X) defined by 
  //' @name glr
  template <typename TX, typename TV>
  Eigen::MatrixXd glr(Eigen::MatrixBase<TX>& X, Eigen::MatrixBase<TV>& V){ // was MatrixBase and SparseMatrixBase
    int P = V.rows();
    int D = V.cols();
    if (X.rows() != D) throw std::invalid_argument("X.rows() != V.cols()");
    
    MatrixXd Y = X.array().log().matrix();
    return V*Y;
  }

  //' @rdname glr
  template <typename TX, typename TV>
  Eigen::MatrixXd glr(Eigen::MatrixBase<TX>& X, Eigen::SparseMatrixBase<TV>& V){ // was MatrixBase and SparseMatrixBase
    int P = V.rows();
    int D = V.cols();
    if (X.rows() != D) throw std::invalid_argument("X.rows() != V.cols()");
    
    MatrixXd Y = X.array().log().matrix();
    return V*Y;
  }


  //' @rdname glr
  template <typename TX, typename TV>
  Eigen::MatrixXd glrInv(Eigen::MatrixBase<TX>& X, Eigen::MatrixBase<TV>& V){
    int P = V.rows();
    int D = V.cols();
    if (X.rows() != P) throw std::invalid_argument("X.rows() != V.rows()");
    
    MatrixXd O;
    O.noalias() = V.transpose()*X;
    O = O.array().exp().matrix();
    return clo(O);
  }


  //' @rdname glr
  template <typename TX, typename TV>
  Eigen::MatrixXd glrInv(Eigen::MatrixBase<TX>& X, Eigen::SparseMatrixBase<TV>& V){
    int P = V.rows();
    int D = V.cols();
    if (X.rows() != P) throw std::invalid_argument("X.rows() != V.cols()");
    
    MatrixXd O;
    O.noalias() = V.transpose()*X;
    O = O.array().exp().matrix();
    return clo(O);
  }

} /* End coda Namespace */




#endif 
