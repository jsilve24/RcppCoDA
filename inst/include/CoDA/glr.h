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

namespace coda{

  //' Generalized LR Transform
  //' @param X data (D x N) (e.g., parts by samples)
  //' @param V1 (P1 X D1, e.g., transfer Contrast matrix)
  //' @param V2 (P2 x D2, e.g., transfer Contrast matrix)
  //' @return Eigen::MatrixXd
  //' @details calculates
  //' 
  //'   let Y = cbind[Y1, Y2] 
  //'   let X = cbind[X1, X2] 
  //'   Calculates Y = GLR(X) or Y=GLRINV(X) defined by 
  //'   Y1 = GLR(X1, V1) or Y1 = GLRINV(X1, V1)
  //'   Y2 = GLR(X2, V2) or Y2 = GLRINV(X2, V2)
  //'   
  //' If V2 is not given then assumes that Y2 = X2. 
  //' In that case determines the size of X2 and Y2 based on on size of V1. 
  //' @name glr
  template <typename TX, typename TV1>
  Eigen::MatrixXd glr(Eigen::MatrixBase<TX>& X, Eigen::EigenBase<TV1>& V1){ // was MatrixBase and SparseMatrixBase
    int P1 = V1.rows();
    int D1 = V1.cols();
    int D = X.rows();
    int D2 = D-D1;
    int P2 = D2;
    int N = X.cols();
    if (X.rows() < D1) throw std::invalid_argument("X.rows() >= V1.cols()");
    
    if (D1==D){
      MatrixXd Y = X.array().log().matrix();
      return V1*Y;
    }
    
    // else 
    MatrixXd O(P1+P2, N);
    MatrixXd Y = X.topRows(D1).array().log().matrix();
    O.topRows(P1).noalias() = V1*Y;
    O.bottomRows(D-D1) = X.bottomRows(D-D1);
    return O;
  }
  
  //' @rdname glr
  template <typename TX, typename TV1, typename TV2>
  Eigen::MatrixXd glr(Eigen::MatrixBase<TX>& X, 
                      Eigen::EigenBase<TV1>& V1, 
                      Eigen::EigenBase<TV2>& V2){ // was MatrixBase and SparseMatrixBase
    int P1 = V1.rows();
    int D1 = V1.cols();
    int P2 = V2.rows();
    int D2 = V2.cols();
    int D = X.rows();
    int N = X.cols();
    if (X.rows() != D1+D2) throw std::invalid_argument("X.rows() must = V1.cols() + V2.cols()");
    
    // else 
    MatrixXd O(P1+P2, N);
    MatrixXd Y = X.array().log().matrix();
    O.topRows(P1).noalias() = V1*Y.topRows(D1);
    O.bottomRows(D-D1) = V2*Y.bottomRows(D-D1);
    return O;
  }
  
  
  //' @rdname glr
  template <typename TX, typename TV1>
  Eigen::MatrixXd glrInv(Eigen::MatrixBase<TX>& X, Eigen::EigenBase<TV1>& V1){
    int P1 = V1.rows();
    int D1 = V1.cols();
    int P = X.rows();
    int P2 = P-P1;
    int D2 = P2;
    int N = X.cols();
    if (X.rows() < P1) throw std::invalid_argument("X.rows() >= V1.rows()");
    
    if (P1 == P){
      MatrixXd O;
      O.noalias() = V1.transpose()*X;
      O = O.array().exp().matrix();
      return clo(O);
    }
    
    // else
    MatrixXd O(D1+D2, N);
    MatrixXd Y;
    Y.noalias() = V1.transpose()*X.topRows(P1);
    Y = Y.array().exp().matrix();
    O.topRows(D1) = clo(Y);
    O.bottomRows(D2) = X.bottomRows(D2);
    return O;
  }
  
  //' @rdname glr
  template <typename TX, typename TV1, typename TV2>
  Eigen::MatrixXd glrInv(Eigen::MatrixBase<TX>& X, Eigen::EigenBase<TV1>& V1, 
                         Eigen::EigenBase<TV2>& V2){
    int P1 = V1.rows();
    int D1 = V1.cols();
    int P2 = V2.rows();
    int D2 = V2.cols();
    int P = X.rows();
    int N = X.cols();
    if (X.rows() != P1+P2) throw std::invalid_argument("X.rows() must = V1.rows() + V2.rows()");
    
    // else
    MatrixXd O(D1+D2, N);
    MatrixXd Y1(D1, N);
    MatrixXd Y2(D2, N);
    Y1.noalias() = V1.transpose()*X.topRows(P1);
    Y2.noalias() = V2.transpose()*X.bottomRows(P2);
    Y1 = Y1.array().exp().matrix();
    Y2 = Y2.array().exp().matrix();
    O.topRows(D1) = clo(Y1);
    O.bottomRows(D2) = clo(Y2);
    return O;
  }

}

#endif 