#ifndef CODA__CODAMATH__H
#define CODA__CODAMATH__H

#include <RcppEigen.h>

using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::Map;
using Eigen::MatrixBase;
using Eigen::SparseMatrix;

namespace coda{

  //' Closure divide each column by column sum
  //' @param X matrix
  //' @return matrix of same dimension as X
  template <typename T>
  Eigen::MatrixXd clo(Eigen::DenseBase<T>& X){
    VectorXd s = X.colwise().sum();
    int c = X.cols();
    MatrixXd Y(X.rows(), X.cols());
    for (int i=0; i<c; i++) Y.col(i) = X.col(i)/s(i);
    return Y;
  }
  
  //' Center opteration - subtract each column by column mean
  //' @param X matrix
  //' @return matrix of same dimension as X
  template <typename T>
  Eigen::MatrixXd center(Eigen::DenseBase<T>& X){
    VectorXd n = X.colwise().mean();
    X.rowwise() -= n.transpose();
    return X;
  }
  
  //' Linear form Y=VX
  //' @param X data (D x N)
  //' @param V1 (P1 X D1, e.g., transfer Contrast matrix)
  //' @param V2 (P2 x D2, e.g., transfer Contrast matrix)
  //' @details Calculates 
  //'      
  //'   Y = bdiag[V1, V2] * X
  //' where bdiag a block diagonal matrix with blocks V1 and V2. 
  //' If V2 is not given then assumes that second block is identity. 
  //' In that case determines the size of second block based on size of 1. 
  //' 
  //' @return Eigen::MatrixXd 
  //' @name linForm
  template <typename TX, typename TV1>
  Eigen::MatrixXd linForm(Eigen::MatrixBase<TX>& X, 
                               Eigen::MatrixBase<TV1>& V1){
    int P1=V1.rows();
    int D1=V1.cols();
    int D = X.rows();
    int N = X.cols();
    if (X.rows() < D1) throw std::invalid_argument("X.rows() >= V1.cols()");
    
    if (X.rows() == D1) {
      return V1*X;
    }
    
    // else 
    MatrixXd O(D-D1+P1, N);
    O.topRows(P1) = V1*X.topRows(D1);
    O.bottomRows(D-D1) = X.bottomRows(D-D1);  
    return O;
  }
  
  //' @rdname linForm
  template <typename TX, typename TV1, typename TV2>
  Eigen::MatrixXd linForm(Eigen::MatrixBase<TX>& X, 
                          Eigen::MatrixBase<TV1>& V1, 
                          Eigen::MatrixBase<TV2>& V2){
    int P1=V1.rows();
    int D1=V1.cols();
    int P2 = V2.rows();
    int D2 = V2.cols();
    int N = X.cols();
    if (X.rows() != D1 + D2) throw std::invalid_argument("X.rows() must equal V1.cols() + V2.cols()");
    MatrixXd O(P2+P1, N);
    O.topRows(D1) = V1*X.topRows(P1);
    O.bottomRows(D2) = V2*X.bottomRows(P2);
    return O;
  }
  
  
  //' Quadratic Form Y = V X V^T
  //' @param X Dx(DN) matrix (e.g., N column bound square matricies)
  //' @param V1 (e.g., transfer Contrast matrix)
  //' @param V2 (e.g., transfer Contrast matrix)
  //' @return Eigen::MatrixXd
  //' @details Calculates:
  //' 
  //' Y_i = bdiag[V1, V2] * X_i * bdiag[V1, V2]^T
  //' 
  //' where Y_i and X_i refer to the i-th set of D-columns of Y and X respectively
  //' and where bdiag a block diagonal matrix with blocks V1 and V2. 
  //' If V2 is not given then assumes that second block is identity. 
  //' In that case determines the size of second block based on size of V1.
  //' 
  //' @name quadForm
  template <typename TX, typename TV1>
  Eigen::MatrixXd quadForm(Eigen::MatrixBase<TX>& X, 
                                     Eigen::MatrixBase<TV1>& V1){
    int P1=V1.rows();
    int D1=V1.cols();
    int D = X.rows();
    int D2 = D-D1;
    int P2 = D2;
    int P = P1+P2;
    int N = X.cols();
    if ( (N % D) != 0 ) throw std::invalid_argument("X must be Dx(DN) see documentation");
    if (N == 0 ) throw std::invalid_argument("X must have columns");
    N = N/D; // safe after above validation
    if (D < D1) throw std::invalid_argument("X row/column dimension must be larger than V1.cols()");
    
    if ( (D == D1) & (N==1) ){
      return V1*X*V1.transpose();
    } 
    
    // else N > 0
    MatrixXd O(P, N*P);
    if ( D == D1 ){ 
      for (int i=0; i<N; i++){
        O.middleCols(i*P1, P1) = V1*X.middleCols(i*D, D)*V1.transpose();
      }
    }
    
    // else
    for (int i=0; i<N; i++){
      O.block(0,i*P, P1, P1).noalias() = V1*X.block(0, i*D, D1, D1)*V1.transpose();
      O.block(0,i*P+P1, P1, P2).noalias() = V1*X.block(0, i*D+D1, D1, D2);
      O.block(P1,i*P, P2, P1) = O.block(0,i*P+P1, P1, P2).transpose();
      O.block(P1,i*P+P1, P2, P2) = X.block(D1, i*D+D1, D2, D2);
    }
    return O;
  }
  
  //' @rdname quadForm
  template <typename TX, typename TV1, typename TV2>
  Eigen::MatrixXd quadForm(Eigen::MatrixBase<TX>& X, 
                           Eigen::MatrixBase<TV1>& V1, 
                           Eigen::MatrixBase<TV2>& V2){
    int P1=V1.rows();
    int D1=V1.cols();
    int P2=V2.rows();
    int D2=V2.cols();
    int P = P1+P2;
    int D = X.rows();

    int N = X.cols();
    if ( (N % D) != 0 ) throw std::invalid_argument("X must be Dx(DN) see documentation");
    if (N == 0 ) throw std::invalid_argument("X must have columns");
    N = N/D; // safe after above validation
    if (D != D1+D2) throw std::invalid_argument("X.rows() must equal V1.cols() + V2.cols()");
    
    MatrixXd O(P, N*P);
    for (int i=0; i<N; i++){
      O.block(0,i*P, P1, P1).noalias() = V1*X.block(0, i*D, D1, D1)*V1.transpose();
      O.block(0,i*P+P1, P1, P2).noalias() = V1*X.block(0, i*D+D1, D1, D2)*V2.transpose();
      O.block(P1,i*P, P2, P1) = O.block(0,i*P+P1, P1, P2).transpose();
      O.block(P1,i*P+P1, P2, P2) = V2*X.block(D1, i*D+D1, D2, D2)*V2.transpose();
    }
    return O;
  }

}

#endif 