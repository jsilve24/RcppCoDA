#ifndef CODA__PROPORTIONALITY__H
#define CODA__PROPORTIONALITY__H

#include <RcppEigen.h>

using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::Map;
using Eigen::MatrixBase;
using Eigen::SparseMatrix;


// Convention here is that samples are columns and parts are rows
// This is notably different than how Driver was originally implemented. Sorry
// for any confusion.

namespace coda {

  namespace internal {
    template <typename TS>
    Eigen::MatrixXd clrvar2vararray_single(Eigen::MatrixBase<TS>& Sigma){
      int P = Sigma.rows();
      MatrixXd res(P,P);
      if (Sigma.cols() != P) throw std::invalid_argument("Sigma must be PxP see documentation");
      for (int i=0; i<P; i++){
        for (int j=0; j<P; j++){
          res(i,j) = Sigma(i,i) + Sigma(j,j) - 2*Sigma(i,j);
        }
      }
      return res;
    }
    
    template <typename TS>
    Eigen::MatrixXd clrvar2phi_single(Eigen::MatrixBase<TS>& Sigma){
      int P = Sigma.rows();
      MatrixXd res(P,P);
      if (Sigma.cols() != P) throw std::invalid_argument("Sigma must be PxP see documentation");
      for (int i=0; i<P; i++){
        for (int j=0; j<P; j++){
          res(i,j) = Sigma(i,i) + Sigma(j,j) - 2*Sigma(i,j);
          res(i,j) = res(i,j)/Sigma(i,i);
        }
      }
      return res;
    }
    
  }
  
  //' Calculate Variation Array from CLR Covariances
  //' @param Sigma Covariance matrix Px(PN) where N is number of 
  //'   covariance matricies in CLR space
  //' See Lovell "Proportionality: A valid alternative ..." for details. 
  template <typename TS>
  Eigen::MatrixXd clrvar2vararray(Eigen::MatrixBase<TS>& Sigma){
    int P1 = Sigma.rows();
    int N = Sigma.cols();
    if ( (N % P1) != 0 ) throw std::invalid_argument("Sigma must be Px(PN) see documentation");
    if (N == 0 ) throw std::invalid_argument("Sigma must have columns");
    N = N/P1; // safe after above validation
    MatrixXd S;
    MatrixXd res(P1, N*P1);
    for (int n=0; n<N; n++){
      S = Sigma.middleCols(n*P1, P1);
      res.middleCols(n*P1,P1) = internal::clrvar2vararray_single(S);
    }
    return res;
  } 
  
  
  //' Calculate Phi Statistics (proportionality) from CLR Covariances
  //' @param Sigma Covariance matrix Px(PN) where N is number of 
  //'   covariance matricies in CLR space
  //' See Lovell "Proportionality: A valid alternative ..." for details. 
  template <typename TS>
  Eigen::MatrixXd clrvar2phi(Eigen::MatrixBase<TS>& Sigma){
    int P1 = Sigma.rows();
    int N = Sigma.cols();
    if ( (N % P1) != 0 ) throw std::invalid_argument("Sigma must be Px(PN) see documentation");
    if (N == 0 ) throw std::invalid_argument("Sigma must have columns");
    N = N/P1; // safe after above validation
    MatrixXd S;
    MatrixXd res(P1, N*P1);
    for (int n=0; n<N; n++){
      S = Sigma.middleCols(n*P1, P1);
      res.middleCols(n*P1,P1) = internal::clrvar2phi_single(S);
    }
    return res;
  } 

}

#endif 
