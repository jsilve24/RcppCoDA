#ifndef RCPPCODA__IQLR__H
#define RCPPCODA__IQLR__H


#include <RcppEigen.h>
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

  //' Calculate Contrast matrix for IQLR from D-dimensional vector of variances
  //' This is an advanced feature and should not be used without prior knowledge of the
  //' IQLR transform. 
  //' @param S D-vector of "variances" for IQLR (e.g., diagonal of CLR covariance matrix)
  //' @param qLow lower quantile (between 0 and 1)
  //' @param qHigh upper quantile (between 0 and 1; must be greater than qLow)
  //' @return DxD contrast matrix 
  template <typename TS>
  Eigen::MatrixXd var2iqlrContrast(Eigen::MatrixBase<TS>& S, double qLow, double qHigh){
    int D = S.rows();
    if  (S.cols() > 1) throw std::invalid_argument("S must be a vector not a matrix");
    if (qLow < 0 | qLow > 1) std::invalid_argument("qLow must be between 0 and 1");
    if (qHigh < 0 | qHigh > 1) std::invalid_argument("qLow must be between 0 and 1");
    if (qHigh <= qLow) throw std::invalid_argument("qHigh must be > qLow");
    
    // Calculate Quantiles
    std::vector<double> sv;
    for (int i=0; i<D; i++){
      sv.push_back(S(i,0));
    }
    std::sort(sv.begin(), sv.end());
    int iLow =  std::floor(qLow*(D-1));
    int iHigh =  std::ceil(qHigh*(D-1));
    double Low = sv[iLow];
    double High = sv[iHigh];
    double nIQVF = iHigh-iLow+1.0;
    MatrixXd V = MatrixXd::Zero(D,D);
    VectorXd c = MatrixXd::Constant(D, 1, 1/nIQVF);
    for (int i=0; i<D; i++){
      if ( (S(i,0) >= Low) & (S(i,0) <= High) ) {
        V.col(i) = c;
      }
    }
    V = MatrixXd::Identity(D,D) - V;
    return V;
  }
  
  //' Transfer CLR Covariance Matricies into IQLR 
  //' 
  //' Uses Diagonal of covariance matricies for calculation of IQLR
  //' @param Sigma (PN)xP covariance matrix in CLR
  //' @param qLow lower quantile (between 0 and 1)
  //' @param qHigh upper quantile (between 0 and 1; must be greater than qLow)
  //' @return (PN)xP Covariance matrix in IQLR coordinates 
  template <typename TS>
  Eigen::MatrixXd clrvar2iqlrvar(Eigen::MatrixBase<TS>& Sigma, double qLow, double qHigh){
    int P1 = Sigma.rows();
    int N = Sigma.cols();
    if ( (N % P1) != 0 ) throw std::invalid_argument("Sigma must be Px(PN) see documentation");
    if (N == 0 ) throw std::invalid_argument("Sigma must have columns");
    N = N/P1; // safe after above validation
    MatrixXd S, V;
    VectorXd s;
    MatrixXd res(P1, N*P1);
    for (int n=0; n<N; n++){
      S = Sigma.middleCols(n*P1, P1);
      s = S.diagonal();
      V = var2iqlrContrast(s, qLow, qHigh);
      res.middleCols(n*P1,P1) = quadForm(S, V);
    }
    return res;
  }

}

#endif 