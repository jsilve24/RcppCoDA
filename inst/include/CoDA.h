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
  
  //' Create alr contrast matrix
  //' @param D total number of parts
  //' @param d reference part for alr
  //' @param inv (is this being used for the alr or alrInv?)
  //' @return D-1 x D matrix
  Eigen::SparseMatrix<double> alrContrastSparse(int d, int D, bool inv){
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
  Eigen::MatrixXd alrContrast(int d, int D, bool inv){
    Eigen::MatrixXd B;
    SparseMatrix<double> Bsp = alrContrastSparse(d, D, inv);
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
  //' @param inv (is this being used for the clr or clrInv?)
  //' @return D x D contrast matrix 
  Eigen::MatrixXd clrContrast(int D, bool inv){
    if (D < 2) throw std::invalid_argument("D must be >1");
    if (!inv){
      MatrixXd B = -MatrixXd::Ones(D,D);
      B.diagonal().array() += D;
      return B.array()/D;  
    } else {
      MatrixXd B = MatrixXd::Identity(D,D);
      return B;
    }
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
    MatrixXd Y;
    Y.noalias() = V.transpose()*X;
    Y = Y.array().exp().matrix();
    return clo(Y);
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
  //' @export
  template <typename TX, typename TV>
  Eigen::MatrixXd glrInv(Eigen::MatrixBase<TX>& X, Eigen::SparseMatrixBase<TV>& V){
    MatrixXd Y = V.transpose()*X;
    Y = Y.array().exp().matrix();
    return clo(Y);
  }

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
    //MatrixXd V = clrContrast(X.cols(), false);
    //return glr(X, V);
  }
  
  //' Inverse Centered LR Transform
  //' @param X data (coords x samples aka P x N) here P = D because its the CLR
  template <typename T>
  Eigen::MatrixXd clrInv(Eigen::MatrixBase<T>& X){
    MatrixXd Y = X.array().exp().matrix();
    return clo(Y);
    //MatrixXd V = clrContrast(X.cols(), true);
    //return glrInv(X, V);
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
  
  
  //' @param V1 contrast matrix of basis in currently
  //' @param V2 contrast matrix of basis to go to
  //' @return V2*V1.transpose()
  template <typename TV1, typename TV2>
  Eigen::MatrixXd transferContrast(Eigen::MatrixBase<TV1>& V1, 
                                   Eigen::MatrixBase<TV2>& V2){
    return V2*V1.transpose();
  }
  
  //' @param V1 contrast matrix of basis in currently
  //' @param V2 contrast matrix of basis to go to
  //' @return V2*V1.transpose()
  template <typename TV1, typename TV2>
  Eigen::MatrixXd transferContrast(Eigen::MatrixBase<TV1>& V1, 
                                   
                                   Eigen::SparseMatrixBase<TV2>& V2){
    return V2*V1.transpose();
  }
  
  //' @param V1 contrast matrix of basis in currently
  //' @param V2 contrast matrix of basis to go to
  //' @return V2*V1.transpose()
  template <typename TV1, typename TV2>
  Eigen::MatrixXd transferContrast(Eigen::SparseMatrixBase<TV1>& V1, 
                                   Eigen::MatrixBase<TV2>& V2){
    return V2*V1.transpose();
  }
  
  //' @param V1 contrast matrix of basis in currently
  //' @param V2 contrast matrix of basis to go to
  //' @return V2*V1.transpose()
  template <typename TV1, typename TV2>
  Eigen::MatrixXd transferContrast(Eigen::SparseMatrixBase<TV1>& V1, 
                                   Eigen::SparseMatrixBase<TV2>& V2){
    return V2*V1.transpose();
  }
  
  //' ILR to ILR (just call transferContrast if you know what your doing)
  //' @param V1 contrast matrix of basis in currently
  //' @param V2 contrast matrix of basis to go to
  template <typename TV1, typename TV2>
  Eigen::MatrixXd iiContrast(Eigen::MatrixBase<TV1>& V1, 
                             Eigen::MatrixBase<TV2>& V2){
    return transferContrast(V1, V2);
  }
  
  //' ILR to CLR (don't bother if you know what your doing)
  //' @param V1 contrast matrix of basis in currently
  template <typename TV1>
  Eigen::MatrixXd icContrast(Eigen::MatrixBase<TV1>& V1){
    MatrixXd V = V1.transpose();
    return V;
  }
  
  //' CLR to ILR (don't bother if you know what your doing)
  //' @param V2 contrast matrix of basis in currently
  template <typename TV2>
  Eigen::MatrixXd ciContrast(Eigen::MatrixBase<TV2>& V2){
    return V2;
  }
  
  
  //' ILR to ALR 
  //' @param V1 contrast matrix of basis in currently
  //' @param d2 alr reference part to go into
  //' @param D number of total parts
  template <typename T>
  Eigen::MatrixXd iaContrast(Eigen::MatrixBase<T>& V1, int d2, int D){
    SparseMatrix<double> V2 = alrContrastSparse(d2, D, false);
    MatrixXd V = transferContrast(V1, V2);
    return V;
  }
  
  
  //' CLR to ALR (don't bother if you know what your doing)
  //' @param d2 contrast matrix of basis in currently
  //' @param D total number of parts
  Eigen::MatrixXd caContrast(int d2, int D){
    MatrixXd V = alrContrast(d2, D, false);
    return V;
  }
  
  //' ALR to CLR
  //' @param d1 contrast matrix of basis in currently
  //' @param D total number of parts
  Eigen::MatrixXd acContrast(int d1, int D){
    MatrixXd G = alrContrast(d1, D, true);
    G.array() -= 1/D;
    return G;
  }
  
  //' ALR to ALR
  //' @param d1 alr reference part is in currently
  //' @param d2 alr reference part to go into
  //' @param D number of total parts
  Eigen::MatrixXd aaContrast(int d1, int d2, int D){
    MatrixXd V1 = acContrast(d1, D);
    SparseMatrix<double> V2 = alrContrastSparse(d2, D, false);
    MatrixXd V = transferContrast(V1, V2);
    return V;
  }
  
  //' ALR to ILR
  //' @param d1 alr reference part is in currently
  //' @param V2 contrast matrix of basis to go into
  //' @param D number of total parts
  template <typename T>
  Eigen::MatrixXd aiContrast(int d1, Eigen::MatrixBase<T>& V2, int D){
    MatrixXd V1 = acContrast(d1, D);
    MatrixXd V = transferContrast(V1, V2);
    return V;
  }
  
  //' transfer data between coordinate systems
  //' @param X data
  //' @param V transferContrast matrix
  template <typename TX, typename TV>
  Eigen::MatrixXd transferData(Eigen::MatrixBase<TX>& X, 
                           Eigen::MatrixBase<TV>& V){
    MatrixXd O = V*X;
    return O;
  }
  
  template <typename TX, typename TV1, typename TV2>
  Eigen::MatrixXd ilr2ilr(Eigen::MatrixBase<TX>& X, 
                          Eigen::MatrixBase<TV1>& V1, 
                          Eigen::MatrixBase<TV2>& V2){
    MatrixXd V = iiContrast(V1, V2);
    MatrixXd Y = transferData(X, V);
    return Y;
  }
  
  template <typename TX, typename TV1>
  Eigen::MatrixXd ilr2clr(Eigen::MatrixBase<TX>& X, 
                          Eigen::MatrixBase<TV1>& V1){
    MatrixXd V = icContrast(V1);
    MatrixXd Y = transferData(X, V);
    return Y;
  }
  
  template <typename TX, typename TV2>
  Eigen::MatrixXd clr2ilr(Eigen::MatrixBase<TX>& X, 
                        Eigen::MatrixBase<TV2>& V2){
    MatrixXd V = ciContrast(V2);
    MatrixXd Y = transferData(X, V);
    return Y;
  }
  
  template <typename TX>
  Eigen::MatrixXd alr2clr(Eigen::MatrixBase<TX>& X, int d1){
    MatrixXd V = acContrast(d1, X.rows()+1);
    MatrixXd Y = transferData(X, V);
    return Y;
  }
  
  template <typename TX>
  Eigen::MatrixXd clr2alr(Eigen::MatrixBase<TX>& X, int d2){
    MatrixXd V = caContrast(d2, X.rows());
    MatrixXd Y = transferData(X, V);
    return Y;
  }
  
  template <typename TX>
  Eigen::MatrixXd alr2alr(Eigen::MatrixBase<TX>& X, int d1, int d2){
    MatrixXd V = aaContrast(d1, d2, X.rows()+1);
    MatrixXd Y = transferData(X, V);
    return Y;
  }
  
  template <typename TX, typename TV>
  Eigen::MatrixXd ilr2alr(Eigen::MatrixBase<TX>& X, 
                          Eigen::MatrixBase<TV>& V1, int d2){
    MatrixXd V = iaContrast(V1, d2, X.rows()+1);
    MatrixXd Y = transferData(X, V);
    return Y;
  }
  
  template <typename TX, typename TV>
  Eigen::MatrixXd alr2ilr(Eigen::MatrixBase<TX>& X,int d1, Eigen::MatrixBase<TV>& V2){
    MatrixXd V = aiContrast(d1, V2, X.rows()+1);
    MatrixXd Y = transferData(X, V);
    return Y;
  }
    
    
  //' transfer Covariance between coordinate systems
  //' @param Sigma covariance matrix - spd
  //' @param V transferContrast matrix
  template <typename TS, typename TV>
  Eigen::MatrixXd transferCovariance(Eigen::MatrixBase<TS>& Sigma, 
                                     Eigen::MatrixBase<TV>& V){
    MatrixXd O =  V*Sigma*V.transpose();
    return O;
  }

  //' common internal to transfer covariance matricies
  //' @param Sigma Covariance matrix Px(PN) where N is number of covariance matricies
  //' @param V contrast matrix to use for transfer
  template <typename TS, typename TV>
  Eigen::MatrixXd transferCovarianceIterate(Eigen::MatrixBase<TS>& Sigma, 
                                              Eigen::MatrixBase<TV>& V){
    int P1 = V.cols();
    int P2 = V.rows();
    int N = Sigma.cols();
    if ( (N % P1) != 0 ) throw std::invalid_argument("Sigma must be Px(PN) see documentation");
    if (N == 0 ) throw std::invalid_argument("Sigma must have columns");
    N = N/P1; // safe after above validation
    MatrixXd res(P2, N*P2);
    MatrixXd S;
    for (int i=0; i<N; i++){
      S = Sigma.middleCols(i*P1, P1); //(&Sigma(i*P1), P1, P1);
      res.middleCols(i*P2, P2)= transferCovariance(S, V);
    }
    return res;
  }
  
  
  template <typename TS, typename TV1, typename TV2>
  Eigen::MatrixXd ilrvar2ilrvar(Eigen::MatrixBase<TS>& X, 
                          Eigen::MatrixBase<TV1>& V1, 
                          Eigen::MatrixBase<TV2>& V2){
    MatrixXd V = iiContrast(V1, V2);
    MatrixXd O = transferCovarianceIterate(X, V);
    return O;
  }
  
  template <typename TX, typename TV1>
  Eigen::MatrixXd ilrvar2clrvar(Eigen::MatrixBase<TX>& X, 
                          Eigen::MatrixBase<TV1>& V1){
    MatrixXd V = icContrast(V1);
    MatrixXd O = transferCovarianceIterate(X, V);
    return O;
  }
  
  template <typename TX, typename TV2>
  Eigen::MatrixXd clrvar2ilrvar(Eigen::MatrixBase<TX>& X, 
                          Eigen::MatrixBase<TV2>& V2){
    MatrixXd V = ciContrast(V2);
    MatrixXd O = transferCovarianceIterate(X, V);
    return O;
  }
  
  template <typename TX>
  Eigen::MatrixXd alrvar2clrvar(Eigen::MatrixBase<TX>& X, int d1){
    MatrixXd V = acContrast(d1, X.rows()+1);
    MatrixXd O = transferCovarianceIterate(X, V);
    return O;
  }
  
  template <typename TX>
  Eigen::MatrixXd clrvar2alrvar(Eigen::MatrixBase<TX>& X, int d2){
    MatrixXd V = caContrast(d2, X.rows());
    MatrixXd O = transferCovarianceIterate(X, V);
    return O;
  }
  
  template <typename TX>
  Eigen::MatrixXd alrvar2alrvar(Eigen::MatrixBase<TX>& X, int d1, int d2){
    MatrixXd V = aaContrast(d1, d2, X.rows()+1);
    MatrixXd O = transferCovarianceIterate(X, V);
    return O;
  }
  
  template <typename TX, typename TV>
  Eigen::MatrixXd ilrvar2alrvar(Eigen::MatrixBase<TX>& X, 
                          Eigen::MatrixBase<TV>& V1, int d2){
    MatrixXd V = iaContrast(V1, d2, X.rows()+1);
    MatrixXd O = transferCovarianceIterate(X, V);
    return O;
  }
  
  template <typename TX, typename TV>
  Eigen::MatrixXd alrvar2ilrvar(Eigen::MatrixBase<TX>& X,int d1, Eigen::MatrixBase<TV>& V2){
    MatrixXd V = aiContrast(d1, V2, X.rows()+1);
    MatrixXd O = transferCovarianceIterate(X, V);
    return O;
  }
  
  //' Calculate Phi statistics (propotionality) from CLR covariance (single matrix)
  //' @param Sigma Covariance matrix PxP in CLR space
  template <typename TS>
  Eigen::MatrixXd clrvar2phi_single(Eigen::MatrixBase<TS>& Sigma){
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
  
  //' Calculate Phi statistics (proportionality) from CLR Covariances
  //' @param Sigma Covariance matrix Px(PN) where N is number of 
  //'   covariance matricies in CLR space
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
        res.middleCols(n*P1,P1) = clrvar2phi_single(S);
      }
      return res;
  } 
  
  
  
  
  
  
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
      res.middleCols(n*P1,P1) = transferCovariance(S, V);
    }
    return res;
  }
  
  
  // BLOCK OPERATIONS ------------
  
  //' Generalized LR Transform
  //' @param X data (parts x samples aka D x N; D = D1+D2)
  //' @param V1 contrast matrix (P1 x D1)
  //' @param V2 contrast matrix (P2 x D2)
  //' @details assumes that if only V1 is given then its the first block to be transformed
  template <typename TX, typename TV1>
  Eigen::MatrixXd glrBlock(Eigen::MatrixBase<TX>& X, Eigen::MatrixBase<TV1>& V1){
    int P1=V1.rows();
    int D1=V1.cols();
    int D = X.rows();
    int N = X.cols();
    if (X.rows() < D1) throw std::invalid_argument("X.rows() >= V1.cols()");
    MatrixXd O(D-D1+P1, N);
    MatrixXd Y = X.topRows(D1).array().log().matrix();
    O.topRows(P1).noalias() = V1*Y;
    O.bottomRows(D-D1) = X.bottomRows(D-D1);
    return O;
  }
  
  template <typename TX, typename TV1, typename TV2>
  Eigen::MatrixXd glrBlock(Eigen::MatrixBase<TX>& X, Eigen::MatrixBase<TV1>& V1, 
                           Eigen::MatrixBase<TV2>& V2){
    int P1=V1.rows();
    int D1=V1.cols();
    int P2 = V2.rows();
    int D2 = V2.cols();
    int N = X.cols();
    if (X.rows() != D1+D2) throw std::invalid_argument("X.rows() must equal V1.cols() + V2.cols()");
    MatrixXd O(P1+P2, N);
    MatrixXd Y = X.array().log().matrix();
    O.topRows(P1).noalias() = V1*Y.bottomRows(D1);
    O.bottomRows(P2) = V2*Y.bottomRows(D2);
    return O;
  }
  
  template <typename TX, typename TV1>
  Eigen::MatrixXd glrInvBlock(Eigen::MatrixBase<TX>& X, Eigen::MatrixBase<TV1>& V1){
    int P1=V1.rows();
    int D1=V1.cols();
    int D = X.rows();
    int N = X.cols();
    if (X.rows() < P1) throw std::invalid_argument("X.rows() >= V1.rows()");
    MatrixXd Y;
    Y.noalias() = V1.transpose()*X.topRows(P1);
    Y = Y.array().exp().matrix();
    MatrixXd O(D-P1+D1, N);
    O.topRows(D1) = clo(Y);
    O.bottomRows(D-P1) = X.bottomRows(D-P1);
    return O;
  }
  
  template <typename TX, typename TV1, typename TV2>
  Eigen::MatrixXd glrInvBlock(Eigen::MatrixBase<TX>& X, Eigen::MatrixBase<TV1>& V1, 
                              Eigen::MatrixBase<TV2>& V2){
    int P1=V1.rows();
    int D1=V1.cols();
    int P2 = V2.rows();
    int D2 = V2.cols();
    int N = X.cols();
    if (X.rows() != P1 + P2) throw std::invalid_argument("X.rows() >= V1.rows() + V2.rows()");
    MatrixXd Y(D1+D2, N);
    Y.topRows(D1).noalias() = V1.transpose()*X.topRows(P1);
    Y.bottomRows(D2).noalias() = V2.transpose()*X.bottomRows(P2);
    Y = Y.array().exp().matrix();
    MatrixXd O(D1+D2, N);
    MatrixXd Ytmp = Y.topRows(D1);
    O.topRows(D1) = clo(Ytmp);
    Ytmp = Y.bottomRows(D2);
    O.bottomRows(D2) = clo(Ytmp);
    return O;
  }
  
  template <typename TX, typename TV1>
  Eigen::MatrixXd transferDataBlock(Eigen::MatrixBase<TX>& X, 
                               Eigen::MatrixBase<TV1>& V1){
    int P1=V1.rows();
    int D1=V1.cols();
    int D = X.rows();
    int N = X.cols();
    if (X.rows() < D1) throw std::invalid_argument("X.rows() >= V1.cols()");
    MatrixXd O(D-D1+P1, N);
    O.topRows(P1) = V1*X.topRows(D1);
    O.bottomRows(D-D1) = X.bottomRows(D-D1);  
    return O;
  }
  
  template <typename TX, typename TV1, typename TV2>
  Eigen::MatrixXd transferDataBlock(Eigen::MatrixBase<TX>& X,
                                    Eigen::MatrixBase<TV1>& V1,
                                    Eigen::MatrixBase<TV2>& V2){
    int P1=V1.rows();
    int D1=V1.cols();
    int P2 = V2.rows();
    int D2 = V2.cols();
    int N = X.cols();
    if (X.rows() != D1 + D2) throw std::invalid_argument("X.rows() >= V1.cols() + V2.cols()");
    MatrixXd O(P2+P1, N);
    O.topRows(D1) = V1*X.topRows(P1);
    O.bottomRows(D2) = V2*X.bottomRows(P2);
    return O;
  }
  
  template <typename TS, typename TV1>
  Eigen::MatrixXd transferCovarianceBlock(Eigen::MatrixBase<TS>& Sigma, 
                                     Eigen::MatrixBase<TV1>& V1){
    int P1=V1.rows();
    int D1=V1.cols();
    int D = Sigma.rows();
    if (D != Sigma.cols()) throw std::invalid_argument("Sigma must be a square matrix");
    if (D < D1) throw std::invalid_argument("Sigma dimension must be larger than V1.cols()");
    MatrixXd O(D-D1+P1, D-D1+P1);
    O.topLeftCorner(P1, P1).noalias() =  V1*Sigma.topLeftCorner(D1,D1)*V1.transpose();
    O.topRightCorner(P1, D-D1).noalias() = V1*Sigma.topRightCorner(D1, D-D1);
    O.bottomLeftCorner(D-D1, P1) = O.topRightCorner(P1, D-D1).transpose();
    O.bottomRightCorner(D-D1, D-D1) = Sigma.bottomRightCorner(D-D1, D-D1);
    return O;
  }
  
  template <typename TS, typename TV1, typename TV2>
  Eigen::MatrixXd transferCovarianceBlock(Eigen::MatrixBase<TS>& Sigma, 
                                          Eigen::MatrixBase<TV1>& V1, 
                                          Eigen::MatrixBase<TV2>& V2){
    int P1=V1.rows();
    int D1=V1.cols();
    int P2=V2.rows();
    int D2=V2.cols();
    int D = Sigma.rows();
    if (D != Sigma.cols()) throw std::invalid_argument("Sigma must be a square matrix");
    if (D != D1+D2) throw std::invalid_argument("Sigma dimension must = V1.cols() + V2.cols()");
    MatrixXd O(P1+P2, P1+P2);
    O.topLeftCorner(P1, P1).noalias() =  V1*Sigma.topLeftCorner(D1,D1)*V1.transpose();
    O.topRightCorner(P1, P2).noalias() = V1*Sigma.topRightCorner(D1, D2)*V2.transpose();
    O.bottomLeftCorner(P2, P1) = O.topRightCorner(P1, P2).transpose();
    O.bottomRightCorner(P2, P2).noalias() = V2*Sigma.bottomRightCorner(D2, D2)*V2.transpose();
    return O;
  }
  
  template <typename TS, typename TV1>
  Eigen::MatrixXd transferCovarianceIterateBlock(Eigen::MatrixBase<TS>& Sigma, 
                                            Eigen::MatrixBase<TV1>& V1){
    int P1=V1.rows();
    int D1=V1.cols();
    int D = Sigma.rows();
    int D2 = D-D1;
    int N = Sigma.cols();
    if (D < D1) throw std::invalid_argument("Sigma dimension must be larger than V1.cols()");
    if ( (N % D ) != 0 ) throw std::invalid_argument("Sigma must be Px(PN) see documentation");
    if (N == 0 ) throw std::invalid_argument("Sigma must have columns");
    N = N/D; // safe after above validation
    MatrixXd O(P1+D2, N*(P1+D2));
    MatrixXd S;
    for (int i=0; i<N; i++){
      S = Sigma.middleCols(i*D, D); 
      O.middleCols(i*(P1+D2), P1+D2)= transferCovarianceBlock(S, V1);
    }
    return O;
  }
  
  template <typename TS, typename TV1, typename TV2>
  Eigen::MatrixXd transferCovarianceIterateBlock(Eigen::MatrixBase<TS>& Sigma, 
                                                 Eigen::MatrixBase<TV1>& V1, 
                                                 Eigen::MatrixBase<TV2>& V2){
    int P1=V1.rows();
    int D1=V1.cols();
    int P2 = V2.rows();
    int D2 = V2.cols();
    int D = Sigma.rows();
    int N = Sigma.cols();
    if (D != D1+D2) throw std::invalid_argument("Sigma dimension must = V1.cols() + V2.cols()");
    if ( (N % D ) != 0 ) throw std::invalid_argument("Sigma must be Px(PN) see documentation");
    if (N == 0 ) throw std::invalid_argument("Sigma must have columns");
    N = N/D; // safe after above validation
    MatrixXd O(P1+P2, N*(P1+P2));
    MatrixXd S;
    for (int i=0; i<N; i++){
      S = Sigma.middleCols(i*D, D); 
      O.middleCols(i*(P1+P2), P1+P2)= transferCovarianceBlock(S, V1, V2);
    }
    return O;
  }
  
}

#endif
