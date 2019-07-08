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
    MatrixXd Y.noalias() = V.transpose()*X;
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
    return V1.transpose();
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
    return transferContrast(V1, V2);
  }
  
  
  //' CLR to ALR (don't bother if you know what your doing)
  //' @param d2 contrast matrix of basis in currently
  //' @param D total number of parts
  Eigen::MatrixXd caContrast(int d2, int D){
    return alrContrast(d2, D, false);
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
    return transferContrast(V1, V2);
  }
  
  //' ALR to ILR
  //' @param d1 alr reference part is in currently
  //' @param V2 contrast matrix of basis to go into
  //' @param D number of total parts
  template <typename T>
  Eigen::MatrixXd aiContrast(int d1, Eigen::MatrixBase<T>& V2, int D){
    MatrixXd V1 = acContrast(d1, D);
    return transferContrast(V1, V2);
  }
  
  //' transfer data between coordinate systems
  //' @param X data
  //' @param V transferContrast matrix
  template <typename TX, typename TV>
  Eigen::MatrixXd transferData(Eigen::MatrixBase<TX>& X, 
                           Eigen::MatrixBase<TV>& V){
    return V*X;
  }
  
  template <typename TX, typename TV1, typename TV2>
  Eigen::MatrixXd ilr2ilr(Eigen::MatrixBase<TX>& X, 
                          Eigen::MatrixBase<TV1>& V1, 
                          Eigen::MatrixBase<TV2>& V2){
    MatrixXd V = iiContrast(V1, V2);
    return transferData(X, V);
  }
  
  template <typename TX, typename TV1>
  Eigen::MatrixXd ilr2clr(Eigen::MatrixBase<TX>& X, 
                          Eigen::MatrixBase<TV1>& V1){
    MatrixXd V = icContrast(V1);
    return transferData(X, V);
  }
  
  template <typename TX, typename TV2>
  Eigen::MatrixXd clr2ilr(Eigen::MatrixBase<TX>& X, 
                          Eigen::MatrixBase<TV2>& V2){
    MatrixXd V = ciContrast(V2);
    return transferData(X, V);
  }
  
  template <typename TX>
  Eigen::MatrixXd alr2clr(Eigen::MatrixBase<TX>& X, int d1){
    MatrixXd V = acContrast(d1, X.rows()+1);
    return transferData(X, V);
  }
  
  template <typename TX>
  Eigen::MatrixXd clr2alr(Eigen::MatrixBase<TX>& X, int d2){
    MatrixXd V = caContrast(d2, X.rows());
    return transferData(X, V);
  }
  
  template <typename TX>
  Eigen::MatrixXd alr2alr(Eigen::MatrixBase<TX>& X, int d1, int d2){
    MatrixXd V = aaContrast(d1, d2, X.rows()+1);
    return transferData(X, V);
  }
  
  template <typename TX, typename TV>
  Eigen::MatrixXd ilr2alr(Eigen::MatrixBase<TX>& X, 
                          Eigen::MatrixBase<TV>& V1, int d2){
    MatrixXd V = iaContrast(V1, d2, X.rows()+1);
    return transferData(X, V);
  }
  
  template <typename TX, typename TV>
  Eigen::MatrixXd alr2ilr(Eigen::MatrixBase<TX>& X,int d1, Eigen::MatrixBase<TV>& V2){
    MatrixXd V = aiContrast(d1, V2, X.rows()+1);
    return transferData(X, V);
  }
    
    
  //' transfer Covariance between coordinate systems
  //' @param Sigma covariance matrix - spd
  //' @param V transferContrast matrix
  template <typename TS, typename TV>
  Eigen::MatrixXd transferCovariance(Eigen::MatrixBase<TS>& Sigma, 
                                     Eigen::MatrixBase<TV>& V){
    return V*Sigma.template selfadjointView<Eigen::Lower>()*V.transpose();
  }

  //' common internal to transfer covariance functions
  //' @param Sigma Covariance matrix Px(PN) where N is number of covariance matricies
  //' @param V contrast matrix to use for transfer
  template <typename TS, typename TV>
  Eigen::MatrixXd transferCovarianceIterate(Eigen::MatrixBase<TS>& Sigma, 
                                              Eigen::MatrixBase<TV>& V){
    int P1 = V.cols();
    int P2 = V.rows();
    int N = Sigma.cols();
    if ( N % P1 != 0 ) throw std::invalid_argument("Sigma must be Px(PN) see documentation");
    if (N == 0 ) throw std::invalid_argument("Sigma must have columns");
    N = N/P1; // safe after above validation
    MatrixXd res(P2, N*P2);
    for (int i=0; i<N; i++){
      Map<MatrixXd> S(&Sigma(i*P1), P2, P2);
      res.middleRows(i*P2,(i+1)*P2-1)= transferCovariance(S, V);
    }
    return res;
  }
  
  
  template <typename TS, typename TV1, typename TV2>
  Eigen::MatrixXd ilrvar2ilrvar(Eigen::MatrixBase<TS>& X, 
                          Eigen::MatrixBase<TV1>& V1, 
                          Eigen::MatrixBase<TV2>& V2){
    MatrixXd V = iiContrast(V1, V2);
    return transferCovarianceIterate(X, V);
  }
  
  template <typename TX, typename TV1>
  Eigen::MatrixXd ilrvar2clrvar(Eigen::MatrixBase<TX>& X, 
                          Eigen::MatrixBase<TV1>& V1){
    MatrixXd V = icContrast(V1);
    return transferCovarianceIterate(X, V);
  }
  
  template <typename TX, typename TV2>
  Eigen::MatrixXd clrvar2ilrvar(Eigen::MatrixBase<TX>& X, 
                          Eigen::MatrixBase<TV2>& V2){
    MatrixXd V = ciContrast(V2);
    return transferCovarianceIterate(X, V);
  }
  
  template <typename TX>
  Eigen::MatrixXd alrvar2clrvar(Eigen::MatrixBase<TX>& X, int d1){
    MatrixXd V = acContrast(d1, X.rows()+1);
    return transferCovarianceIterate(X, V);
  }
  
  template <typename TX>
  Eigen::MatrixXd clrvar2alrvar(Eigen::MatrixBase<TX>& X, int d2){
    MatrixXd V = caContrast(d2, X.rows());
    return transferCovarianceIterate(X, V);
  }
  
  template <typename TX>
  Eigen::MatrixXd alrvar2alrvar(Eigen::MatrixBase<TX>& X, int d1, int d2){
    MatrixXd V = aaContrast(d1, d2, X.rows()+1);
    return transferCovarianceIterate(X, V);
  }
  
  template <typename TX, typename TV>
  Eigen::MatrixXd ilrvar2alrvar(Eigen::MatrixBase<TX>& X, 
                          Eigen::MatrixBase<TV>& V1, int d2){
    MatrixXd V = iaContrast(V1, d2, X.rows()+1);
    return transferCovarianceIterate(X, V);
  }
  
  template <typename TX, typename TV>
  Eigen::MatrixXd alrvar2ilrvar(Eigen::MatrixBase<TX>& X,int d1, Eigen::MatrixBase<TV>& V2){
    MatrixXd V = aiContrast(d1, V2, X.rows()+1);
    return transferCovarianceIterate(X, V);
  }
  
}



#endif
