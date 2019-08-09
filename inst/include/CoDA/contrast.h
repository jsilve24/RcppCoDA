#ifndef CODA__CONTRAST__H
#define CODA__CONTRAST__H

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

  //' Create alr contrast matrix
  //' @param D total number of parts
  //' @param d reference part for alr
  //' @param inv (is this being used for the alr or alrInv?)
  //' @return D-1 x D matrix
  typedef Eigen::Triplet<double> T;
  Eigen::SparseMatrix<double> alrContrastSparse(int d, int D, bool inv){
    if ((d<1) || (d>D)) throw std::invalid_argument("d must be within [1,...,D]");
    std::vector<T> tripletList;
    if (inv) tripletList.reserve(D);
    if (!inv) tripletList.reserve(2*D);
    int pos=0;
    for (int i=0; i<D; i++){
      if (i == (d-1)) pos++;
      if (pos >= D) break;
      tripletList.push_back(T(i,pos, 1));
      //B.insert(i, pos) = 1;
      pos++;
    }
    if (!inv) {
      for (int i=0; i<(D-1); i++){
        tripletList.push_back(T(i,d-1, -1));
        //B.insert(i,d-1) = -1;
      }
    }
    SparseMatrix<double> B(D-1, D);
    B.setFromTriplets(tripletList.begin(), tripletList.end());
    return B;
  }
  // Eigen::SparseMatrix<double> alrContrastSparse(int d, int D, bool inv){
  //   if ((d<1) || (d>D)) throw std::invalid_argument("d must be within [1,...,D]");
  //   SparseMatrix<double> B(D-1, D);// = MatrixXd::Zero(D-1, D);
  //   B.reserve(Eigen::VectorXi::Constant(D, 2)); // reserve 2 non-zero per column
  //   int pos=0;
  //   for (int i=0; i<D; i++){
  //     if (i == (d-1)) pos++;
  //     B.insert(i, pos) = 1;
  //     pos++;
  //   }
  //   if (!inv) {
  //     for (int i=0; i<(D-1); i++){
  //       B.insert(i,d-1) = -1;
  //     }
  //   }
  //   return B;
  // }
  
  //' Create alr contrast matrix
  //' @param D total number of parts
  //' @param d reference part for alr
  //' @param inv (is this being used for the alr or alrInv?)
  //' @return D-1 x D matrix
  Eigen::MatrixXd alrContrast(int d, int D, bool inv){
    Eigen::MatrixXd B;
    SparseMatrix<double> Bsp = alrContrastSparse(d, D, inv);
    B = MatrixXd(Bsp);
    return B;
  }
  
  
  // //' Create ilr default contrast matrix
  // //' @param D total number of parts
  // //' @return D-1xD contrast matrix (a basis)
  // Eigen::MatrixXd ilrContrast(int D){
  //   if (D < 2) throw std::invalid_argument("D must be > 1");
  //   MatrixXd B = alrContrast(D, D, false);
  //   MatrixXd Q = B.transpose().householderQr().householderQ().transpose();
  //   return Q.topRows(D-1);
  // }
  //' Create ilr default contrast matrix
  //' @param D total number of parts
  //' @return D-1xD contrast matrix (a basis)
  Eigen::MatrixXd ilrContrast(int D){
    MatrixXd V(D-1, D);
    double ap, an;
    for (int i=1; i<=(D-1); i++){
      ap = sqrt(1.0/( (D-i)*(D-i+1.0)  ) );
      an = -sqrt((D-i)/ (D-i+1.0) );
      for (int j=1; j<=D; j++){
        if (j <= D-i) V(i-1,j-1) = ap;
        else if (j == D-i+1) V(i-1,j-1) = an;
        else V(i-1,j-1)=0;
      }
    }
    return V;
  }
  
  
  //' Create ilr contrast from sign matrix
  //' @param S D-1 x D sign matrix (all entries either 0, 1, or -1)
  //' @return D-1 x D contrast matrix (a basis)
  template <typename TS>
  Eigen::MatrixXd ilrContrast(Eigen::MatrixBase<TS>& S){
    int D = S.cols();
    MatrixXd V(D-1, D);
    double np, nn; // number positive / negative
    double ap, an, a;
    for (int i=0; i<(D-1); i++){
      np=0;
      nn=0;
      for (int j=0; j<D; j++){
        if (S(i,j)==1) np += 1.0;
        else if (S(i,j)==-1) nn += 1.0;
        else if (S(i,j) != 0) throw std::invalid_argument("S must be a valid sign matrix");
      }
      a = sqrt(nn*np/(nn+np));
      ap = a/np;
      an = -a/nn;
      for (int j=0; j<D; j++){
        if (S(i,j)==1) V(i,j) = ap;
        else if (S(i,j)==-1) V(i,j) = an;
        else V(i,j)=0;
      }
    }
    return V;
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
  
  namespace internal {

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
  
  }
  
  //' ILR to ILR (just call transferContrast if you know what your doing)
  //' @param V1 contrast matrix of basis in currently
  //' @param V2 contrast matrix of basis to go to
  template <typename TV1, typename TV2>
  Eigen::MatrixXd iiTransfer(Eigen::MatrixBase<TV1>& V1, 
                             Eigen::MatrixBase<TV2>& V2){
    return internal::transferContrast(V1, V2);
  }
  
  //' ILR to CLR (don't bother if you know what your doing)
  //' @param V1 contrast matrix of basis in currently
  template <typename TV1>
  Eigen::MatrixXd icTransfer(Eigen::MatrixBase<TV1>& V1){
    MatrixXd V = V1.transpose();
    return V;
  }
  
  //' CLR to ILR (don't bother if you know what your doing)
  //' @param V2 contrast matrix of basis in currently
  template <typename TV2>
  Eigen::MatrixXd ciTransfer(Eigen::MatrixBase<TV2>& V2){
    return V2;
  }
  
  
  //' ILR to ALR 
  //' @param V1 contrast matrix of basis in currently
  //' @param d2 alr reference part to go into
  //' @param D number of total parts
  template <typename T>
  Eigen::MatrixXd iaTransfer(Eigen::MatrixBase<T>& V1, int d2, int D){
    SparseMatrix<double> V2 = alrContrastSparse(d2, D, false);
    MatrixXd V = internal::transferContrast(V1, V2);
    return V;
  }
  
  
  //' CLR to ALR (don't bother if you know what your doing)
  //' @param d2 contrast matrix of basis in currently
  //' @param D total number of parts
  Eigen::MatrixXd caTransfer(int d2, int D){
    MatrixXd V = alrContrast(d2, D, false);
    return V;
  }
  
  //' ALR to CLR
  //' @param d1 contrast matrix of basis in currently
  //' @param D total number of parts
  Eigen::MatrixXd acTransfer(int d1, int D){
    MatrixXd G = alrContrast(d1, D, true);
    G.array() -= 1.0/D;
    return G.transpose();
  }
  
  //' ALR to ALR
  //' @param d1 alr reference part is in currently
  //' @param d2 alr reference part to go into
  //' @param D number of total parts
  Eigen::MatrixXd aaTransfer(int d1, int d2, int D){
    MatrixXd V1 = acTransfer(d1, D);
    SparseMatrix<double> V2 = alrContrastSparse(d2, D, false);
    MatrixXd V = V2*V1; //transferContrast(V1, V2);
    return V;
  }
  
  //' ALR to ILR
  //' @param d1 alr reference part is in currently
  //' @param V2 contrast matrix of basis to go into
  //' @param D number of total parts
  template <typename T>
  Eigen::MatrixXd aiTransfer(int d1, Eigen::MatrixBase<T>& V2, int D){
    MatrixXd V1 = acTransfer(d1, D);
    MatrixXd V = V2*V1;//transferContrast(V1, V2);
    return V;
  }


}


#endif 
