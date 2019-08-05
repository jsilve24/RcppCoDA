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
  
  
  //' Create ilr default contrast matrix
  //' @param D total number of parts
  //' @return D-1xD contrast matrix (a basis)
  Eigen::MatrixXd ilrContrast(int D){
    if (D < 2) throw std::invalid_argument("D must be > 1");
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


}


#endif 