
#ifndef __femus_GramSchmidt_hpp__
#define __femus_GramSchmidt_hpp__

//#include "CutFemWeight.hpp"

#include <iostream>
#include <iomanip>

//using namespace boost;

//This file Gives the lower triangular matrix A, which gives the polynomial orthonormal basis in terms of the standard basis, or specified basis. All neccassry functions were built in order to complete the task.

template <class Type>
void ModifiedGramSchmidt(const std::vector<std::vector<Type>> &M, std::vector<std::vector<Type>> &A, bool &testIdentity);

template <class Type>
std::vector<std::vector<Type>> Transpose(const std::vector<std::vector<Type>> &M);

template <class TypeO, class TypeA>
std::vector<std::vector<TypeO>> MatrixMatrixMultiply(const std::vector<std::vector<TypeA>> &M, const std::vector<std::vector<TypeA>> &N);

template <class Type>
std::vector<Type> MatrixVectorMultiply(const std::vector<std::vector<Type>> &M, const std::vector<Type> &Nv);

template <class Type>
Type DotProduct(const std::vector<Type> &M, const std::vector<Type> &Nv);

template <class Type>
std::vector<Type> ExtractVector(const std::vector<std::vector<Type>> &M, const unsigned &n, const bool &ROW);

template <class Type>
void GetMassMatrix(const unsigned &element, const unsigned &degree, std::vector<std::vector<Type>> &MM);

template <class Type>
void Get_GS_ATA_Matrix(const GeomElType &geom, const unsigned &d, std::vector<std::vector<Type>> &ATA, bool testIdentity = 0);

template <class Type>
void ModifiedGramSchmidt(const std::vector<std::vector<Type>> &M, std::vector<std::vector<Type>> &A, bool &testIdentity) {

  int size = M.size();
  A.assign(size, std::vector<Type>(size, 0));

  Type DP = 0;
  unsigned c = 0;

  for(int i = 0; i < size; i++) {
    A[i][i] = Type(1);
  }

  Type vN = pow((DotProduct(ExtractVector(A, 0, true), MatrixVectorMultiply(M, ExtractVector(A, 0, true)))), 0.5);
  Type vNd = (DotProduct(ExtractVector(A, 0, true), MatrixVectorMultiply(M, ExtractVector(A, 0, true))));

  for(int i = 0; i < size; i++) {
    if(vN != 0) {
      A[0][i] = A[0][i] / vN;
    }
    else if(vN == 0) {
      std::cout << "division by zero" << std::endl;
    }
  }

  for(unsigned i = 0; i < size - 1; i++) {
    for(unsigned j = i + 1; j < size; j++) {
      DP = DotProduct(ExtractVector(A, j, true), MatrixVectorMultiply(M, ExtractVector(A, i, true)));
      for(unsigned k = 0; k < size; k++) {
        A[j][k] = A[j][k] - (DP * A[i][k]);
      }
    }
    c = i + 1;
    vNd = (DotProduct(ExtractVector(A, c, true), MatrixVectorMultiply(M, ExtractVector(A, c, true))));
    vN = (vNd > 0) ? pow((DotProduct(ExtractVector(A, c, true), MatrixVectorMultiply(M, ExtractVector(A, c, true)))), 0.5) : Type(-1);

    if(vN == -1) std::cout <<  " negative number under sqrt " << std::endl;

    for(int l = 0; l < size; l++) {
      A[c][l] = A[c][l] / vN;
    }
  }

  if(testIdentity) {

    std::cout.precision(14);

    std::vector<std::vector<Type>> Atest(size, std::vector<Type>(size, 0));
    Atest = MatrixMatrixMultiply<Type, Type>(A, MatrixMatrixMultiply<Type, Type>(M, Transpose(A)));

    for(unsigned i = 0; i < size; i++) {
      for(unsigned j = 0; j < size; j++) {
        std::cout << M[i][j] << "  " ;
        if(j ==  M.size() - 1) {
          std::cout << std::endl;
        }
      }
    }

    std::cout << std::endl;

    for(unsigned i = 0; i < size; i++) {
      for(unsigned j = 0; j < size; j++) {
        std::cout << A[i][j] << "  " ;
        if(j ==  A.size() - 1) {
          std::cout << std::endl;
        }
      }
    }

    std::cout << std::endl;

    for(unsigned i = 0; i < size; i++) {
      for(unsigned j = 0; j < size; j++) {
        std::cout << Atest[i][j] << "  " ;
        if(j ==  Atest.size() - 1) {
          std::cout << std::endl;
        }
      }
    }
  }

}

template <class TypeO, class TypeA>
std::vector<std::vector<TypeO>> MatrixMatrixMultiply(const std::vector<std::vector<TypeA>> &M, const std::vector<std::vector<TypeA>> &N) {

  const unsigned &row = M.size();
  const unsigned &columnM = M[0].size();
  const unsigned &rowN = N.size();
  const unsigned &column = N[0].size();

  if(rowN == columnM) {
    TypeA temp;
    std::vector<std::vector<TypeO>> MN(row, std::vector<TypeO>(column, 0));

    for(unsigned i = 0; i < row; i++) {
      for(unsigned j = 0; j < column; j++) {
        temp = 0;
        for(unsigned k = 0; k < rowN; k++) {
          temp += M[i][k] * N[k][j];
        }
        MN[i][j] = static_cast<TypeO>(temp);

      }
    }
    return MN;
  }

  else {
    std::cout << " Dimension mismatch!!!!" << std::endl;
    abort();
  }

}

template <class TypeO, class TypeA>
void MatrixMatrixMultiply(const std::vector<std::vector<TypeA>> &M, const std::vector<std::vector<TypeA>> &N, std::vector<std::vector<TypeO>> &MN) {

  const unsigned &row = M.size();
  const unsigned &columnM = M[0].size();
  const unsigned &rowN = N.size();
  const unsigned &column = N[0].size();

  if(rowN == columnM) {
    TypeA temp;
    MN.assign(row, std::vector<TypeO>(column, 0));
    for(unsigned i = 0; i < row; i++) {
      for(unsigned j = 0; j < column; j++) {
        temp = 0;
        for(unsigned k = 0; k < rowN; k++) {
          temp += M[i][k] * N[k][j];
        }
        MN[i][j] = static_cast<TypeO>(temp);
      }
    }
    return;
  }
  else {
    std::cout << " Dimension mismatch!!!!" << std::endl;
    abort();
  }

}



#include <Eigen/Dense>
template <class TypeA>
void MatrixVectorMultiply(const Eigen::Matrix <TypeA, Eigen::Dynamic, Eigen::Dynamic> &M,
                          const Eigen::Matrix <TypeA, Eigen::Dynamic, 1> &v,
                          Eigen::Matrix <TypeA, Eigen::Dynamic, 1> &u) {
  u = M * v;
}

template <class TypeA>
void MatrixTransposeVectorMultiply(const Eigen::Matrix <TypeA, Eigen::Dynamic, Eigen::Dynamic> &M,
                                   const Eigen::Matrix <TypeA, Eigen::Dynamic, 1> &v,
                                   Eigen::Matrix <TypeA, Eigen::Dynamic, 1> &u) {
  u = M.Transpose() * v;
}


template <class TypeO, class TypeA>
void MatrixVectorMultiply(const std::vector<std::vector<TypeA>> &M, const std::vector<TypeA> &v, std::vector<TypeO> &u) {

  const unsigned &row = M.size();
  const unsigned &column = M[0].size();

  if(column == v.size()) {
    TypeA temp;
    u.assign(row, 0);
    for(unsigned i = 0; i < row; i++) {
      temp = 0;
      for(unsigned j = 0; j < column; j++) {
        temp += M[i][j] * v[j];
      }
      u[i] = static_cast<TypeO>(temp);
    }
    return;
  }
  else {
    std::cout << " Dimension mismatch!!!!" << std::endl;
    abort();
  }

}

template <class TypeO, class TypeA>
void MatrixTransposeVectorMultiply(const std::vector<std::vector<TypeA>> &M, const std::vector<TypeA> &v, std::vector<TypeO> &u) {

  const unsigned &column = M.size();
  const unsigned &row = M[0].size();

  if(column == v.size()) {
    TypeA temp;
    u.assign(row, 0);
    for(unsigned i = 0; i < row; i++) {
      temp = 0;
      for(unsigned j = 0; j < column; j++) {
        temp += M[j][i] * v[j];
      }
      u[i] = static_cast<TypeO>(temp);
    }
    return;
  }
  else {
    std::cout << " Dimension mismatch!!!!" << std::endl;
    abort();
  }

}


template <class Type>
std::vector<Type> MatrixVectorMultiply(const std::vector<std::vector<Type>> &M, const std::vector<Type> &Nv) {

  int row = M.size();
  int column = M[0].size();
  int rowN = Nv.size();
  Type temp = Type(0);
  std::vector<Type> bad(row, Type(0));

  if(column == rowN) {

    std::vector<Type> MNv(row, Type(0));

    for(unsigned i = 0; i < row; i++) {

      for(unsigned j = 0; j < column; j++) {

        temp += M[i][j] * Nv[j];

      }
      MNv[i] = temp;
      temp = Type(0);
    }

    return MNv;
  }


  else {
    std::cout << " Dimension mismatch!!!!" << std::endl;
    return  bad;
  }

}

template <class Type>
Type DotProduct(const std::vector<Type> &M, const std::vector<Type> &Nv) {

  int row = M.size();
  int column = M.size();
  int rowN = Nv.size();
  Type temp = Type(0);
  Type bad = Type(0);

  if(column == rowN) {

    Type MNv;

    for(unsigned i = 0; i < column; i++) {
      temp += M[i] * Nv[i];
    }

    MNv = temp;
    //std::cout << MNv << "  ";

    return MNv;
  }


  else {
    std::cout << " Dimension mismatch!!!!" << std::endl;
    return  bad;
  }

}

template <class Type>
std::vector<Type> ExtractVector(const std::vector<std::vector<Type>> &M, const unsigned &n, const bool &ROW) {

  int row = M.size();
  //std::cout << row << "this is row " << std::endl;
  int column = M[0].size();
  //std::cout << n <<"this is coumn " << std::endl;
  std::vector<Type> bad(row, 0);

  if(n >= row || n >= column) {
    std::cout << "index out of bounds" << std::endl;
    return bad;
  }

  if(ROW) {
    std::vector<Type> result(column, Type(0));

    for(int i = 0; i < column; i++) {
      result[i] = M[n][i];

    }
    return result;

  }

  else {

    std::vector<Type> result(row, 0);

    for(int i = 0; i < row; i++) {

      result[i] = M[i][n];

    }
    return result;

  }

}

template <class Type>
std::vector<std::vector<Type>> Transpose(const std::vector<std::vector<Type>> &M) {

  unsigned rowsize = M.size();
  unsigned columnsize = M[0].size();
  std::vector<std::vector<Type>> MT(columnsize, std::vector<Type>(rowsize, Type(0)));

  for(unsigned i = 0; i < columnsize; i++) {

    for(unsigned j = 0; j < rowsize; j++) {

      MT[i][j] = M[j][i];

    }
  }

  return MT;

}

template <class Type>
void GetMassMatrix(const unsigned &element, const unsigned &degree, std::vector<std::vector<Type>> &MM) {

  //TODO hexahedron = 0, tet = 1, wedge =2, quad = 3, tri = 4, line = 5, point = 6

  if(element == 0 || element == 1 || element == 2) { //3D

    unsigned size = ((degree + 1) * (degree + 2) * (degree + 3)) / 6;

    MM.assign(size, std::vector<Type>(size, Type(0)));
    std::vector<unsigned> m(size, 0);
    std::vector<unsigned> n(size, 0);
    std::vector<unsigned> o(size, 0);

    unsigned count = 0;
    for(int l = 0; l <= degree; l++) {
      for(int i = l; i >= 0; i--) {
        for(int j = l - i; j >= 0; j--) {
          unsigned k = l - i - j;
          m[count] = i;
          n[count] = j;
          o[count] = k;
          count++;
          //std::cout << "i = " << i << " j = " << j << " k = " << k << std::endl;
        }
      }
    }

    if(element == 0) {
      for(int h = 0; h < size; h++) {
        count = 0;
        for(int l = 0; l <= degree; l++) {
          for(int i = l; i >= 0; i--) {
            for(int j = l - i; j >= 0; j--) {
              MM[h][count] = Type(1) / ((m[h] + m[count] + 1) * (n[h] + n[count] + 1) * (o[h] + o[count] + 1));
              count++;
              //std::cout << "i - j = " << i - j << " j = " << j <<  std::endl;
            }
          }
        }
      }
    }

    else if(element == 1) {
      unsigned temp1 = 0;
      unsigned temp2 = 0;
      for(int h = 0; h < size; h++) {
        count = 0;
        for(int l = 0; l <= degree; l++) {
          for(int i = l; i >= 0; i--) {
            for(int j = l - i; j >= 0; j--) {
              MM[h][count] = (Type)(1) / ((3 + m[h] + m[count] + n[h] + n[count] + o[h] + o[count]) * (2 + n[h] + n[count] + o[h] + o[count]) * (1 + o[h] + o[count]));

              count++;
            }
          }
        }
      }
    }

    else if(element == 2) {
      for(int h = 0; h < size; h++) {
        count = 0;
        for(int l = 0; l <= degree; l++) {
          for(int i = l; i >= 0; i--) {
            for(int j = l - i; j >= 0; j--) {
              MM[h][count] = Type(1) / ((m[h] + m[count] + n[h] + n[count] + 2) * (n[h] + n[count] + 1) * (o[h] + o[count] + 1));
              count++;
              //std::cout << "i - j = " << i - j << " j = " << j <<  std::endl;
            }
          }
        }
      }
    }
  } //3D

  else if(element == 3 || element == 4) {  //2D
    unsigned size = ((degree + 1) * (degree + 2)) / 2;
    MM.assign(size, std::vector<Type>(size, Type(0)));
    std::vector<unsigned> m(size, 0);
    std::vector<unsigned> n(size, 0);

    unsigned count = 0;
    for(unsigned q = 0; q <= degree; q++) {
      for(unsigned j = 0; j <= q; j++) {
        unsigned i = q - j;
        m[count] = i;
        n[count] = j;
        count++;
        //std::cout << "i = " << i << " j = " << j <<  std::endl;
      }
    }

    if(element == 3) {
      for(int k = 0; k < size; k++) {
        count = 0;

        for(unsigned q = 0; q <= degree; q++) {
          for(unsigned j = 0; j <= q; j++) {
            unsigned i = q - j;
//         for(int i = 0; i <= degree; i++) {
//           for(int j = i; j >= 0; j--) {
            MM[k][count] = Type(1) / ((m[k] + m[count] + 1) * (n[k] + n[count] + 1));
            count++;
            //std::cout << "i - j = " << i - j << " j = " << j <<  std::endl;
          }
        }
      }
    }


    else if(element == 4) {
      for(int k = 0; k < size; k++) {
        count = 0;
        for(unsigned q = 0; q <= degree; q++) {
          for(unsigned j = 0; j <= q; j++) {
            unsigned i = q - j;

//         for(int i = 0; i <= degree; i++) {
//           for(int j = i; j >= 0; j--) {
            MM[k][count] = Type(1) / ((m[k] + m[count] + n[k] + n[count] + 2) * (n[k] + n[count] + 1));
            count++;
            //std::cout << "i - j = " << i - j << " j = " << j <<  std::endl;
          }
        }
      }
    }
  }  //2D

  else if(element == 5) { //1D
    unsigned size = degree + 1;
    MM.assign(size, std::vector<Type>(size, Type(0)));
    unsigned count = 0;
    for(int i = 0; i < size; i++) {
      for(int j = 0; j < size; j++) {
        MM[i][j] = Type(1) / (i + j + 1);
      }
    }
  }  //1D

  else {
    std::cout << " Not an appropriate element " << std::endl;
    abort;
  }

}

template <class Type>
void Get_GS_ATA_Matrix(const GeomElType &geom, const unsigned &d, std::vector<std::vector<Type>> &ATA, bool testIdentity) {
  std::vector<std::vector<Type>> M;
  GetMassMatrix <Type>(geom, d, M);
  std::vector<std::vector<Type>> A;
  ModifiedGramSchmidt(M, A, testIdentity);
  ATA = MatrixMatrixMultiply<Type, Type> (Transpose(A), A);

  if(testIdentity) {
    std::cout.precision(20);
    std::cout << std::endl;
    for(unsigned i = 0; i < ATA.size(); i++) {
      for(unsigned j = 0; j < ATA[i].size(); j++) {
        std::cout << ATA[i][j] << " ";
      }
      std::cout << std::endl;
    }
  }

  return;
}

#endif



















