#ifndef FEMUS_GRAM_SMITH_LS_HPP
#define FEMUS_GRAM_SMITH_LS_HPP

#include "GeomElTypeEnum.hpp"
#include <cmath>
#include <stdexcept>
#include <vector>

namespace CutFemWeightLSBasis {

template <class Type>
void Orthonormalize(const std::vector<std::vector<Type>> &M, std::vector<std::vector<Type>> &A);

template <class TypeO, class TypeA>
std::vector<std::vector<TypeO>> MatrixMatrixMultiply(const std::vector<std::vector<TypeA>> &M, const std::vector<std::vector<TypeA>> &N);

template <class Type>
std::vector<Type> MatrixVectorMultiply(const std::vector<std::vector<Type>> &M, const std::vector<Type> &Nv);

template <class Type>
Type DotProduct(const std::vector<Type> &M, const std::vector<Type> &Nv);

template <class Type>
std::vector<Type> ExtractRow(const std::vector<std::vector<Type>> &M, const unsigned &n);

template <class Type>
std::vector<std::vector<Type>> Transpose(const std::vector<std::vector<Type>> &M);

template <class Type>
void GetMassMatrix(const unsigned &element, const unsigned &degree, std::vector<std::vector<Type>> &MM);

template <class Type>
void Orthonormalize(const std::vector<std::vector<Type>> &M, std::vector<std::vector<Type>> &A) {

  int size = M.size();
  A.assign(size, std::vector<Type>(size, 0));

  Type DP = 0;
  unsigned c = 0;

  for(int i = 0; i < size; i++) {
    A[i][i] = Type(1);
  }

  Type vN = pow((CutFemWeightLSBasis::DotProduct(CutFemWeightLSBasis::ExtractRow(A, 0), CutFemWeightLSBasis::MatrixVectorMultiply(M, CutFemWeightLSBasis::ExtractRow(A, 0)))), 0.5);
  Type vNd = (CutFemWeightLSBasis::DotProduct(CutFemWeightLSBasis::ExtractRow(A, 0), CutFemWeightLSBasis::MatrixVectorMultiply(M, CutFemWeightLSBasis::ExtractRow(A, 0))));

  for(int i = 0; i < size; i++) {
    if(vN != 0) {
      A[0][i] = A[0][i] / vN;
    }

  }

  for(unsigned i = 0; i < size - 1; i++) {
    for(unsigned j = i + 1; j < size; j++) {
      DP = CutFemWeightLSBasis::DotProduct(CutFemWeightLSBasis::ExtractRow(A, j), CutFemWeightLSBasis::MatrixVectorMultiply(M, CutFemWeightLSBasis::ExtractRow(A, i)));
      for(unsigned k = 0; k < size; k++) {
        A[j][k] = A[j][k] - (DP * A[i][k]);
      }
    }
    c = i + 1;
    vNd = (CutFemWeightLSBasis::DotProduct(CutFemWeightLSBasis::ExtractRow(A, c), CutFemWeightLSBasis::MatrixVectorMultiply(M, CutFemWeightLSBasis::ExtractRow(A, c))));
    vN = (vNd > 0) ? pow((CutFemWeightLSBasis::DotProduct(CutFemWeightLSBasis::ExtractRow(A, c), CutFemWeightLSBasis::MatrixVectorMultiply(M, CutFemWeightLSBasis::ExtractRow(A, c)))), 0.5) : Type(-1);

    for(int l = 0; l < size; l++) {
      A[c][l] = A[c][l] / vN;
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
    throw std::invalid_argument("Incompatible basis matrix dimensions.");
  }

}

template <class Type>
std::vector<Type> MatrixVectorMultiply(const std::vector<std::vector<Type>> &M, const std::vector<Type> &Nv) {

  int row = M.size();
  int column = M[0].size();
  int rowN = Nv.size();
  Type temp = Type(0);

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
    throw std::invalid_argument("Incompatible basis matrix dimensions.");
  }

}

template <class Type>
Type DotProduct(const std::vector<Type> &M, const std::vector<Type> &Nv) {

  int column = M.size();
  int rowN = Nv.size();
  Type temp = Type(0);

  if(column == rowN) {

    Type MNv;

    for(unsigned i = 0; i < column; i++) {
      temp += M[i] * Nv[i];
    }

    MNv = temp;

    return MNv;
  }

  else {
    throw std::invalid_argument("Incompatible basis matrix dimensions.");
  }

}

template <class Type>
std::vector<Type> ExtractRow(const std::vector<std::vector<Type>> &M, const unsigned &n) {

  int row = M.size();

  if(n >= row) {
    throw std::out_of_range("Invalid basis matrix index.");
  }

  int column = M[0].size();
  std::vector<Type> result(column, Type(0));
  for(int i = 0; i < column; i++) {
    result[i] = M[n][i];
  }
  return result;
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

  if(element == 0 || element == 1) {

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

            }
          }
        }
      }
    }

    else if(element == 1) {
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

  }

  else if(element == 3 || element == 4) {
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

      }
    }

    if(element == 3) {
      for(int k = 0; k < size; k++) {
        count = 0;

        for(unsigned q = 0; q <= degree; q++) {
          for(unsigned j = 0; j <= q; j++) {
            MM[k][count] = Type(1) / ((m[k] + m[count] + 1) * (n[k] + n[count] + 1));
            count++;

          }
        }
      }
    }

    else if(element == 4) {
      for(int k = 0; k < size; k++) {
        count = 0;
        for(unsigned q = 0; q <= degree; q++) {
          for(unsigned j = 0; j <= q; j++) {
            MM[k][count] = Type(1) / ((m[k] + m[count] + n[k] + n[count] + 2) * (n[k] + n[count] + 1));
            count++;

          }
        }
      }
    }
  }

  else if(element == 5) {
    unsigned size = degree + 1;
    MM.assign(size, std::vector<Type>(size, Type(0)));
    for(int i = 0; i < size; i++) {
      for(int j = 0; j < size; j++) {
        MM[i][j] = Type(1) / (i + j + 1);
      }
    }
  }

  else {
    throw std::invalid_argument("Unsupported level-set element.");
  }

}

template <class Type>
void Get_GS_ATA_Matrix(const GeomElType &geom, const unsigned &d, std::vector<std::vector<Type>> &ATA) {
  std::vector<std::vector<Type>> M;
  CutFemWeightLSBasis::GetMassMatrix <Type>(geom, d, M);
  std::vector<std::vector<Type>> A;
  CutFemWeightLSBasis::Orthonormalize(M, A);
  ATA = CutFemWeightLSBasis::MatrixMatrixMultiply<Type, Type> (CutFemWeightLSBasis::Transpose(A), A);

  return;
}

}

#endif
