#ifndef FEMUS_TET15_LEVEL_SET_HPP
#define FEMUS_TET15_LEVEL_SET_HPP

#include <vector>
#include <cmath>
#include <stdexcept>

namespace FELevelSetMoments {

inline unsigned Tet15Index(unsigned i, unsigned j, unsigned k) {
  return 25U*i+5U*j+k;
}

inline void Tet15AddProduct(const std::vector<std::vector<long double>>& barycentric, const std::vector<unsigned>& factors, long double multiplier, std::vector<long double>& result) {
  if(multiplier == 0.0L) {
    return;
  }
  std::vector<long double> product(125,0.0L);
  product[0] = 1.0L;
  for(unsigned n = 0; n < factors.size(); ++n) {
    std::vector<long double> next(125,0.0L);
    for(unsigned i = 0; i < 5; ++i) {
      for(unsigned j = 0; j < 5; ++j) {
        for(unsigned k = 0; k < 5; ++k) {
          long double value = product[Tet15Index(i,j,k)];
          if(value == 0.0L) {
            continue;
          }
          for(unsigned a = 0; a+i < 5; ++a) {
            for(unsigned b = 0; b+j < 5; ++b) {
              for(unsigned c = 0; c+k < 5; ++c) {
                next[Tet15Index(i+a,j+b,k+c)] += value*barycentric[factors[n]][Tet15Index(a,b,c)];
              }
            }
          }
        }
      }
    }
    product.swap(next);
  }
  for(unsigned i = 0; i < 125; ++i) {
    result[i] += multiplier*product[i];
  }
}

// Expand the actual FEMuS fifteen-node enriched basis directly in reference or Duffy coordinates.
inline void Tet15PolynomialFromNodalValues(const std::vector<double>& nodalValues, bool duffy, std::vector<long double>& coefficient) {
  if(nodalValues.size() != 15) {
    throw std::invalid_argument("TET15 requires fifteen nodal values in FEMuS order.");
  }
  for(unsigned i = 0; i < 15; ++i) {
    if(!std::isfinite(nodalValues[i])) {
      throw std::invalid_argument("A TET15 nodal value is not finite.");
    }
  }
  std::vector<std::vector<long double>> barycentric(4,std::vector<long double>(125,0.0L));
  barycentric[0][0] = 1.0L;
  barycentric[0][Tet15Index(1,0,0)] = -1.0L;
  if(duffy) {
    barycentric[1][Tet15Index(1,0,0)] = 1.0L;
    barycentric[1][Tet15Index(1,1,0)] = -1.0L;
    barycentric[2][Tet15Index(1,1,0)] = 1.0L;
    barycentric[2][Tet15Index(1,1,1)] = -1.0L;
    barycentric[3][Tet15Index(1,1,1)] = 1.0L;
  }
  else {
    barycentric[0][Tet15Index(0,1,0)] = -1.0L;
    barycentric[0][Tet15Index(0,0,1)] = -1.0L;
    barycentric[1][Tet15Index(1,0,0)] = 1.0L;
    barycentric[2][Tet15Index(0,1,0)] = 1.0L;
    barycentric[3][Tet15Index(0,0,1)] = 1.0L;
  }
  const unsigned edges[6][2] = {{0,1},{1,2},{0,2},{0,3},{1,3},{2,3}};
  const unsigned faces[4][3] = {{0,1,2},{0,1,3},{1,2,3},{0,2,3}};
  coefficient.assign(125,0.0L);
  long double interior = 256.0L*nodalValues[14];
  for(unsigned i = 0; i < 4; ++i) {
    Tet15AddProduct(barycentric,{i,i},2.0L*nodalValues[i],coefficient);
    Tet15AddProduct(barycentric,{i},-static_cast<long double>(nodalValues[i]),coefficient);
    interior -= 4.0L*nodalValues[i];
  }
  for(unsigned edge = 0; edge < 6; ++edge) {
    Tet15AddProduct(barycentric,{edges[edge][0],edges[edge][1]},4.0L*nodalValues[4U+edge],coefficient);
    interior += 32.0L*nodalValues[4U+edge];
  }
  for(unsigned face = 0; face < 4; ++face) {
    long double value = 27.0L*nodalValues[10U+face];
    for(unsigned vertex = 0; vertex < 3; ++vertex) {
      value += 3.0L*nodalValues[faces[face][vertex]];
    }
    for(unsigned edge = 0; edge < 6; ++edge) {
      unsigned contained = 0;
      for(unsigned vertex = 0; vertex < 3; ++vertex) {
        contained += edges[edge][0] == faces[face][vertex] || edges[edge][1] == faces[face][vertex];
      }
      if(contained == 2) {
        value -= 12.0L*nodalValues[4U+edge];
      }
    }
    Tet15AddProduct(barycentric,{faces[face][0],faces[face][1],faces[face][2]},value,coefficient);
    interior -= 108.0L*nodalValues[10U+face];
  }
  Tet15AddProduct(barycentric,{0,1,2,3},interior,coefficient);
}

inline void Tet15ReferenceCoefficientsFromNodalValues(const std::vector<double>& nodalValues, std::vector<std::vector<std::vector<double>>>& coefficient) {
  std::vector<long double> polynomial;
  Tet15PolynomialFromNodalValues(nodalValues,false,polynomial);
  std::vector<std::vector<std::vector<double>>> result(5,std::vector<std::vector<double>>(5,std::vector<double>(5,0.0)));
  for(unsigned i = 0; i < 5; ++i) {
    for(unsigned j = 0; j < 5; ++j) {
      for(unsigned k = 0; k < 5; ++k) {
        result[i][j][k] = static_cast<double>(polynomial[Tet15Index(i,j,k)]);
        if(!std::isfinite(result[i][j][k])) {
          throw std::overflow_error("A TET15 coefficient exceeds double precision range.");
        }
      }
    }
  }
  coefficient.swap(result);
}

}

#endif
