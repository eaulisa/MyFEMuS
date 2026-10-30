#ifndef FEMUS_FE_FIELD_EXTRACTION_HPP
#define FEMUS_FE_FIELD_EXTRACTION_HPP

#include "GeomElTypeEnum.hpp"
#include <limits>
#include <stdexcept>
#include <vector>

// C[i][j][k] is the coefficient of xi^i * eta^j * zeta^k in the original parent coordinates.

inline unsigned FEFieldDimension(const unsigned geom) {
  if(geom == LINE) return 1;
  if(geom == QUAD || geom == TRI) return 2;
  if(geom == HEX || geom == TET || geom == WEDGE) return 3;
  throw std::runtime_error("FEField: unsupported element geometry");
}

inline unsigned FEFieldNodeNumber(const unsigned geom, const unsigned psiType) {
  FEFieldDimension(geom);
  if(psiType > 2) throw std::runtime_error("FEField: only continuous FE types 0, 1 and 2 are supported");
  if(geom == LINE) return psiType == 0 ? 2 : 3;
  if(geom == QUAD) return psiType == 0 ? 4 : (psiType == 1 ? 8 : 9);
  if(geom == TRI) return psiType == 0 ? 3 : (psiType == 1 ? 6 : 7);
  if(geom == HEX) return psiType == 0 ? 8 : (psiType == 1 ? 20 : 27);
  if(geom == TET) return psiType == 0 ? 4 : (psiType == 1 ? 10 : 15);
  return psiType == 0 ? 6 : (psiType == 1 ? 15 : 21);
}

inline void FEFieldReferenceNodes(const unsigned geom, const unsigned psiType, std::vector<std::vector<double>>& nodes) {
  const unsigned count = FEFieldNodeNumber(geom, psiType);
  if(geom == LINE) nodes = {{-1.}, {1.}, {0.}};
  else if(geom == QUAD) nodes = {{-1., -1.}, {1., -1.}, {1., 1.}, {-1., 1.}, {0., -1.}, {1., 0.}, {0., 1.}, {-1., 0.}, {0., 0.}};
  else if(geom == TRI) nodes = {{0., 0.}, {1., 0.}, {0., 1.}, {0.5, 0.}, {0.5, 0.5}, {0., 0.5}, {1. / 3., 1. / 3.}};
  else if(geom == HEX) nodes = {{-1., -1., -1.}, {1., -1., -1.}, {1., 1., -1.}, {-1., 1., -1.}, {-1., -1., 1.}, {1., -1., 1.}, {1., 1., 1.}, {-1., 1., 1.}, {0., -1., -1.}, {1., 0., -1.}, {0., 1., -1.}, {-1., 0., -1.}, {0., -1., 1.}, {1., 0., 1.}, {0., 1., 1.}, {-1., 0., 1.}, {-1., -1., 0.}, {1., -1., 0.}, {1., 1., 0.}, {-1., 1., 0.}, {0., -1., 0.}, {1., 0., 0.}, {0., 1., 0.}, {-1., 0., 0.}, {0., 0., -1.}, {0., 0., 1.}, {0., 0., 0.}};
  else if(geom == TET) nodes = {{0., 0., 0.}, {1., 0., 0.}, {0., 1., 0.}, {0., 0., 1.}, {0.5, 0., 0.}, {0.5, 0.5, 0.}, {0., 0.5, 0.}, {0., 0., 0.5}, {0.5, 0., 0.5}, {0., 0.5, 0.5}, {1. / 3., 1. / 3., 0.}, {1. / 3., 0., 1. / 3.}, {1. / 3., 1. / 3., 1. / 3.}, {0., 1. / 3., 1. / 3.}, {0.25, 0.25, 0.25}};
  else nodes = {{0., 0., -1.}, {1., 0., -1.}, {0., 1., -1.}, {0., 0., 1.}, {1., 0., 1.}, {0., 1., 1.}, {0.5, 0., -1.}, {0.5, 0.5, -1.}, {0., 0.5, -1.}, {0.5, 0., 1.}, {0.5, 0.5, 1.}, {0., 0.5, 1.}, {0., 0., 0.}, {1., 0., 0.}, {0., 1., 0.}, {0.5, 0., 0.}, {0.5, 0.5, 0.}, {0., 0.5, 0.}, {1. / 3., 1. / 3., -1.}, {1. / 3., 1. / 3., 1.}, {1. / 3., 1. / 3., 0.}};
  nodes.resize(count);
}

template <class Type> void FEFieldRequireFinite(const Type& value) {
  if(value != value) throw std::runtime_error("FEField: nonfinite scalar");
  if(std::numeric_limits<Type>::has_infinity && (value == std::numeric_limits<Type>::infinity() || value == -std::numeric_limits<Type>::infinity())) throw std::runtime_error("FEField: nonfinite scalar");
}

template <class Type> void FEFieldAddBarycentricProduct(const unsigned dim, const Type& factor, const std::vector<unsigned>& vertices, std::vector<Type>& coefficients) {
  if(factor == Type(0)) return;
  const unsigned stride[3] = {9, 3, 1};
  std::vector<Type> product(27, Type(0));
  product[0] = factor;
  for(unsigned v = 0; v < vertices.size(); v++) {
    std::vector<Type> next(27, Type(0));
    for(unsigned p = 0; p < 27; p++) {
      if(product[p] == Type(0)) continue;
      if(vertices[v] == 0) next[p] += product[p];
      for(unsigned d = 0; d < dim; d++) {
        if(vertices[v] != 0 && vertices[v] != d + 1) continue;
        if((p / stride[d]) % 3 == 2) throw std::runtime_error("FEField: polynomial exceeds coordinate degree two");
        if(vertices[v] == 0) next[p + stride[d]] -= product[p];
        else next[p + stride[d]] += product[p];
      }
    }
    product.swap(next);
  }
  for(unsigned p = 0; p < 27; p++) coefficients[p] += product[p];
}

template <class Type> void FEFieldSimplexCoefficients(const unsigned dim, const unsigned psiType, const std::vector<Type>& psi, std::vector<Type>& coefficients) {
  const unsigned edges[6][2] = {{0, 1}, {1, 2}, {2, 0}, {0, 3}, {1, 3}, {2, 3}};
  const unsigned faces[4][3] = {{0, 1, 2}, {0, 1, 3}, {1, 2, 3}, {2, 0, 3}};
  const unsigned faceEdges[4][3] = {{0, 1, 2}, {0, 4, 3}, {1, 5, 4}, {2, 3, 5}};
  coefficients.assign(27, Type(0));
  for(unsigned v = 0; v <= dim; v++) {
    if(psiType == 0) FEFieldAddBarycentricProduct<Type>(dim, psi[v], {v}, coefficients);
    else {
      FEFieldAddBarycentricProduct<Type>(dim, Type(2) * psi[v], {v, v}, coefficients);
      FEFieldAddBarycentricProduct<Type>(dim, -psi[v], {v}, coefficients);
    }
  }
  if(psiType == 0) return;
  const unsigned numberOfEdges = dim == 2 ? 3 : 6;
  for(unsigned e = 0; e < numberOfEdges; e++) FEFieldAddBarycentricProduct<Type>(dim, Type(4) * psi[dim + 1 + e], {edges[e][0], edges[e][1]}, coefficients);
  if(psiType == 1) return;
  if(dim == 2) {
    Type bubble = Type(27) * psi[6];
    for(unsigned v = 0; v < 3; v++) bubble += Type(3) * psi[v] - Type(12) * psi[3 + v];
    FEFieldAddBarycentricProduct<Type>(dim, bubble, {0, 1, 2}, coefficients);
  }
  else {
    for(unsigned f = 0; f < 4; f++) {
      Type bubble = Type(27) * psi[10 + f];
      for(unsigned v = 0; v < 3; v++) bubble += Type(3) * psi[faces[f][v]] - Type(12) * psi[4 + faceEdges[f][v]];
      FEFieldAddBarycentricProduct<Type>(dim, bubble, {faces[f][0], faces[f][1], faces[f][2]}, coefficients);
    }
    Type bubble = Type(256) * psi[14];
    for(unsigned v = 0; v < 4; v++) bubble -= Type(4) * psi[v] + Type(108) * psi[10 + v];
    for(unsigned e = 0; e < 6; e++) bubble += Type(32) * psi[4 + e];
    FEFieldAddBarycentricProduct<Type>(dim, bubble, {0, 1, 2, 3}, coefficients);
  }
}

template <class Type> void FEFieldTensorCoefficients(const unsigned geom, const unsigned psiType, const std::vector<Type>& psi, std::vector<Type>& coefficients) {
  const unsigned dim = FEFieldDimension(geom);
  const unsigned corners = 1u << dim;
  const unsigned stride[3] = {9, 3, 1};
  std::vector<std::vector<double>> nodes;
  FEFieldReferenceNodes(geom, 2, nodes);
  std::vector<Type> values(nodes.size(), Type(0));
  for(unsigned n = 0; n < nodes.size(); n++) {
    if(n < psi.size()) {
      values[n] = psi[n];
      continue;
    }
    for(unsigned v = 0; v < psi.size(); v++) {
      Type shape = 1;
      if(v < corners) {
        Type sum = 0;
        for(unsigned d = 0; d < dim; d++) {
          const Type coordinate = Type(nodes[v][d]) * Type(nodes[n][d]);
          shape *= (Type(1) + coordinate) / Type(2);
          sum += coordinate;
        }
        if(psiType == 1) shape *= sum - Type(dim - 1);
      }
      else {
        for(unsigned d = 0; d < dim; d++) {
          if(nodes[v][d] == 0.) shape *= Type(1) - Type(nodes[n][d]) * Type(nodes[n][d]);
          else shape *= (Type(1) + Type(nodes[v][d]) * Type(nodes[n][d])) / Type(2);
        }
      }
      values[n] += psi[v] * shape;
    }
  }
  coefficients.assign(27, Type(0));
  for(unsigned n = 0; n < nodes.size(); n++) {
    unsigned index = 0;
    for(unsigned d = 0; d < dim; d++) index += static_cast<unsigned>(nodes[n][d] + 1.) * stride[d];
    coefficients[index] = values[n];
  }
  const Type transform[3][3] = {{Type(0), Type(1), Type(0)}, {-Type(1) / Type(2), Type(0), Type(1) / Type(2)}, {Type(1) / Type(2), -Type(1), Type(1) / Type(2)}};
  for(unsigned d = 0; d < dim; d++) {
    std::vector<Type> next(27, Type(0));
    for(unsigned p = 0; p < 27; p++) {
      const unsigned row = (p / stride[d]) % 3;
      const unsigned base = p - row * stride[d];
      for(unsigned q = 0; q < 3; q++) next[p] += transform[row][q] * coefficients[base + q * stride[d]];
    }
    coefficients.swap(next);
  }
}

template <class Type> void FEFieldWedgeCoefficients(const unsigned psiType, const std::vector<Type>& psi, std::vector<Type>& coefficients) {
  const unsigned layers[3][7] = {{0, 1, 2, 6, 7, 8, 18}, {3, 4, 5, 9, 10, 11, 19}, {12, 13, 14, 15, 16, 17, 20}};
  const unsigned edges[3][2] = {{0, 1}, {1, 2}, {2, 0}};
  std::vector<Type> values(21, Type(0));
  for(unsigned i = 0; i < psi.size(); i++) values[i] = psi[i];
  if(psiType == 1) {
    for(unsigned e = 0; e < 3; e++) values[15 + e] = -(psi[edges[e][0]] + psi[edges[e][1]] + psi[3 + edges[e][0]] + psi[3 + edges[e][1]]) / Type(4) + (psi[6 + e] + psi[9 + e] + psi[12 + edges[e][0]] + psi[12 + edges[e][1]]) / Type(2);
  }
  coefficients.assign(27, Type(0));
  const unsigned count = psiType == 0 ? 2 : 3;
  const unsigned triangleNodes = psiType == 0 ? 3 : (psiType == 1 ? 6 : 7);
  for(unsigned layer = 0; layer < count; layer++) {
    std::vector<Type> triangleValues(triangleNodes, Type(0));
    for(unsigned i = 0; i < triangleNodes; i++) triangleValues[i] = values[layers[layer][i]];
    std::vector<Type> triangleCoefficients;
    FEFieldSimplexCoefficients<Type>(2, psiType, triangleValues, triangleCoefficients);
    std::vector<Type> vertical(3, Type(0));
    if(psiType == 0) {
      vertical[0] = Type(1) / Type(2);
      vertical[1] = Type(layer == 0 ? -1 : 1) / Type(2);
    }
    else if(layer == 2) {
      vertical[0] = 1;
      vertical[2] = -1;
    }
    else {
      vertical[1] = Type(layer == 0 ? -1 : 1) / Type(2);
      vertical[2] = Type(1) / Type(2);
    }
    for(unsigned i = 0; i < 3; i++) {
      for(unsigned j = 0; j < 3; j++) {
        for(unsigned k = 0; k < 3; k++) coefficients[9 * i + 3 * j + k] += triangleCoefficients[9 * i + 3 * j] * vertical[k];
      }
    }
  }
}

template <class Type> void ExtractFEFieldPolynomial(const unsigned geom, const unsigned psiType, const std::vector<Type>& psi, std::vector<std::vector<std::vector<Type>>>& C) {
  static_assert(std::numeric_limits<Type>::is_specialized && !std::numeric_limits<Type>::is_integer, "FEField requires a real floating-point scalar");
  if(psi.size() != FEFieldNodeNumber(geom, psiType)) throw std::runtime_error("FEField: nodal value count does not match geometry and FE type");
  for(unsigned i = 0; i < psi.size(); i++) FEFieldRequireFinite<Type>(psi[i]);
  std::vector<Type> coefficients;
  if(geom == TRI || geom == TET) FEFieldSimplexCoefficients<Type>(FEFieldDimension(geom), psiType, psi, coefficients);
  else if(geom == WEDGE) FEFieldWedgeCoefficients<Type>(psiType, psi, coefficients);
  else FEFieldTensorCoefficients<Type>(geom, psiType, psi, coefficients);
  for(unsigned i = 0; i < coefficients.size(); i++) FEFieldRequireFinite<Type>(coefficients[i]);
  C.assign(3, std::vector<std::vector<Type>>(3, std::vector<Type>(3, Type(0))));
  for(unsigned i = 0; i < 3; i++) {
    for(unsigned j = 0; j < 3; j++) {
      for(unsigned k = 0; k < 3; k++) C[i][j][k] = coefficients[9 * i + 3 * j + k];
    }
  }
}

template <class Type> void FEFieldCheckPolynomial(const std::vector<std::vector<std::vector<Type>>>& C, const std::vector<Type>& xi) {
  if(xi.empty() || xi.size() > 3 || C.size() != 3) throw std::runtime_error("FEField: invalid polynomial or coordinate dimensions");
  for(unsigned d = 0; d < xi.size(); d++) FEFieldRequireFinite<Type>(xi[d]);
  for(unsigned i = 0; i < 3; i++) {
    if(C[i].size() != 3) throw std::runtime_error("FEField: coefficient array must be 3 x 3 x 3");
    for(unsigned j = 0; j < 3; j++) {
      if(C[i][j].size() != 3) throw std::runtime_error("FEField: coefficient array must be 3 x 3 x 3");
      for(unsigned k = 0; k < 3; k++) {
        FEFieldRequireFinite<Type>(C[i][j][k]);
        if(((xi.size() == 1 && j != 0) || (xi.size() < 3 && k != 0)) && C[i][j][k] != Type(0)) throw std::runtime_error("FEField: insufficient coordinates for this polynomial");
      }
    }
  }
}

template <class Type> void EvaluateFEFieldPolynomial(const std::vector<std::vector<std::vector<Type>>>& C, const std::vector<Type>& xi, Type& value, std::vector<Type>& gradient, std::vector<std::vector<Type>>& hessian) {
  FEFieldCheckPolynomial<Type>(C, xi);
  const unsigned dim = static_cast<unsigned>(xi.size());
  std::vector<std::vector<Type>> powers(3, std::vector<Type>(3, Type(1)));
  for(unsigned d = 0; d < 3; d++) {
    powers[d][1] = d < dim ? xi[d] : Type(0);
    powers[d][2] = powers[d][1] * powers[d][1];
  }
  value = 0;
  gradient.assign(dim, Type(0));
  hessian.assign(dim, std::vector<Type>(dim, Type(0)));
  for(unsigned i = 0; i < 3; i++) {
    for(unsigned j = 0; j < 3; j++) {
      for(unsigned k = 0; k < 3; k++) {
        if(C[i][j][k] == Type(0)) continue;
        const unsigned exponent[3] = {i, j, k};
        value += C[i][j][k] * powers[0][i] * powers[1][j] * powers[2][k];
        for(unsigned d = 0; d < dim; d++) {
          if(exponent[d] == 0) continue;
          Type first = C[i][j][k] * Type(exponent[d]);
          for(unsigned a = 0; a < 3; a++) first *= powers[a][exponent[a] - (a == d ? 1 : 0)];
          gradient[d] += first;
          for(unsigned e = 0; e < dim; e++) {
            const unsigned remaining = exponent[e] - (e == d ? 1 : 0);
            if(remaining == 0) continue;
            Type second = C[i][j][k] * Type(exponent[d]) * Type(remaining);
            for(unsigned a = 0; a < 3; a++) second *= powers[a][exponent[a] - (a == d ? 1 : 0) - (a == e ? 1 : 0)];
            hessian[d][e] += second;
          }
        }
      }
    }
  }
  FEFieldRequireFinite<Type>(value);
  for(unsigned d = 0; d < dim; d++) {
    FEFieldRequireFinite<Type>(gradient[d]);
    for(unsigned e = 0; e < dim; e++) FEFieldRequireFinite<Type>(hessian[d][e]);
  }
}

template <class Type> void RestrictFEFieldToCoordinate(const std::vector<std::vector<std::vector<Type>>>& C, const std::vector<Type>& xi, const unsigned direction, std::vector<Type>& line) {
  FEFieldCheckPolynomial<Type>(C, xi);
  if(direction >= xi.size()) throw std::runtime_error("FEField: invalid restriction direction");
  std::vector<std::vector<Type>> powers(3, std::vector<Type>(3, Type(1)));
  for(unsigned d = 0; d < 3; d++) {
    powers[d][1] = d < xi.size() ? xi[d] : Type(0);
    powers[d][2] = powers[d][1] * powers[d][1];
  }
  line.assign(3, Type(0));
  for(unsigned i = 0; i < 3; i++) {
    for(unsigned j = 0; j < 3; j++) {
      for(unsigned k = 0; k < 3; k++) {
        const unsigned exponent[3] = {i, j, k};
        Type term = C[i][j][k];
        for(unsigned d = 0; d < 3; d++) {
          if(d != direction) term *= powers[d][exponent[d]];
        }
        line[exponent[direction]] += term;
      }
    }
  }
  for(unsigned i = 0; i < 3; i++) FEFieldRequireFinite<Type>(line[i]);
}

#endif
