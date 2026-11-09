#ifndef FEMUS_TET15_LEVEL_SET_MOMENTS_HPP
#define FEMUS_TET15_LEVEL_SET_MOMENTS_HPP

#include "Tet15LevelSet.hpp"
#include "AdaptiveLineSliceMoments.hpp"
#include "Hex27LevelSetMoments.hpp"
#include "GraphCellIntegration.hpp"

namespace FELevelSetMoments {

inline unsigned Tet15MomentIndex(unsigned i, unsigned j, unsigned k) {
  unsigned total = i+j+k;
  return total*(total+1U)*(total+2U)/6U+(total-i)*(total-i+1U)/2U+k;
}

inline unsigned Tet15TriangleIndex(unsigned degree, unsigned i, unsigned j) {
  return i*(2U*degree+3U-i)/2U+j;
}

inline void Tet15SliceReferenceMap(const std::vector<long double>& reference, unsigned vertex, std::vector<long double>& field) {
  const unsigned inputSize = SurfacePolynomialCubeSize(reference.size());
  unsigned size = inputSize;
  for(unsigned i = 0; i < inputSize; ++i) {
    for(unsigned j = 0; j < inputSize; ++j) {
      for(unsigned k = 0; k < inputSize; ++k) {
        if(reference[inputSize*inputSize*i+inputSize*j+k] != 0.0L) size = std::max(size,i+j+k+1U);
      }
    }
  }
  field.assign(size*size*size,0.0L);
  for(unsigned i = 0; i < inputSize; ++i) {
    for(unsigned j = 0; j < inputSize; ++j) {
      for(unsigned k = 0; k < inputSize; ++k) {
        long double value = reference[inputSize*inputSize*i+inputSize*j+k];
        if(value == 0.0L) continue;
        if(vertex == 0) {
          for(unsigned a = 0; a <= i; ++a) {
            for(unsigned b = 0; b <= i-a; ++b) {
              for(unsigned c = 0; c <= j; ++c) {
                long double multiplier = Binomial(i,a)*Binomial(i-a,b)*Binomial(j,c);
                if((a+b+c)%2U != 0) multiplier = -multiplier;
                field[size*size*a+size*(b+j-c)+k+c] += multiplier*value;
              }
            }
          }
        }
        else {
          unsigned height = vertex == 1 ? i : (vertex == 2 ? j : k);
          unsigned first = vertex == 1 ? j : i;
          unsigned second = vertex == 3 ? j : k;
          for(unsigned a = 0; a <= first; ++a) {
            long double multiplier = Binomial(first,a);
            if(a%2U != 0) multiplier = -multiplier;
            field[size*size*height+size*(first-a)+second+a] += multiplier*value;
          }
        }
      }
    }
  }
}

inline void Tet15ReferenceSurfaceGradient(const std::vector<double>& nodalValues, std::vector<std::vector<long double>>& gradient, const std::vector<std::vector<double>>* inverseJacobian = nullptr, const std::vector<std::vector<double>>* geometryCoordinates = nullptr) {
  std::vector<long double> coefficient;
  Tet15PolynomialFromNodalValues(nodalValues,false,coefficient);
  gradient.assign(3,std::vector<long double>(125,0.0L));
  for(unsigned i = 0; i < 5; ++i) {
    for(unsigned j = 0; i+j < 5; ++j) {
      for(unsigned k = 0; i+j+k < 5; ++k) {
        const long double value = coefficient[Tet15Index(i,j,k)];
        if(i != 0) gradient[0][Tet15Index(i-1U,j,k)] = i*value;
        if(j != 0) gradient[1][Tet15Index(i,j-1U,k)] = j*value;
        if(k != 0) gradient[2][Tet15Index(i,j,k-1U)] = k*value;
      }
    }
  }
  if(inverseJacobian != nullptr) SurfaceGradientPhysical(gradient,*inverseJacobian);
  if(geometryCoordinates != nullptr) {
    std::vector<std::vector<long double>> coordinates(3);
    for(unsigned component = 0; component < 3; ++component) {
      std::vector<double> nodal(15);
      for(unsigned node = 0; node < 15; ++node) nodal[node] = (*geometryCoordinates)[component][node]-(*geometryCoordinates)[component][0];
      Tet15PolynomialFromNodalValues(nodal,false,coordinates[component]);
    }
    SurfaceGradientCofactor3D(gradient,coordinates,1.0L);
  }
}

inline void Tet15ReferenceSurfaceGradientFields(const std::vector<double>& nodalValues, std::vector<std::vector<std::vector<long double>>>& fields, const std::vector<std::vector<double>>* inverseJacobian = nullptr, const std::vector<std::vector<double>>* geometryCoordinates = nullptr) {
  std::vector<std::vector<long double>> gradient;
  Tet15ReferenceSurfaceGradient(nodalValues,gradient,inverseJacobian,geometryCoordinates);
  fields.assign(4,std::vector<std::vector<long double>>(3));
  for(unsigned vertex = 0; vertex < 4; ++vertex) {
    for(unsigned component = 0; component < 3; ++component) Tet15SliceReferenceMap(gradient[component],vertex,fields[vertex][component]);
  }
}

inline void Tet15SliceReferenceSurfaceGradient(const std::vector<std::vector<long double>>& field, long double height, std::vector<std::vector<std::vector<long double>>>& gradient) {
  gradient.resize(field.size());
  const long double width = 1.0L-height;
  for(unsigned component = 0; component < field.size(); ++component) {
    const unsigned size = SurfacePolynomialCubeSize(field[component].size());
    gradient[component].assign(size,std::vector<long double>(size,0.0L));
    for(unsigned i = 0; i < size; ++i) {
      for(unsigned j = 0; j < size; ++j) {
        long double value = 0.0L;
        for(unsigned h = size; h > 0; --h) value = std::fma(value,height,field[component][size*size*(h-1U)+size*i+j]);
        gradient[component][i][j] = value*Power(width,i+j);
      }
    }
  }
}

inline std::vector<std::vector<std::vector<long double>>> Tet15BuildSliceBasis() {
  std::vector<std::vector<std::vector<long double>>> result(4,std::vector<std::vector<long double>>(15,std::vector<long double>(35,0.0L)));
  const unsigned first[7] = {0,1,0,2,1,0,2};
  const unsigned second[7] = {0,0,1,0,1,2,1};
  for(unsigned node = 0; node < 15; ++node) {
    std::vector<double> values(15,0.0);
    values[node] = 1.0;
    std::vector<long double> reference;
    Tet15PolynomialFromNodalValues(values,false,reference);
    for(unsigned vertex = 0; vertex < 4; ++vertex) {
      std::vector<long double> field;
      Tet15SliceReferenceMap(reference,vertex,field);
      for(unsigned term = 0; term < 7; ++term) {
        for(unsigned h = 0; h < 5; ++h) result[vertex][node][5U*term+h] = field[Tet15Index(h,first[term],second[term])];
      }
      for(unsigned h = 0; h < 5; ++h) {
        for(unsigned i = 0; i < 5; ++i) {
          for(unsigned j = 0; j < 5; ++j) {
            long double expected = 0.0L;
            for(unsigned term = 0; term < 7; ++term) {
              if(i == first[term] && j == second[term]) expected = result[vertex][node][5U*term+h];
            }
            if(i == 1 && j == 2) expected = -result[vertex][node][30U+h];
            if(field[Tet15Index(h,i,j)] != expected) throw std::logic_error("A TET15 slice is outside the enriched TRI7 space.");
          }
        }
      }
    }
  }
  return result;
}

inline const std::vector<std::vector<std::vector<long double>>>& Tet15SliceBasis() {
  static const std::vector<std::vector<std::vector<long double>>> basis = Tet15BuildSliceBasis();
  return basis;
}

inline void Tet15SliceField(const std::vector<double>& nodalValues, unsigned vertex, long double scale, std::vector<std::vector<long double>>& field, std::vector<std::vector<long double>>& unit) {
  const std::vector<std::vector<std::vector<long double>>>& basis = Tet15SliceBasis();
  const unsigned degree[7] = {0,1,1,2,2,2,3};
  field.assign(7,std::vector<long double>(5,0.0L));
  unit.resize(7);
  for(unsigned term = 0; term < 7; ++term) {
    for(unsigned h = 0; h < 5; ++h) {
      for(unsigned node = 0; node < 15; ++node) field[term][h] += basis[vertex][node][5U*term+h]*(static_cast<long double>(nodalValues[node])/scale);
    }
    Trim(field[term]);
    unit[term] = field[term];
    for(unsigned power = 0; power < degree[term]; ++power) unit[term] = PolynomialProduct(unit[term],{1.0L,-1.0L});
  }
}

inline std::vector<long double> Tet15SliceBernstein(const std::vector<std::vector<long double>>& unit) {
  std::vector<std::vector<long double>> power(10,std::vector<long double>(5,0.0L));
  power[Tet15TriangleIndex(3,0,0)] = unit[0];
  power[Tet15TriangleIndex(3,1,0)] = unit[1];
  power[Tet15TriangleIndex(3,0,1)] = PolynomialSum(unit[1],unit[2]);
  power[Tet15TriangleIndex(3,2,0)] = unit[3];
  power[Tet15TriangleIndex(3,1,1)] = PolynomialSum(PolynomialSum(unit[3],unit[3]),unit[4]);
  power[Tet15TriangleIndex(3,0,2)] = PolynomialSum(PolynomialSum(unit[3],unit[4]),unit[5]);
  power[Tet15TriangleIndex(3,2,1)] = unit[6];
  power[Tet15TriangleIndex(3,1,2)] = unit[6];
  std::vector<long double> result(50,0.0L);
  for(unsigned h = 0; h < 5; ++h) {
    for(unsigned i = 0; i <= 3; ++i) {
      for(unsigned j = 0; j <= 3-i; ++j) {
        for(unsigned a = 0; a <= i; ++a) {
          for(unsigned b = 0; b <= j; ++b) {
            const std::vector<long double>& polynomial = power[Tet15TriangleIndex(3,a,b)];
            long double multiplier = Binomial(i,a)*Binomial(j,b)/(Binomial(3,a+b)*Binomial(a+b,a));
            for(unsigned k = 0; k <= h && k < polynomial.size(); ++k) result[10U*h+Tet15TriangleIndex(3,i,j)] += multiplier*polynomial[k]*Binomial(h,k)/Binomial(4,k);
          }
        }
      }
    }
  }
  long double collapsedValue = PolynomialValue(unit[0],1.0L);
  for(unsigned n = 0; n < 10; ++n) result[40U+n] = collapsedValue;
  return result;
}

inline bool Tet15SliceGraph(const std::vector<long double>& bernstein) {
  long double norm = 0.0L;
  for(unsigned n = 0; n < bernstein.size(); ++n) norm = std::max(norm,std::abs(bernstein[n]));
  long double padding = 512.0L*std::numeric_limits<long double>::epsilon()*norm;
  for(unsigned direction = 0; direction < 2; ++direction) {
    bool positive = true;
    bool negative = true;
    for(unsigned i = 0; i <= 2; ++i) {
      for(unsigned j = 0; j <= 2-i; ++j) {
        bool positiveHeight = false;
        bool negativeHeight = false;
        for(unsigned h = 0; h < 5; ++h) {
          long double value = 3.0L*(bernstein[10U*h+Tet15TriangleIndex(3,i+(direction == 0),j+(direction == 1))]-bernstein[10U*h+Tet15TriangleIndex(3,i,j)]);
          positive = positive && value >= -padding;
          negative = negative && value <= padding;
          positiveHeight = positiveHeight || value > padding;
          negativeHeight = negativeHeight || value < -padding;
        }
        positive = positive && positiveHeight;
        negative = negative && negativeHeight;
      }
    }
    if(positive || negative) return true;
  }
  return false;
}

inline void Tet15SplitHeight(const std::vector<long double>& input, long double point, std::vector<long double>& left, std::vector<long double>& right) {
  left.resize(50);
  right.resize(50);
  for(unsigned node = 0; node < 10; ++node) {
    long double value[5];
    for(unsigned h = 0; h < 5; ++h) value[h] = input[10U*h+node];
    left[node] = value[0];
    right[40U+node] = value[4];
    for(unsigned step = 1; step < 5; ++step) {
      for(unsigned h = 0; h < 5-step; ++h) value[h] = (1.0L-point)*value[h]+point*value[h+1U];
      left[10U*step+node] = value[0];
      right[10U*(4U-step)+node] = value[4U-step];
    }
  }
}

inline void Tet15SplitTriangle(const std::vector<long double>& input, unsigned edge, unsigned child, std::vector<long double>& output) {
  unsigned first = edge == 0 ? 0U : (edge == 1 ? 1U : 0U);
  unsigned second = edge == 0 ? 1U : 2U;
  long double vertices[3][3] = {{1.0L,0.0L,0.0L},{0.0L,1.0L,0.0L},{0.0L,0.0L,1.0L}};
  unsigned replaced = child == 0 ? second : first;
  for(unsigned d = 0; d < 3; ++d) vertices[replaced][d] = d == first || d == second ? 0.5L : 0.0L;
  output.resize(50);
  for(unsigned h = 0; h < 5; ++h) {
    for(unsigned i = 0; i <= 3; ++i) {
      for(unsigned j = 0; j <= 3-i; ++j) {
        long double values[10];
        for(unsigned n = 0; n < 10; ++n) values[n] = input[10U*h+n];
        unsigned multiplicity[3] = {3U-i-j,i,j};
        unsigned degree = 3;
        for(unsigned vertex = 0; vertex < 3; ++vertex) {
          for(unsigned repetition = 0; repetition < multiplicity[vertex]; ++repetition) {
            long double next[10] = {0.0L};
            for(unsigned a = 0; a < degree; ++a) {
              for(unsigned b = 0; b < degree-a; ++b) next[Tet15TriangleIndex(degree-1U,a,b)] = vertices[vertex][0]*values[Tet15TriangleIndex(degree,a,b)]+vertices[vertex][1]*values[Tet15TriangleIndex(degree,a+1U,b)]+vertices[vertex][2]*values[Tet15TriangleIndex(degree,a,b+1U)];
            }
            for(unsigned n = 0; n < 10; ++n) values[n] = next[n];
            --degree;
          }
        }
        output[10U*h+Tet15TriangleIndex(3,i,j)] = values[0];
      }
    }
  }
}

inline int Tet15BernsteinPhase(const std::vector<long double>& bernstein, bool computeDelta, unsigned depth, unsigned& budget) {
  long double low = bernstein[0];
  long double high = low;
  for(unsigned n = 1; n < bernstein.size(); ++n) {
    low = std::min(low,bernstein[n]);
    high = std::max(high,bernstein[n]);
  }
  long double padding = 256.0L*std::numeric_limits<long double>::epsilon()*std::max(std::abs(low),std::abs(high));
  if(low > padding || (!computeDelta && low >= 0.0L && high > 0.0L)) return 1;
  if(high < -padding || (!computeDelta && high <= 0.0L && low < 0.0L)) return -1;
  bool strictlyPositiveInside = low >= -padding;
  bool strictlyNegativeInside = high <= padding;
  for(unsigned node = 0; node < 10; ++node) {
    bool positiveHeight = false;
    bool negativeHeight = false;
    for(unsigned h = 0; h < 5; ++h) {
      positiveHeight = positiveHeight || bernstein[10U*h+node] > padding;
      negativeHeight = negativeHeight || bernstein[10U*h+node] < -padding;
    }
    strictlyPositiveInside = strictlyPositiveInside && positiveHeight;
    strictlyNegativeInside = strictlyNegativeInside && negativeHeight;
  }
  if(strictlyPositiveInside) return 1;
  if(strictlyNegativeInside) return -1;
  const unsigned corners[3] = {0,3,9};
  bool positive = false;
  bool negative = false;
  for(unsigned h = 0; h < 5; h += 4) {
    for(unsigned corner = 0; corner < 3; ++corner) {
      positive = positive || bernstein[10U*h+corners[corner]] > padding;
      negative = negative || bernstein[10U*h+corners[corner]] < -padding;
    }
  }
  if((positive && negative) || depth >= 18 || budget == 0) return 0;
  --budget;
  std::vector<long double> left;
  std::vector<long double> right;
  Tet15SplitTriangle(bernstein,depth%3U,0,left);
  int phase = Tet15BernsteinPhase(left,computeDelta,depth+1U,budget);
  if(phase == 0) return 0;
  Tet15SplitTriangle(bernstein,depth%3U,1,right);
  return Tet15BernsteinPhase(right,computeDelta,depth+1U,budget) == phase ? phase : 0;
}

inline int Tet15SlabPhase(const std::vector<long double>& bernstein, long double low, long double high, bool computeDelta) {
  std::vector<long double> discarded;
  std::vector<long double> tail;
  std::vector<long double> slab;
  Tet15SplitHeight(bernstein,low,discarded,tail);
  Tet15SplitHeight(tail,(high-low)/(1.0L-low),slab,discarded);
  unsigned budget = 512;
  return Tet15BernsteinPhase(slab,computeDelta,0,budget);
}

// Lift local triangular moments directly to monomials in the standard FEMuS tetrahedron coordinates.
inline void Tet15LiftMoments(unsigned degree, unsigned vertex, long double height, long double jacobian, const std::vector<long double>& source, const std::vector<long double>& sourceError, std::vector<long double>& values, std::vector<long double>& errors) {
  unsigned count = (degree+1U)*(degree+2U)*(degree+3U)/6U;
  unsigned triangleCount = (degree+1U)*(degree+2U)/2U;
  std::vector<long double> polynomial(count*triangleCount,0.0L);
  polynomial[0] = 1.0L;
  long double barycentric[4][3] = {{0.0L}};
  unsigned remaining[3];
  unsigned position = 0;
  for(unsigned n = 0; n < 4; ++n) {
    if(n != vertex) remaining[position++] = n;
  }
  long double width = 1.0L-height;
  barycentric[vertex][0] = height;
  barycentric[remaining[0]][0] = width;
  barycentric[remaining[0]][1] = -width;
  barycentric[remaining[1]][1] = width;
  barycentric[remaining[1]][2] = -width;
  barycentric[remaining[2]][2] = width;
  long double coordinate[3][3] = {{barycentric[1][0],barycentric[1][1],barycentric[1][2]},{barycentric[2][0],barycentric[2][1],barycentric[2][2]},{barycentric[3][0],barycentric[3][1],barycentric[3][2]}};
  values.assign(2U*count,0.0L);
  errors.assign(2U*count,0.0L);
  for(unsigned total = 0; total <= degree; ++total) {
    for(int first = static_cast<int>(total); first >= 0; --first) {
      for(int second = static_cast<int>(total)-first; second >= 0; --second) {
        unsigned i = static_cast<unsigned>(first);
        unsigned j = static_cast<unsigned>(second);
        unsigned k = total-i-j;
        unsigned index = Tet15MomentIndex(i,j,k);
        if(total != 0) {
          unsigned axis = i != 0 ? 0U : (j != 0 ? 1U : 2U);
          unsigned parent = Tet15MomentIndex(i-(axis == 0),j-(axis == 1),k-(axis == 2));
          for(unsigned a = 0; a < total; ++a) {
            for(unsigned b = 0; b < total-a; ++b) {
              long double value = polynomial[parent*triangleCount+Tet15TriangleIndex(degree,a,b)];
              polynomial[index*triangleCount+Tet15TriangleIndex(degree,a,b)] += coordinate[axis][0]*value;
              polynomial[index*triangleCount+Tet15TriangleIndex(degree,a+1U,b)] += coordinate[axis][1]*value;
              polynomial[index*triangleCount+Tet15TriangleIndex(degree,a,b+1U)] += coordinate[axis][2]*value;
            }
          }
        }
        long double absoluteSum[2] = {0.0L,0.0L};
        long double correction[2] = {0.0L,0.0L};
        for(unsigned a = 0; a <= total; ++a) {
          for(unsigned b = 0; b <= total-a; ++b) {
            const unsigned triangleIndex = Tet15TriangleIndex(degree,a,b);
            const long double multiplier = jacobian*polynomial[index*triangleCount+triangleIndex];
            for(unsigned pair = 0; pair < 2; ++pair) {
              const unsigned destination = pair*count+index;
              const unsigned origin = pair*triangleCount+triangleIndex;
              const long double term = multiplier*source[origin];
              const long double adjusted = term-correction[pair];
              const long double next = values[destination]+adjusted;
              correction[pair] = (next-values[destination])-adjusted;
              values[destination] = next;
              errors[destination] += std::abs(multiplier)*sourceError[origin];
              absoluteSum[pair] += std::abs(term);
            }
          }
        }
        for(unsigned pair = 0; pair < 2; ++pair) errors[pair*count+index] += (8.0L+4.0L*total+2.0L*triangleCount)*std::numeric_limits<long double>::epsilon()*absoluteSum[pair];
      }
    }
  }
}

// Bound the triangular moment lift by its largest absolute row sum.
inline long double Tet15LiftNorm(unsigned degree, long double height, long double jacobian) {
  const long double bound = std::abs(jacobian)*Power(std::max(1.0L,2.0L*(1.0L-height)),degree);
  if(!(bound > 0.0L) || !std::isfinite(bound)) throw std::runtime_error("The TET15 slice has an invalid moment-lift bound.");
  return bound;
}

inline void Tet15AddHeightRoots(std::vector<long double> polynomial, std::vector<long double>& events) {
  Trim(polynomial);
  if(PolynomialNorm(polynomial) == 0.0L) return;
  for(unsigned side = 0; side < 2; ++side) {
    long double root = side == 0 ? 0.0L : 1.0L;
    while(polynomial.size() > 1 && PolynomialValue(polynomial,root) == 0.0L) polynomial = DivideLinear(polynomial,root);
  }
  long double norm = PolynomialNorm(polynomial);
  if(norm == 0.0L) return;
  for(unsigned n = 0; n < polynomial.size(); ++n) polynomial[n] /= norm;
  PolynomialRoots(polynomial,0.0L,1.0L,events);
}

inline std::vector<long double> Tet15HeightEvents(const std::vector<std::vector<long double>>& field, bool regularGraph = false, bool* projectionComplete = nullptr) {
  if(projectionComplete != nullptr) *projectionComplete = true;
  std::vector<long double> events = {0.0L,1.0L};
  const std::vector<long double> width = {1.0L,-1.0L};
  std::vector<long double> square = PolynomialProduct(width,width);
  std::vector<std::vector<std::vector<long double>>> edge(3,std::vector<std::vector<long double>>(3));
  edge[0] = {field[0],field[1],field[3]};
  edge[1][0] = PolynomialSum(PolynomialSum(field[0],PolynomialProduct(field[1],width)),PolynomialProduct(field[3],square));
  edge[1][1] = PolynomialSum(PolynomialSum(field[2],PolynomialProduct(field[4],width)),PolynomialProduct(field[6],square));
  edge[1][2] = PolynomialSum(field[5],PolynomialProduct(field[6],width),-1.0L);
  edge[2] = {field[0],PolynomialSum(field[1],field[2]),PolynomialSum(PolynomialSum(field[3],field[4]),field[5])};
  for(unsigned side = 0; side < 3; ++side) {
    Tet15AddHeightRoots(edge[side][0],events);
    std::vector<long double> candidates;
    Tet15AddHeightRoots(PolynomialSum(PolynomialProduct(edge[side][1],edge[side][1]),PolynomialProduct(edge[side][0],edge[side][2]),-4.0L),candidates);
    for(unsigned n = 0; n < candidates.size(); ++n) {
      long double a = PolynomialValue(edge[side][2],candidates[n]);
      long double b = PolynomialValue(edge[side][1],candidates[n]);
      if(a == 0.0L) continue;
      long double point = -b/(2.0L*a);
      if(point >= 0.0L && point <= 1.0L-candidates[n]) events.push_back(candidates[n]);
    }
    Tet15AddHeightRoots(PolynomialSum(PolynomialSum(edge[side][0],PolynomialProduct(edge[side][1],width)),PolynomialProduct(edge[side][2],square)),events);
  }
  if(regularGraph) {
    SortDistinct(events);
    return events;
  }
  std::vector<std::vector<long double>> a = {field[5],PolynomialSum({0.0L},field[6],-1.0L)};
  std::vector<std::vector<long double>> b = {field[2],field[4],field[6]};
  std::vector<std::vector<long double>> c = {field[0],field[1],field[3]};
  bool quadratic = PolynomialNorm(a[0]) != 0.0L || PolynomialNorm(a[1]) != 0.0L;
  bool linear = PolynomialNorm(b[0]) != 0.0L || PolynomialNorm(b[1]) != 0.0L || PolynomialNorm(b[2]) != 0.0L;
  std::vector<std::vector<long double>> discriminant(5,std::vector<long double>(1,0.0L));
  if(quadratic) {
    for(unsigned i = 0; i < b.size(); ++i) {
      for(unsigned j = 0; j < b.size(); ++j) discriminant[i+j] = PolynomialSum(discriminant[i+j],PolynomialProduct(b[i],b[j]));
    }
    for(unsigned i = 0; i < a.size(); ++i) {
      for(unsigned j = 0; j < c.size(); ++j) discriminant[i+j] = PolynomialSum(discriminant[i+j],PolynomialProduct(a[i],c[j]),-4.0L);
    }
  }
  else if(linear) {
    std::vector<long double> resultant = Hex27Resultant(b,c);
    if(projectionComplete != nullptr && b.size() > 1 && c.size() > 1 && PolynomialNorm(resultant) == 0.0L) *projectionComplete = false;
    Tet15AddHeightRoots(resultant,events);
    discriminant = b;
  }
  else discriminant = c;
  Hex27TrimBivariate(discriminant);
  for(unsigned n = 0; n < discriminant.size(); ++n) Tet15AddHeightRoots(discriminant[n],events);
  if(discriminant.size() > 2) {
    std::vector<std::vector<long double>> derivative(discriminant.size()-1U);
    for(unsigned n = 1; n < discriminant.size(); ++n) {
      derivative[n-1U] = discriminant[n];
      for(unsigned k = 0; k < derivative[n-1U].size(); ++k) derivative[n-1U][k] *= n;
    }
    std::vector<long double> resultant = Hex27Resultant(discriminant,derivative);
    if(projectionComplete != nullptr && PolynomialNorm(resultant) == 0.0L) *projectionComplete = false;
    Tet15AddHeightRoots(resultant,events);
  }
  SortDistinct(events);
  return events;
}

inline std::vector<long double> Tet15FlatHeights(const std::vector<std::vector<long double>>& field) {
  unsigned selected = 7;
  for(unsigned term = 0; term < 7; ++term) {
    if(PolynomialNorm(field[term]) != 0.0L && (selected == 7 || field[term].size() < field[selected].size())) selected = term;
  }
  std::vector<long double> roots;
  std::vector<long double> result;
  if(selected == 7) return result;
  PolynomialRoots(field[selected],0.0L,1.0L,roots);
  for(unsigned n = 0; n < roots.size(); ++n) {
    bool common = true;
    for(unsigned term = 0; term < 7; ++term) common = common && std::abs(PolynomialValue(field[term],roots[n])) <= 1024.0L*std::numeric_limits<long double>::epsilon()*PolynomialBound(field[term],roots[n]);
    if(common && roots[n] < 1.0L) result.push_back(roots[n]);
  }
  SortDistinct(result);
  return result;
}

inline void Tet15RemoveZeroFace(std::vector<std::vector<long double>>& field, unsigned face) {
  std::vector<std::vector<long double>> copy = field;
  field.assign(7,std::vector<long double>(1,0.0L));
  if(face == 0) {
    field[0] = PolynomialSum(PolynomialSum({0.0L},copy[1],-1.0L),PolynomialProduct(copy[3],{1.0L,-1.0L}),-1.0L);
    field[1] = PolynomialSum({0.0L},copy[3],-1.0L);
    field[2] = PolynomialSum(PolynomialSum({0.0L},copy[4],-1.0L),PolynomialProduct(copy[6],{1.0L,-1.0L}),-1.0L);
    field[4] = PolynomialSum({0.0L},copy[6],-1.0L);
    field[5] = copy[6];
  }
  if(face == 1) {
    field[0] = PolynomialSum({0.0L},copy[2],-1.0L);
    field[1] = PolynomialSum(PolynomialSum({0.0L},copy[4],-1.0L),copy[5],-1.0L);
    field[2] = PolynomialSum({0.0L},copy[5],-1.0L);
    field[4] = copy[6];
  }
  if(face == 2) {
    field[0] = copy[2];
    field[1] = copy[4];
    field[2] = copy[5];
    field[3] = copy[6];
    field[4] = PolynomialSum({0.0L},copy[6],-1.0L);
  }
}

inline void Tet15SliceSample(const std::vector<std::vector<long double>>& unit, const std::vector<long double>& heightFactors, const std::vector<unsigned>& zeroBoundary, unsigned degree, unsigned vertex, long double height, long double jacobian, long double absoluteH, long double absoluteDelta, long double relativeTolerance, bool computeDelta, unsigned& evaluations, unsigned maximumEvaluations, std::vector<long double>& values, std::vector<long double>& errors, bool& cut, const std::vector<std::vector<long double>>* surfaceGradient = nullptr, bool* precisionLimited = nullptr) {
  if(evaluations >= maximumEvaluations) throw std::runtime_error("TET15 slice evaluation budget exceeded.");
  if(precisionLimited != nullptr) *precisionLimited = false;
  ++evaluations;
  std::vector<long double> coefficient(7);
  long double factor = 1.0L;
  for(unsigned n = 0; n < heightFactors.size(); ++n) factor *= height-heightFactors[n];
  for(unsigned term = 0; term < 7; ++term) coefficient[term] = factor*PolynomialValue(unit[term],height);
  long double norm = 0.0L;
  for(unsigned term = 0; term < 7; ++term) norm = std::max(norm,std::abs(coefficient[term]));
  int exponent = 0;
  std::frexp(norm,&exponent);
  norm = norm == 0.0L ? 1.0L : std::ldexp(1.0L,exponent-1);
  std::vector<std::vector<long double>> normalized(3,std::vector<long double>(3,0.0L));
  normalized[0][0] = coefficient[0]/norm;
  normalized[1][0] = coefficient[1]/norm;
  normalized[0][1] = coefficient[2]/norm;
  normalized[2][0] = coefficient[3]/norm;
  normalized[1][1] = coefficient[4]/norm;
  normalized[0][2] = coefficient[5]/norm;
  normalized[2][1] = coefficient[6]/norm;
  normalized[1][2] = -coefficient[6]/norm;
  const unsigned triangleCount = (degree+1U)*(degree+2U)/2U;
  std::vector<long double> triangle(2U*triangleCount,0.0L);
  std::vector<long double> triangleError(2U*triangleCount,0.0L);
  std::vector<std::vector<std::vector<long double>>> sliceGradient;
  if(surfaceGradient != nullptr) Tet15SliceReferenceSurfaceGradient(*surfaceGradient,height,sliceGradient);
  const long double measure = jacobian*(1.0L-height)*(1.0L-height);
  long double innerTolerance = absoluteH/Tet15LiftNorm(degree,height,measure);
  const long double deltaBudgetScale = norm*(absoluteDelta/absoluteH);
  if(!(innerTolerance > 0.0L) || !std::isfinite(innerTolerance) || !(deltaBudgetScale > 0.0L) || !std::isfinite(deltaBudgetScale)) throw std::runtime_error("The TET15 inner error allocation exceeds floating point range.");
  int phase = 0;
  LineSliceUtilities::CanonicalGraphOrder(normalized,0,phase);
  if(computeDelta) {
    for(unsigned side = 0; side < zeroBoundary.size(); ++side) {
      if(zeroBoundary[side] != 0) phase = 0;
    }
  }
  unsigned innerEvaluations = 0;
  unsigned innerAccepted = 0;
  unsigned selected = 0;
  long double innerRelative = 0.02L*relativeTolerance;
  long double previousRatio = std::numeric_limits<long double>::infinity();
  long double bestRatio = previousRatio;
  std::vector<long double> bestValues;
  std::vector<long double> bestErrors;
  bool bestCut = false;
  const unsigned count = (degree+1U)*(degree+2U)*(degree+3U)/6U;
  for(unsigned refinement = 0; refinement < 4; ++refinement) {
    if(phase > 0) {
      for(unsigned i = 0; i <= degree; ++i) {
        for(unsigned j = 0; j <= degree-i; ++j) triangle[Tet15TriangleIndex(degree,i,j)] = 1.0L/((i+j+2U)*(j+1U));
      }
    }
    else if(phase == 0) {
      try { TriangleIntegration::Integrate(normalized,degree,{{0.0L,0.0L},{1.0L,0.0L},{1.0L,1.0L}},computeDelta,innerTolerance,innerTolerance*deltaBudgetScale,innerRelative,32,200000,0,10,triangle,triangleError,innerEvaluations,innerAccepted,selected,surfaceGradient == nullptr ? nullptr : &sliceGradient); }
      catch(const std::runtime_error&) {
        if(bestValues.empty()) throw;
        break;
      }
    }
    for(unsigned index = 0; index < triangleCount; ++index) {
      triangle[triangleCount+index] /= norm;
      triangleError[triangleCount+index] /= norm;
      if(!std::isfinite(triangle[triangleCount+index]) || !std::isfinite(triangleError[triangleCount+index])) throw std::overflow_error("A TET15 inner moment exceeds floating point range.");
    }
    Tet15LiftMoments(degree,vertex,height,measure,triangle,triangleError,values,errors);
    long double ratio = 0.0L;
    // Reserve one quarter of the native relative budget for propagated inner errors.
    for(unsigned index = 0; index < values.size(); ++index) {
      if(!std::isfinite(values[index]) || !std::isfinite(errors[index])) throw std::overflow_error("A TET15 lifted moment exceeds floating point range.");
      const long double tolerance = (index < count ? absoluteH : absoluteDelta)+0.25L*relativeTolerance*std::abs(values[index]);
      ratio = std::max(ratio,errors[index]/tolerance);
    }
    if(ratio <= 1.0L) {
      cut = (triangle[0] > 0.0L && triangle[0] < 0.5L) || (computeDelta && triangle[triangleCount] != 0.0L);
      return;
    }
    if(ratio < bestRatio) {
      bestValues = values;
      bestErrors = errors;
      bestRatio = ratio;
      bestCut = (triangle[0] > 0.0L && triangle[0] < 0.5L) || (computeDelta && triangle[triangleCount] != 0.0L);
    }
    if(phase != 0 || ratio >= previousRatio) break;
    std::vector<long double> zeroError(triangleError.size(),0.0L);
    std::vector<long double> rounded;
    std::vector<long double> roundoff;
    Tet15LiftMoments(degree,vertex,height,measure,triangle,zeroError,rounded,roundoff);
    bool limited = false;
    for(unsigned index = 0; index < values.size(); ++index) {
      const long double tolerance = (index < count ? absoluteH : absoluteDelta)+0.25L*relativeTolerance*std::abs(values[index]);
      if(roundoff[index] > tolerance) limited = true;
    }
    if(limited) break;
    previousRatio = ratio;
    innerTolerance /= 2.0L*ratio;
    innerRelative /= 2.0L*ratio;
  }
  if(bestValues.empty()) throw std::runtime_error("The TET15 native slice did not produce a bounded result.");
  values.swap(bestValues);
  errors.swap(bestErrors);
  cut = bestCut;
  if(precisionLimited != nullptr) *precisionLimited = true;
}

inline void Tet15IntegrateInterval(const std::vector<std::vector<long double>>& unit, const std::vector<long double>& heightFactors, const std::vector<long double>& bernstein, const std::vector<unsigned>& zeroBoundary, unsigned degree, unsigned vertex, long double left, long double right, long double parameterLeft, long double parameterRight, long double absoluteH, long double absoluteDelta, long double relativeTolerance, bool computeDelta, bool projectionComplete, unsigned depth, unsigned maximumDepth, unsigned& evaluations, unsigned maximumEvaluations, unsigned& accepted, std::vector<long double>& total, std::vector<long double>& totalError, const std::vector<std::vector<long double>>* surfaceGradient = nullptr) {
  const long double nodes[8] = {0.991455371120812639206854697526329L,0.949107912342758524526189684047851L,0.864864423359769072789712788640926L,0.741531185599394439863864773280788L,0.586087235467691130294144838258730L,0.405845151377397166906606412076961L,0.207784955007898467600689403773245L,0.0L};
  const long double kronrod[8] = {0.022935322010529224963732008058970L,0.063092092629978553290700663189204L,0.104790010322250183839876322541518L,0.140653259715525918745189590510238L,0.169004726639267902826583426598550L,0.190350578064785409913256402421014L,0.204432940075298892414161999234649L,0.209482141084727828012999174891714L};
  const long double gauss[4] = {0.129484966168869693270611432679082L,0.279705391489276667901467771423780L,0.381830050505118944950369775488975L,0.417959183673469387755102040816327L};
  long double width = right-left;
  long double heightLeft = left+width*Power(std::sin(parameterLeft),2);
  long double heightRight = left+width*Power(std::sin(parameterRight),2);
  if(!(heightRight > heightLeft)) throw std::runtime_error("TET15 slicing reached floating point resolution.");
  unsigned count = (degree+1U)*(degree+2U)*(degree+3U)/6U;
  int phase = Tet15SlabPhase(bernstein,heightLeft,heightRight,computeDelta);
  if(phase < 0) {
    ++accepted;
    return;
  }
  long double middle = 0.5L*(parameterLeft+parameterRight);
  long double halfWidth = 0.5L*(parameterRight-parameterLeft);
  const long double fraction = std::max(1.0e-3L,width*(parameterRight-parameterLeft)/(0.5L*std::acos(-1.0L)));
  // Reserve one quarter of this interval's native absolute budget for inner errors.
  const long double innerFraction = 0.25L*fraction/(parameterRight-parameterLeft);
  std::vector<long double> high(2U*count,0.0L);
  std::vector<long double> low(2U*count,0.0L);
  std::vector<long double> innerError(2U*count,0.0L);
  std::vector<long double> absolute(2U*count,0.0L);
  std::vector<std::vector<long double>> samples;
  std::vector<long double> weights;
  bool failed = false;
  bool anyCut = false;
  bool innerLimited = false;
  int uniformConsensus = 2;
  long double split = middle;
  for(unsigned node = 0; node < 8 && !failed; ++node) {
    unsigned sides = node == 7 ? 1U : 2U;
    for(unsigned side = 0; side < sides; ++side) {
      long double parameter = middle+(side == 0 ? -1.0L : 1.0L)*halfWidth*nodes[node];
      long double sine = std::sin(parameter);
      long double height = left+width*sine*sine;
      long double jacobian = 2.0L*width*sine*std::cos(parameter);
      std::vector<long double> value;
      std::vector<long double> error;
      bool sampleCut = false;
      bool sampleLimited = false;
      try {
        if(phase > 0) {
          const unsigned triangleCount = (degree+1U)*(degree+2U)/2U;
          std::vector<long double> triangle(2U*triangleCount,0.0L);
          std::vector<long double> triangleError(2U*triangleCount,0.0L);
          for(unsigned i = 0; i <= degree; ++i) {
            for(unsigned j = 0; j <= degree-i; ++j) triangle[Tet15TriangleIndex(degree,i,j)] = 1.0L/((i+j+2U)*(j+1U));
          }
          Tet15LiftMoments(degree,vertex,height,jacobian*(1.0L-height)*(1.0L-height),triangle,triangleError,value,error);
        }
        else Tet15SliceSample(unit,heightFactors,zeroBoundary,degree,vertex,height,jacobian,absoluteH*innerFraction,absoluteDelta*innerFraction,relativeTolerance,computeDelta,evaluations,maximumEvaluations,value,error,sampleCut,surfaceGradient,&sampleLimited);
      }
      catch(const std::domain_error&) {
        failed = true;
        split = parameter;
        break;
      }
      anyCut = anyCut || sampleCut;
      innerLimited = innerLimited || sampleLimited;
      if(!sampleCut) {
        int samplePhase = value[0] == 0.0L ? -1 : 1;
        if(uniformConsensus == 2) uniformConsensus = samplePhase;
        else if(uniformConsensus != samplePhase) uniformConsensus = 0;
      }
      samples.push_back(value);
      weights.push_back(kronrod[node]);
      for(unsigned k = 0; k < high.size(); ++k) {
        high[k] += kronrod[node]*value[k];
        absolute[k] += kronrod[node]*std::abs(value[k]);
        innerError[k] += kronrod[node]*error[k];
        if(node%2U != 0) low[k] += gauss[node/2U]*value[k];
      }
    }
  }
  // Height events isolate topology changes before a uniformly sampled open interval is accepted.
  bool converged = !failed && (phase != 0 || anyCut || (projectionComplete && (uniformConsensus == 1 || uniformConsensus == -1)));
  bool innerUnresolved = false;
  std::vector<long double> errors(2U*count,0.0L);
  if(!failed) {
    for(unsigned k = 0; k < high.size(); ++k) {
      long double variation = 0.0L;
      for(unsigned node = 0; node < samples.size(); ++node) variation += weights[node]*std::abs(samples[node][k]-0.5L*high[k]);
      long double error = std::abs(high[k]-low[k]);
      if(variation > 0.0L && error > 0.0L) error = variation*std::min(1.0L,std::pow(200.0L*error/variation,1.5L));
      error = std::max(error,50.0L*std::numeric_limits<long double>::epsilon()*absolute[k]);
      errors[k] = halfWidth*(error+innerError[k]);
      high[k] *= halfWidth;
      long double tolerance = (k < count ? absoluteH : absoluteDelta)*fraction+relativeTolerance*std::abs(high[k]);
      if(innerLimited && halfWidth*innerError[k] > tolerance) innerUnresolved = true;
      if(errors[k] > tolerance) converged = false;
    }
  }
  if(converged) {
    for(unsigned k = 0; k < high.size(); ++k) {
      total[k] += high[k];
      totalError[k] += errors[k];
    }
    ++accepted;
    return;
  }
  if(innerUnresolved) throw std::runtime_error("The TET15 inner error exceeds the interval budget at working precision.");
  if(depth >= maximumDepth || !(split > parameterLeft && split < parameterRight)) throw std::runtime_error("TET15 integration did not resolve the requested tolerance or slice topology.");
  Tet15IntegrateInterval(unit,heightFactors,bernstein,zeroBoundary,degree,vertex,left,right,parameterLeft,split,absoluteH,absoluteDelta,relativeTolerance,computeDelta,projectionComplete,depth+1U,maximumDepth,evaluations,maximumEvaluations,accepted,total,totalError,surfaceGradient);
  Tet15IntegrateInterval(unit,heightFactors,bernstein,zeroBoundary,degree,vertex,left,right,split,parameterRight,absoluteH,absoluteDelta,relativeTolerance,computeDelta,projectionComplete,depth+1U,maximumDepth,evaluations,maximumEvaluations,accepted,total,totalError,surfaceGradient);
}

inline void Tet15PlaneSample(const std::vector<long double>& derivative, unsigned degree, unsigned vertex, long double height, long double point, long double absoluteTolerance, std::vector<long double>& value, std::vector<long double>& error, const std::vector<std::vector<std::vector<long double>>>* surfaceGradient = nullptr, long double relativeTolerance = 0.0L, unsigned maximumDepth = 24, unsigned maximumEvaluations = 200000, unsigned* evaluations = nullptr) {
  long double a = derivative[5]-derivative[6]*point;
  long double b = derivative[2]+derivative[4]*point+derivative[6]*point*point;
  long double c = derivative[0]+derivative[1]*point+derivative[3]*point*point;
  long double weight = (1.0L-height)*(1.0L-height);
  if(height == 0.0L) weight *= 0.5L;
  const long double innerTolerance = 0.25L*absoluteTolerance/Tet15LiftNorm(degree,height,weight);
  if(!(innerTolerance > 0.0L) || !std::isfinite(innerTolerance)) throw std::runtime_error("The TET15 plane error allocation exceeds floating point range.");
  std::vector<long double> segment;
  std::vector<long double> segmentError;
  if(surfaceGradient == nullptr) LineSliceUtilities::RationalSegmentMoments(a,b,c,0.0L,point,degree,innerTolerance,0,segment,segmentError);
  else {
    segment.assign(degree+1U,0.0L);
    segmentError.assign(degree+1U,0.0L);
    unsigned innerEvaluations = 0;
    const unsigned budget = evaluations == nullptr ? maximumEvaluations : maximumEvaluations-*evaluations;
    LineSliceSurfaceFiberMoments(a,b,c,point,0.0L,point,degree,*surfaceGradient,innerTolerance,0.02L*relativeTolerance,0,maximumDepth,budget,innerEvaluations,segment,segmentError);
    if(evaluations != nullptr) *evaluations += innerEvaluations;
  }
  const unsigned triangleCount = (degree+1U)*(degree+2U)/2U;
  std::vector<long double> triangle(2U*triangleCount,0.0L);
  std::vector<long double> triangleError(2U*triangleCount,0.0L);
  for(unsigned i = 0; i <= degree; ++i) {
    for(unsigned j = 0; j <= degree-i; ++j) {
      const unsigned index = triangleCount+Tet15TriangleIndex(degree,i,j);
      triangle[index] = Power(point,i)*segment[j];
      triangleError[index] = Power(point,i)*segmentError[j];
    }
  }
  Tet15LiftMoments(degree,vertex,height,weight,triangle,triangleError,value,error);
}

inline void Tet15PlaneInterval(const std::vector<long double>& derivative, unsigned degree, unsigned vertex, long double height, long double left, long double right, long double absoluteTolerance, long double relativeTolerance, unsigned depth, unsigned maximumDepth, unsigned& evaluations, unsigned maximumEvaluations, std::vector<long double>& total, std::vector<long double>& totalError, const std::vector<std::vector<std::vector<long double>>>* surfaceGradient = nullptr) {
  const long double nodes[8] = {0.991455371120812639206854697526329L,0.949107912342758524526189684047851L,0.864864423359769072789712788640926L,0.741531185599394439863864773280788L,0.586087235467691130294144838258730L,0.405845151377397166906606412076961L,0.207784955007898467600689403773245L,0.0L};
  const long double kronrod[8] = {0.022935322010529224963732008058970L,0.063092092629978553290700663189204L,0.104790010322250183839876322541518L,0.140653259715525918745189590510238L,0.169004726639267902826583426598550L,0.190350578064785409913256402421014L,0.204432940075298892414161999234649L,0.209482141084727828012999174891714L};
  const long double gauss[4] = {0.129484966168869693270611432679082L,0.279705391489276667901467771423780L,0.381830050505118944950369775488975L,0.417959183673469387755102040816327L};
  unsigned count = (degree+1U)*(degree+2U)*(degree+3U)/6U;
  std::vector<long double> high(2U*count,0.0L);
  std::vector<long double> low(2U*count,0.0L);
  std::vector<long double> errors(2U*count,0.0L);
  long double middle = 0.5L*(left+right);
  long double halfWidth = 0.5L*(right-left);
  for(unsigned node = 0; node < 8; ++node) {
    for(unsigned side = 0; side < (node == 7 ? 1U : 2U); ++side) {
      if(evaluations >= maximumEvaluations) throw std::runtime_error("TET15 zero-triangle evaluation budget exceeded.");
      ++evaluations;
      long double point = middle+(side == 0 ? -1.0L : 1.0L)*halfWidth*nodes[node];
      std::vector<long double> value;
      std::vector<long double> error;
      Tet15PlaneSample(derivative,degree,vertex,height,point,absoluteTolerance,value,error,surfaceGradient,relativeTolerance,maximumDepth,maximumEvaluations,&evaluations);
      for(unsigned k = count; k < high.size(); ++k) {
        high[k] += kronrod[node]*value[k];
        errors[k] += kronrod[node]*error[k];
        if(node%2U != 0) low[k] += gauss[node/2U]*value[k];
      }
    }
  }
  bool complete = true;
  for(unsigned k = count; k < high.size(); ++k) {
    errors[k] = halfWidth*(errors[k]+4.0L*std::abs(high[k]-low[k]));
    high[k] *= halfWidth;
    if(errors[k] > absoluteTolerance*(right-left)+relativeTolerance*std::abs(high[k])) complete = false;
  }
  if(complete) {
    for(unsigned k = count; k < high.size(); ++k) {
      total[k] += high[k];
      totalError[k] += errors[k];
    }
    return;
  }
  if(depth >= maximumDepth) throw std::runtime_error("TET15 zero-triangle integration did not resolve its requested tolerance.");
  Tet15PlaneInterval(derivative,degree,vertex,height,left,middle,absoluteTolerance,relativeTolerance,depth+1U,maximumDepth,evaluations,maximumEvaluations,total,totalError,surfaceGradient);
  Tet15PlaneInterval(derivative,degree,vertex,height,middle,right,absoluteTolerance,relativeTolerance,depth+1U,maximumDepth,evaluations,maximumEvaluations,total,totalError,surfaceGradient);
}

inline void Tet15PlaneDelta(const std::vector<std::vector<long double>>& field, long double height, unsigned degree, unsigned vertex, long double absoluteTolerance, long double relativeTolerance, unsigned maximumDepth, unsigned& evaluations, unsigned maximumEvaluations, std::vector<long double>& values, std::vector<long double>& errors, bool referenceSurface = false, long double scale = 1.0L, const std::vector<std::vector<double>>* inverseJacobian = nullptr, const std::vector<std::vector<long double>>* surfaceGradient = nullptr) {
  const unsigned powers[7] = {0,1,1,2,2,2,3};
  std::vector<long double> derivative(7);
  std::vector<std::vector<long double>> constant(7);
  for(unsigned term = 0; term < 7; ++term) {
    derivative[term] = PolynomialValue(PolynomialDerivative(field[term]),height)*Power(1.0L-height,powers[term]);
    constant[term] = {derivative[term]};
  }
  std::vector<long double> bernstein = Tet15SliceBernstein(constant);
  unsigned budget = 4096;
  int phase = Tet15BernsteinPhase(bernstein,true,0,budget);
  if(phase == 0) throw std::domain_error("A zero TET15 triangle has a critical or unresolved transverse derivative.");
  if(phase < 0) {
    for(unsigned term = 0; term < 7; ++term) derivative[term] = -derivative[term];
  }
  if(referenceSurface && surfaceGradient == nullptr) {
    const unsigned triangleCount = (degree+1U)*(degree+2U)/2U;
    std::vector<long double> triangle(2U*triangleCount,0.0L);
    std::vector<long double> triangleError(2U*triangleCount,0.0L);
    for(unsigned i = 0; i <= degree; ++i) {
      for(unsigned j = 0; j <= degree-i; ++j) triangle[triangleCount+Tet15TriangleIndex(degree,i,j)] = 1.0L/((i+j+2U)*(j+1U));
    }
    long double normalNorm = vertex == 0 ? std::sqrt(3.0L) : 1.0L;
    if(inverseJacobian != nullptr) {
      std::vector<std::vector<long double>> normal(3,std::vector<long double>(1,vertex == 0 ? -1.0L : 0.0L));
      if(vertex != 0) normal[vertex-1U][0] = 1.0L;
      SurfaceGradientPhysical(normal,*inverseJacobian);
      normalNorm = 0.0L;
      for(unsigned component = 0; component < 3; ++component) normalNorm = std::hypot(normalNorm,normal[component][0]);
    }
    long double weight = scale*normalNorm*(1.0L-height)*(1.0L-height);
    if(height == 0.0L) weight *= 0.5L;
    std::vector<long double> value;
    std::vector<long double> error;
    Tet15LiftMoments(degree,vertex,height,weight,triangle,triangleError,value,error);
    for(unsigned k = 0; k < values.size(); ++k) {
      values[k] += value[k];
      errors[k] += error[k]+128.0L*std::numeric_limits<long double>::epsilon()*std::abs(value[k]);
    }
    return;
  }
  std::vector<std::vector<std::vector<long double>>> sliceGradient;
  if(surfaceGradient != nullptr) Tet15SliceReferenceSurfaceGradient(*surfaceGradient,height,sliceGradient);
  Tet15PlaneInterval(derivative,degree,vertex,height,0.0L,1.0L,absoluteTolerance,relativeTolerance,0,maximumDepth,evaluations,maximumEvaluations,values,errors,surfaceGradient == nullptr ? nullptr : &sliceGradient);
}


inline void Tet15DirectionScores(const std::vector<double>& nodalValues, long double scale, long double* lower, long double* quality) {
  std::vector<long double> polynomial;
  Tet15PolynomialFromNodalValues(nodalValues,false,polynomial);
  for(unsigned n = 0; n < polynomial.size(); ++n) polynomial[n] /= scale;
  const int vertices[4][3] = {{0,0,0},{1,0,0},{0,1,0},{0,0,1}};
  long double edgeLower[4][4] = {{0.0L}};
  long double edgeEnergy[4][4] = {{0.0L}};
  long double totalEnergy = 0.0L;
  for(unsigned first = 0; first < 4; ++first) {
    for(unsigned second = first+1U; second < 4; ++second) {
      std::vector<long double> derivative(125,0.0L);
      long double norm = 0.0L;
      for(unsigned i = 0; i <= 4; ++i) {
        for(unsigned j = 0; j <= 4-i; ++j) {
          for(unsigned k = 0; k <= 4-i-j; ++k) {
            const long double value = polynomial[Tet15Index(i,j,k)];
            if(i != 0) derivative[Tet15Index(i-1U,j,k)] += value*i*(vertices[second][0]-vertices[first][0]);
            if(j != 0) derivative[Tet15Index(i,j-1U,k)] += value*j*(vertices[second][1]-vertices[first][1]);
            if(k != 0) derivative[Tet15Index(i,j,k-1U)] += value*k*(vertices[second][2]-vertices[first][2]);
          }
        }
      }
      for(unsigned n = 0; n < derivative.size(); ++n) norm += std::abs(derivative[n]);
      long double minimum = std::numeric_limits<long double>::infinity();
      long double maximum = -minimum;
      long double energy = 0.0L;
      for(unsigned i = 0; i <= 3; ++i) {
        for(unsigned j = 0; j <= 3-i; ++j) {
          for(unsigned k = 0; k <= 3-i-j; ++k) {
            long double value = 0.0L;
            for(unsigned a = 0; a <= i; ++a) {
              for(unsigned b = 0; b <= j; ++b) {
                for(unsigned c = 0; c <= k; ++c) value += derivative[Tet15Index(a,b,c)]*Binomial(i,a)*Binomial(j,b)*Binomial(k,c)/(Binomial(3,a+b+c)*Binomial(a+b+c,a)*Binomial(b+c,b));
              }
            }
            minimum = std::min(minimum,value);
            maximum = std::max(maximum,value);
            energy += value*value;
          }
        }
      }
      const long double padding = 4096.0L*std::numeric_limits<long double>::epsilon()*norm;
      edgeLower[first][second] = minimum > padding ? minimum-padding : (maximum < -padding ? -maximum-padding : 0.0L);
      edgeEnergy[first][second] = energy;
      totalEnergy += energy;
    }
  }
  for(unsigned vertex = 0; vertex < 4; ++vertex) {
    lower[vertex] = 0.0L;
    quality[vertex] = 0.0L;
    for(unsigned first = 0; first < 4; ++first) {
      for(unsigned second = first+1U; second < 4; ++second) {
        if(first != vertex && second != vertex) {
          lower[vertex] = std::max(lower[vertex],edgeLower[first][second]);
          quality[vertex] += edgeEnergy[first][second];
        }
      }
    }
    quality[vertex] /= totalEnergy+std::numeric_limits<long double>::min();
  }
}

inline bool Tet15DirectionsAgree(const std::vector<long double>& first, const std::vector<long double>& firstError, const std::vector<long double>& second, const std::vector<long double>& secondError, unsigned count, long double absoluteTolerance, long double relativeTolerance, long double scale, bool computeDelta, long double& discrepancy) {
  discrepancy = 0.0L;
  for(unsigned k = 0; k < (computeDelta ? 2U*count : count); ++k) {
    const long double tolerance = (k < count ? absoluteTolerance : absoluteTolerance*scale)+relativeTolerance*std::max(std::abs(first[k]),std::abs(second[k]));
    const long double difference = std::abs(first[k]-second[k]);
    discrepancy = std::max(discrepancy,difference/tolerance);
    if(difference+firstError[k]+secondError[k] > tolerance) return false;
  }
  return true;
}

inline void IntegrateTet15NodalByTriangleSlices(const std::vector<double>& nodalValues, unsigned degree, std::vector<std::vector<std::vector<double>>>& positive, std::vector<std::vector<std::vector<double>>>& delta, double absoluteTolerance = 1.0e-10, double relativeTolerance = 1.0e-10, bool computeDelta = true, unsigned maximumDepth = 24, unsigned integrationAxis = 0, unsigned maximumEvaluations = 200000, std::vector<double>* estimatedError = nullptr, unsigned directionValidation = 1, std::vector<double>* directionInfo = nullptr, bool referenceSurface = false, const std::vector<std::vector<double>>* inverseJacobian = nullptr, const std::vector<std::vector<double>>* geometryCoordinates = nullptr) {
  if(nodalValues.size() != 15 || &positive == &delta || degree > 64 || !(absoluteTolerance > 0.0) || !(relativeTolerance >= 0.0) || !std::isfinite(absoluteTolerance) || !std::isfinite(relativeTolerance) || maximumDepth > 48 || integrationAxis > 4 || maximumEvaluations < 15 || directionValidation > 2 || (directionInfo != nullptr && directionInfo == estimatedError)) throw std::invalid_argument("Invalid adaptive TET15 integration input.");
  if(inverseJacobian != nullptr && !referenceSurface) throw std::invalid_argument("A physical TET15 metric requires surface moments.");
  if(inverseJacobian != nullptr) ValidateSurfaceInverseJacobian(*inverseJacobian,3);
  if(geometryCoordinates != nullptr && (!referenceSurface || inverseJacobian != nullptr)) throw std::invalid_argument("TET15 surface geometry requires surface moments without an inverse Jacobian.");
  if(geometryCoordinates != nullptr) ValidateSurfaceGeometryCoordinates(*geometryCoordinates,3,15);
  long double scale = 0.0L;
  for(unsigned node = 0; node < 15; ++node) {
    if(!std::isfinite(nodalValues[node])) throw std::invalid_argument("A TET15 nodal value is not finite.");
    scale = std::max(scale,std::abs(static_cast<long double>(nodalValues[node])));
  }
  if(scale == 0.0L && computeDelta) throw std::domain_error("The Dirac delta of an identically zero TET15 field is undefined.");
  int exponent = 0;
  std::frexp(scale,&exponent);
  scale = scale == 0.0L ? 1.0L : std::ldexp(1.0L,exponent-1);
  const unsigned boundaryNodes[4][7] = {{1,2,3,5,9,8,12},{0,2,3,6,9,7,13},{0,1,3,4,8,7,11},{0,1,2,4,5,6,10}};
  bool hasZeroFace = false;
  for(unsigned face = 0; face < 4; ++face) {
    bool zero = true;
    for(unsigned node = 0; node < 7; ++node) zero = zero && nodalValues[boundaryNodes[face][node]] == 0.0;
    hasZeroFace = hasZeroFace || zero;
  }
  if(integrationAxis == 0 && directionValidation != 2 && !hasZeroFace) {
    GraphCellData data;
    std::vector<long double> original;
    Tet15PolynomialFromNodalValues(nodalValues,false,original);
    data.coefficient.assign(27,0.0L);
    for(unsigned i = 0; i < 3; ++i) for(unsigned j = 0; j < 3; ++j) for(unsigned k = 0; k < 3; ++k) data.coefficient[Hex27Index(i,j,k)] = original[Tet15Index(i,j,k)]/scale;
    data.degree = degree;
    data.maximumDepth = maximumDepth;
    data.maximumEvaluations = maximumEvaluations;
    data.computeDelta = computeDelta;
    data.absoluteH = absoluteTolerance;
    data.absoluteD = static_cast<long double>(absoluteTolerance)*scale;
    data.relative = relativeTolerance;
    if(referenceSurface && computeDelta) Tet15ReferenceSurfaceGradient(nodalValues,data.gradient,inverseJacobian,geometryCoordinates);
    std::vector<long double> graphValues;
    std::vector<long double> graphErrors;
    GraphCellIntegration::Integrate(data,true,graphValues,graphErrors);
    GraphCellIntegration::Output(data,scale,graphValues,graphErrors,positive,delta,estimatedError);
    if(directionInfo != nullptr) {
      const unsigned attempts = data.innerCounts[0]+data.innerCounts[1]+data.innerCounts[2];
      *directionInfo = {static_cast<double>(attempts),static_cast<double>(attempts),0.0};
    }
    return;
  }
  unsigned count = (degree+1U)*(degree+2U)*(degree+3U)/6U;
  std::vector<std::vector<std::vector<long double>>> fields(4);
  std::vector<std::vector<std::vector<long double>>> units(4);
  std::vector<std::vector<std::vector<long double>>> surfaceGradientFields;
  if(referenceSurface && computeDelta) Tet15ReferenceSurfaceGradientFields(nodalValues,surfaceGradientFields,inverseJacobian,geometryCoordinates);
  long double lower[4];
  long double quality[4];
  unsigned axes[4] = {0,1,2,3};
  Tet15DirectionScores(nodalValues,scale,lower,quality);
  for(unsigned vertex = 0; vertex < 4; ++vertex) Tet15SliceField(nodalValues,vertex,scale,fields[vertex],units[vertex]);
  for(unsigned i = 0; i < 4; ++i) {
    for(unsigned j = i+1U; j < 4; ++j) {
      if(lower[axes[j]] > lower[axes[i]] || (lower[axes[j]] == lower[axes[i]] && quality[axes[j]] > quality[axes[i]])) std::swap(axes[i],axes[j]);
    }
  }
  if(integrationAxis != 0) axes[0] = integrationAxis-1U;
  const unsigned faceNodes[4][7] = {{1,2,3,5,9,8,12},{0,2,3,6,9,7,13},{0,1,3,4,8,7,11},{0,1,2,4,5,6,10}};
  bool zeroFace[4] = {true,true,true,true};
  for(unsigned face = 0; face < 4; ++face) {
    for(unsigned node = 0; node < 7; ++node) zeroFace[face] = zeroFace[face] && nodalValues[faceNodes[face][node]] == 0.0;
  }
  unsigned accepted = 0;
  unsigned evaluations = 0;
  unsigned selected = 0;
  std::vector<long double> values(2U*count,0.0L);
  std::vector<long double> errors(2U*count,0.0L);
  bool complete = false;
  const bool validate = integrationAxis == 0 && (directionValidation == 2 || (directionValidation == 1 && computeDelta && lower[axes[0]] == 0.0L));
  std::vector<std::vector<long double>> previousValues;
  std::vector<std::vector<long double>> previousErrors;
  unsigned successfulDirections = 0;
  unsigned attemptedDirections = 0;
  unsigned totalEvaluations = 0;
  unsigned totalAccepted = 0;
  long double crossDiscrepancy = 0.0L;
  std::string failure;
  for(unsigned attempt = 0; attempt < (integrationAxis == 0 ? 4U : 1U) && !complete; ++attempt) {
    unsigned vertex = axes[attempt];
    ++attemptedDirections;
    evaluations = 0;
    accepted = 0;
    try {
      unsigned remaining[3];
      unsigned position = 0;
      for(unsigned node = 0; node < 4; ++node) {
        if(node != vertex) remaining[position++] = node;
      }
      std::vector<unsigned> zeroBoundary = {static_cast<unsigned>(zeroFace[remaining[2]]),static_cast<unsigned>(zeroFace[remaining[0]]),static_cast<unsigned>(zeroFace[remaining[1]])};
      std::vector<long double> flat = Tet15FlatHeights(fields[vertex]);
      std::vector<std::vector<long double>> residual = fields[vertex];
      for(unsigned n = 0; n < flat.size(); ++n) {
        for(unsigned term = 0; term < 7; ++term) {
          if(residual[term].size() > 1) residual[term] = DivideLinear(residual[term],flat[n]);
          else residual[term].assign(1,0.0L);
        }
      }
      const unsigned slicePowers[7] = {0,1,1,2,2,2,3};
      std::vector<std::vector<long double>> residualUnit = residual;
      for(unsigned term = 0; term < 7; ++term) {
        for(unsigned power = 0; power < slicePowers[term]; ++power) residualUnit[term] = PolynomialProduct(residualUnit[term],{1.0L,-1.0L});
      }
      std::vector<long double> bernstein = Tet15SliceBernstein(residualUnit);
      std::vector<std::vector<long double>> topologyField = residual;
      for(unsigned face = 0; face < 3; ++face) {
        if(zeroFace[remaining[face]]) Tet15RemoveZeroFace(topologyField,face);
      }
      std::vector<std::vector<long double>> topologyUnit = topologyField;
      for(unsigned term = 0; term < 7; ++term) {
        for(unsigned power = 0; power < slicePowers[term]; ++power) topologyUnit[term] = PolynomialProduct(topologyUnit[term],{1.0L,-1.0L});
      }
      std::vector<long double> topologyBernstein = Tet15SliceBernstein(topologyUnit);
      unsigned phaseBudget = 512;
      int topologyPhase = Tet15BernsteinPhase(topologyBernstein,true,0,phaseBudget);
      bool projectionComplete = true;
      std::vector<long double> events = topologyPhase != 0 ? std::vector<long double>{0.0L,1.0L} : Tet15HeightEvents(topologyField,Tet15SliceGraph(topologyBernstein),&projectionComplete);
      for(unsigned n = 0; n < flat.size(); ++n) events.push_back(flat[n]);
      SortDistinct(events);
      std::fill(values.begin(),values.end(),0.0L);
      std::fill(errors.begin(),errors.end(),0.0L);
      accepted = 0;
      evaluations = 0;
      if(computeDelta) {
        for(unsigned n = 0; n < flat.size(); ++n) Tet15PlaneDelta(fields[vertex],flat[n],degree,vertex,0.05L*absoluteTolerance*scale,0.05L*relativeTolerance,maximumDepth,evaluations,maximumEvaluations,values,errors,referenceSurface,scale,inverseJacobian,geometryCoordinates == nullptr ? nullptr : &surfaceGradientFields[vertex]);
      }
      for(unsigned event = 1; event < events.size(); ++event) {
        if(events[event] > events[event-1U]) {
          long double middle = 0.5L*(events[event-1U]+events[event]);
          long double factor = 1.0L;
          for(unsigned n = 0; n < flat.size(); ++n) factor *= middle-flat[n];
          std::vector<long double> signedBernstein = bernstein;
          if(factor < 0.0L) {
            for(unsigned n = 0; n < signedBernstein.size(); ++n) signedBernstein[n] = -signedBernstein[n];
          }
          Tet15IntegrateInterval(residualUnit,flat,signedBernstein,zeroBoundary,degree,vertex,events[event-1U],events[event],0.0L,0.5L*std::acos(-1.0L),absoluteTolerance,static_cast<long double>(absoluteTolerance)*scale,relativeTolerance,computeDelta,projectionComplete,0,maximumDepth,evaluations,maximumEvaluations,accepted,values,errors,referenceSurface && computeDelta ? &surfaceGradientFields[vertex] : nullptr);
        }
      }
      for(unsigned k = 0; k < values.size(); ++k) {
        long double tolerance = (k < count ? absoluteTolerance : static_cast<long double>(absoluteTolerance)*scale)+relativeTolerance*std::abs(values[k]);
        if(errors[k] > tolerance) throw std::runtime_error("The accumulated TET15 error estimate exceeds the requested tolerance.");
      }
      ++successfulDirections;
      totalEvaluations += evaluations;
      totalAccepted += accepted;
      selected = vertex+1U;
      complete = !validate;
      for(unsigned previous = 0; validate && previous < previousValues.size() && !complete; ++previous) {
        long double discrepancy = 0.0L;
        if(Tet15DirectionsAgree(previousValues[previous],previousErrors[previous],values,errors,count,absoluteTolerance,relativeTolerance,scale,computeDelta,discrepancy)) {
          for(unsigned k = 0; k < errors.size(); ++k) errors[k] += std::abs(values[k]-previousValues[previous][k])+previousErrors[previous][k];
          crossDiscrepancy = discrepancy;
          complete = true;
        }
      }
      if(!complete) {
        previousValues.push_back(values);
        previousErrors.push_back(errors);
        failure = "Independent tetrahedral directions did not agree within the requested tolerance.";
      }
    }
    catch(const std::runtime_error& exception) {
      totalEvaluations += evaluations;
      totalAccepted += accepted;
      failure = exception.what();
    }
  }
  if(!complete) throw std::runtime_error("TET15 triangle-slice integration failed: "+failure);
  std::vector<std::vector<std::vector<double>>> positiveOutput(degree+1U);
  std::vector<std::vector<std::vector<double>>> deltaOutput(degree+1U);
  long double maxH = 0.0L;
  long double maxDelta = 0.0L;
  for(unsigned i = 0; i <= degree; ++i) {
    positiveOutput[i].resize(degree-i+1U);
    deltaOutput[i].resize(degree-i+1U);
    for(unsigned j = 0; j <= degree-i; ++j) {
      positiveOutput[i][j].resize(degree-i-j+1U);
      deltaOutput[i][j].resize(degree-i-j+1U);
      for(unsigned k = 0; k <= degree-i-j; ++k) {
        unsigned index = Tet15MomentIndex(i,j,k);
        const long double h = values[index];
        const long double d = computeDelta ? values[count+index]/scale : 0.0L;
        positiveOutput[i][j][k] = static_cast<double>(h);
        deltaOutput[i][j][k] = static_cast<double>(d);
        if(!std::isfinite(positiveOutput[i][j][k]) || !std::isfinite(deltaOutput[i][j][k])) throw std::overflow_error("A TET15 moment exceeds double precision range.");
        const long double errorH = errors[index]+std::abs(static_cast<long double>(positiveOutput[i][j][k])-h);
        const long double errorDelta = errors[count+index]/scale+std::abs(static_cast<long double>(deltaOutput[i][j][k])-d);
        if(!std::isfinite(errorH) || !std::isfinite(errorDelta) || errorH > absoluteTolerance+relativeTolerance*std::abs(h) || errorDelta > absoluteTolerance+relativeTolerance*std::abs(d)) throw std::runtime_error("TET15 output rounding exceeds the requested error budget.");
        maxH = std::max(maxH,errorH);
        maxDelta = std::max(maxDelta,errorDelta);
      }
    }
  }
  positive.swap(positiveOutput);
  delta.swap(deltaOutput);
  if(estimatedError != nullptr) *estimatedError = {static_cast<double>(maxH),static_cast<double>(maxDelta),static_cast<double>(totalAccepted),static_cast<double>(totalEvaluations),static_cast<double>(selected)};
  if(directionInfo != nullptr) *directionInfo = {static_cast<double>(attemptedDirections),static_cast<double>(successfulDirections),static_cast<double>(crossDiscrepancy)};
}

}

#endif

