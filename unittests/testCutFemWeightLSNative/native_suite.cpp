#include "CutFemWeight_LS.hpp"
#include "Basis.hpp"
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

femus::basis* MakeBasis(unsigned element) {
  if(element == 0) return new femus::LineBiquadratic();
  if(element == 1) return new femus::QuadBiquadratic();
  if(element == 2) return new femus::TriBiquadratic();
  if(element == 3) return new femus::HexBiquadratic();
  return new femus::TetBiquadratic();
}

double Field(const std::string& field, unsigned element, const std::vector<double>& point, std::vector<double>& gradient) {
  const unsigned dimension = point.size();
  const bool simplex = element == 2 || element == 4;
  const double x = point[0];
  const double y = dimension > 1 ? point[1] : 0.;
  const double z = dimension > 2 ? point[2] : 0.;
  gradient.assign(dimension,0.);
  if(field == "full") return 1.;
  if(field == "empty") return -1.;
  if(field == "plane" || field == "plane_y" || field == "boundary_low" || field == "boundary_high" || field == "tiny") {
    const unsigned axis = field == "plane_y" && dimension > 1 ? 1U : 0U;
    double location = simplex ? 0.25 : 0.125;
    if(field == "boundary_low") location = simplex ? 0. : -1.;
    if(field == "boundary_high") location = 1.;
    if(field == "tiny") location = 1. - 1./256.;
    gradient[axis] = 1.;
    return point[axis] - location;
  }
  if(field == "diagonal") {
    double value = -1.5;
    for(unsigned d = 0; d < dimension; ++d) { value += 3. * point[d]; gradient[d] = 3.; }
    return value;
  }
  if(field == "monotone") {
    gradient[0] = 1. + 0.2 * x;
    return x + 0.1 * x * x - 0.25;
  }
  if(field == "two_roots") {
    gradient[0] = 2. * x;
    return x * x - 0.25;
  }
  if(field == "near_affine") {
    const double epsilon = 1./1073741824.;
    gradient[dimension-1U] = 1.;
    gradient[0] += 2. * epsilon * x;
    return point[dimension-1U] - 0.25 + epsilon * x * x;
  }
  if(field == "zero_face") {
    gradient[dimension-1U] = 1. + 0.125 * x;
    gradient[0] += 0.125 * point[dimension-1U];
    return point[dimension-1U] * (1. + 0.125 * x);
  }
  if(field == "sphere") {
    const double center = simplex ? (dimension == 2 ? 1./3. : 0.25) : 0.;
    const double radius = simplex ? (dimension == 2 ? 0.2 : 0.125) : (dimension == 1 ? 0.5 : 0.6);
    double value = radius * radius;
    for(unsigned d = 0; d < dimension; ++d) { value -= (point[d]-center) * (point[d]-center); gradient[d] = -2. * (point[d]-center); }
    return value;
  }
  if(field == "graph") {
    if(element == 0) { gradient[0] = 1. + 0.2 * x; return x + 0.1 * x * x - 0.25; }
    if(element == 1) {
      gradient[0] = -0.4 * x + 0.08 * x * y * y;
      gradient[1] = 1. + 0.08 * x * x * y;
      return y - 0.1 - 0.2 * x * x + 0.04 * x * x * y * y;
    }
    if(element == 2) {
      gradient[0] = -0.3 * x + 0.81 * y * (1. - 2. * x - y);
      gradient[1] = 1. + 0.81 * x * (1. - x - 2. * y);
      return y - 0.2 - 0.15 * x * x + 0.81 * x * y * (1. - x - y);
    }
    if(element == 3) {
      gradient[0] = -0.24 * x + 0.06 * x * z * z;
      gradient[1] = -0.14 * y;
      gradient[2] = 1. + 0.06 * x * x * z;
      return z + 0.05 - 0.12 * x * x - 0.07 * y * y + 0.03 * x * x * z * z;
    }
    gradient[0] = -0.16 * x + 0.02 * y + 0.768 * y * z * (1. - 2. * x - y - z);
    gradient[1] = 0.02 * x + 0.768 * x * z * (1. - x - 2. * y - z);
    gradient[2] = 1. + 0.768 * x * y * (1. - x - y - 2. * z);
    return z - 0.2 - 0.08 * x * x + 0.02 * x * y + 0.768 * x * y * z * (1. - x - y - z);
  }
  throw std::invalid_argument("Unknown verification field.");
}

void Geometry(const std::string& mapping, const std::vector<double>& point, std::vector<double>& coordinates, std::vector<std::vector<double>>& jacobian) {
  const unsigned dimension = point.size();
  coordinates = point;
  jacobian.assign(dimension,std::vector<double>(dimension,0.));
  for(unsigned d = 0; d < dimension; ++d) jacobian[d][d] = 1.;
  if(mapping == "affine") {
    for(unsigned d = 0; d < dimension; ++d) { coordinates[d] = (1.25 + 0.25 * d) * point[d] + 0.2 * (d+1U); jacobian[d][d] = 1.25 + 0.25 * d; }
  }
  else if(mapping == "shear" && dimension > 1) {
    coordinates[1] += 0.15 * point[0] * point[0];
    jacobian[1][0] = 0.3 * point[0];
  }
  else if(mapping == "stretch" || (mapping == "shear" && dimension == 1)) {
    coordinates[0] += 0.05 * point[0] * point[0];
    jacobian[0][0] += 0.1 * point[0];
  }
  else if(mapping != "identity") throw std::invalid_argument("Unknown verification geometry.");
}

double Cofactor(const std::vector<std::vector<double>>& jacobian, const std::vector<double>& gradient, double& determinant) {
  const unsigned dimension = gradient.size();
  if(dimension == 1) { determinant = jacobian[0][0]; return std::abs(gradient[0]); }
  std::vector<std::vector<double>> cofactor(dimension,std::vector<double>(dimension,0.));
  if(dimension == 2) {
    cofactor[0][0] = jacobian[1][1];
    cofactor[0][1] = -jacobian[1][0];
    cofactor[1][0] = -jacobian[0][1];
    cofactor[1][1] = jacobian[0][0];
  }
  else {
    for(unsigned i = 0; i < 3; ++i) {
      for(unsigned j = 0; j < 3; ++j) cofactor[i][j] = jacobian[(i+1U)%3U][(j+1U)%3U] * jacobian[(i+2U)%3U][(j+2U)%3U] - jacobian[(i+1U)%3U][(j+2U)%3U] * jacobian[(i+2U)%3U][(j+1U)%3U];
    }
  }
  determinant = 0.;
  for(unsigned j = 0; j < dimension; ++j) determinant += jacobian[0][j] * cofactor[0][j];
  double norm = 0.;
  for(unsigned i = 0; i < dimension; ++i) {
    double value = 0.;
    for(unsigned j = 0; j < dimension; ++j) value += cofactor[i][j] * gradient[j];
    norm = std::hypot(norm,value);
  }
  return norm;
}

void String(const std::string& value) {
  std::cout << '"';
  for(unsigned i = 0; i < value.size(); ++i) {
    if(value[i] == '"' || value[i] == '\\') std::cout << '\\';
    if(value[i] == '\n') std::cout << "\\n";
    else std::cout << value[i];
  }
  std::cout << '"';
}

void Number(long double value) {
  if(std::isfinite(value)) std::cout << value;
  else std::cout << "null";
}

void Vector(const std::vector<double>& values) {
  std::cout << '[';
  for(unsigned i = 0; i < values.size(); ++i) { if(i != 0) std::cout << ','; Number(values[i]); }
  std::cout << ']';
}

long double Power(long double value, unsigned exponent) {
  long double result = 1.;
  for(unsigned i = 0; i < exponent; ++i) result *= value;
  return result;
}

int main(int argc, char** argv) {
  if(argc != 8) return 64;
  const unsigned element = std::strtoul(argv[1],nullptr,10);
  const unsigned degree = std::strtoul(argv[2],nullptr,10);
  const std::string field = argv[3];
  const std::string mapping = argv[4];
  const double scale = std::strtod(argv[5],nullptr);
  const double tolerance = std::strtod(argv[6],nullptr);
  const unsigned validation = std::strtoul(argv[7],nullptr,10);
  const GeomElType elements[5] = {LINE,QUAD,TRI,HEX,TET};
  try {
    std::cout << std::setprecision(21);
    std::unique_ptr<femus::basis> basis(MakeBasis(element));
    const unsigned dimension = element == 0 ? 1U : (element < 3 ? 2U : 3U);
    std::vector<double> nodal(basis->_nc);
    std::vector<std::vector<double>> coordinates(dimension,std::vector<double>(basis->_nc));
    for(unsigned node = 0; node < nodal.size(); ++node) {
      std::vector<double> point(dimension);
      for(unsigned d = 0; d < dimension; ++d) point[d] = basis->GetXcoarse(node)[d];
      std::vector<double> gradient;
      nodal[node] = scale * Field(field,element,point,gradient);
      std::vector<double> mapped;
      std::vector<std::vector<double>> jacobian;
      Geometry(mapping,point,mapped,jacobian);
      for(unsigned d = 0; d < dimension; ++d) coordinates[d][node] = mapped[d];
    }
    CutFemWeight_LS<double,long double> rule(elements[element],degree,"legendre");
    rule.SetIntegrationDirection(0);
    rule.SetDirectionValidation(validation);
    std::vector<std::vector<double>> h(3), interface(3);
    std::vector<std::string> errors(3);
    for(unsigned mode = 0; mode < 3; ++mode) {
      try {
        if(mode < 2) rule.GetWeightsFromNodalValues(nodal,h[mode],interface[mode],tolerance,tolerance,mode == 1);
        else rule.GetPhysicalSurfaceWeightsFromNodalGeometry(nodal,coordinates,h[mode],interface[mode],tolerance,tolerance);
      }
      catch(const std::exception& error) { errors[mode] = error.what(); h[mode].clear(); interface[mode].clear(); }
    }
    double cacheError = 0.;
    if(errors[2].empty()) {
      std::vector<double> cachedH,cachedS;
      rule.GetPhysicalSurfaceWeightsFromNodalGeometry(nodal,coordinates,cachedH,cachedS,tolerance,tolerance);
      for(unsigned i = 0; i < cachedH.size(); ++i) cacheError = std::max(cacheError,std::max(std::abs(cachedH[i]-h[2][i]),std::abs(cachedS[i]-interface[2][i])));
    }
    const unsigned count = rule.GetGaussQuadraturePointNumber();
    std::vector<std::vector<double>> points(count,std::vector<double>(dimension));
    std::vector<double> norm(count),physicalNorm(count),determinants(count);
    double fieldError = 0.,gradientError = 0.,geometryError = 0.,jacobianError = 0.;
    for(unsigned ig = 0; ig < count; ++ig) {
      for(unsigned d = 0; d < dimension; ++d) points[ig][d] = rule.GetGaussCoordinatePointer(d)[ig];
      std::vector<double> gradient(dimension,0.),position(dimension,0.);
      std::vector<std::vector<double>> jacobian(dimension,std::vector<double>(dimension,0.));
      double value = 0.;
      for(unsigned node = 0; node < nodal.size(); ++node) {
        const double phi = basis->eval_phi(node,points[ig]);
        value += nodal[node] * phi;
        for(unsigned d = 0; d < dimension; ++d) {
          const double derivative = basis->eval_dphidxyz(d,basis->GetIND(node),points[ig].data());
          gradient[d] += nodal[node] * derivative;
          position[d] += coordinates[d][node] * phi;
          for(unsigned physical = 0; physical < dimension; ++physical) jacobian[physical][d] += coordinates[physical][node] * derivative;
        }
      }
      std::vector<double> analyticGradient,analyticPosition;
      std::vector<std::vector<double>> analyticJacobian;
      const double analyticValue = Field(field,element,points[ig],analyticGradient);
      Geometry(mapping,points[ig],analyticPosition,analyticJacobian);
      fieldError = std::max(fieldError,std::abs(value / scale - analyticValue));
      for(unsigned d = 0; d < dimension; ++d) {
        norm[ig] = std::hypot(norm[ig],gradient[d]);
        gradientError = std::max(gradientError,std::abs(gradient[d] / scale - analyticGradient[d]));
        geometryError = std::max(geometryError,std::abs(position[d]-analyticPosition[d]));
        for(unsigned j = 0; j < dimension; ++j) jacobianError = std::max(jacobianError,std::abs(jacobian[d][j]-analyticJacobian[d][j]));
      }
      physicalNorm[ig] = Cofactor(jacobian,gradient,determinants[ig]);
    }
    std::cout << "{\"element\":" << element << ",\"degree\":" << degree << ",\"field\":";
    String(field);
    std::cout << ",\"mapping\":";
    String(mapping);
    std::cout << ",\"scale\":";
    Number(scale);
    std::cout << ",\"gauss_points\":" << count << ",\"errors\":[";
    for(unsigned mode = 0; mode < 3; ++mode) { if(mode != 0) std::cout << ','; String(errors[mode]); }
    std::cout << "],\"cache_error\":";
    Number(cacheError);
    std::cout << ",\"fe_errors\":[";
    Number(fieldError); std::cout << ','; Number(gradientError); std::cout << ','; Number(geometryError); std::cout << ','; Number(jacobianError);
    std::cout << "],\"weights\":{\"h\":"; Vector(h[0]); std::cout << ",\"delta\":"; Vector(interface[0]); std::cout << ",\"surface\":"; Vector(interface[1]); std::cout << ",\"physical\":"; Vector(interface[2]);
    std::cout << "},\"moments\":[";
    bool firstMoment = true;
    for(unsigned total = 0; total <= degree; ++total) {
      for(int ii = total; ii >= 0; --ii) {
        const unsigned i = ii;
        if(dimension == 1 && i != total) continue;
        for(int jj = total-i; jj >= 0; --jj) {
          const unsigned j = jj;
          const unsigned k = total-i-j;
          if(dimension < 3 && k != 0) continue;
          long double sums[7] = {};
          for(unsigned ig = 0; ig < count; ++ig) {
            long double factor = rule.GetGaussWeightPointer()[ig] * Power(points[ig][0],i);
            if(dimension > 1) factor *= Power(points[ig][1],j);
            if(dimension > 2) factor *= Power(points[ig][2],k);
            if(errors[0].empty()) { sums[0] += factor * h[0][ig]; sums[1] += factor * interface[0][ig]; sums[3] += factor * interface[0][ig] * norm[ig]; sums[5] += factor * interface[0][ig] * physicalNorm[ig]; sums[6] += factor * h[0][ig] * std::abs(determinants[ig]); }
            if(errors[1].empty()) sums[2] += factor * interface[1][ig];
            if(errors[2].empty()) sums[4] += factor * interface[2][ig];
          }
          if(!firstMoment) std::cout << ',';
          firstMoment = false;
          std::cout << "{\"powers\":[" << i << ',' << j << ',' << k << "],\"values\":[";
          for(unsigned method = 0; method < 7; ++method) {
            if(method != 0) std::cout << ',';
            const unsigned mode = method == 2 ? 1U : (method == 4 ? 2U : 0U);
            if(errors[mode].empty()) Number(sums[method]);
            else std::cout << "null";
          }
          std::cout << "]}";
        }
      }
    }
    std::cout << "]}\n";
  }
  catch(const std::exception& error) { std::cerr << error.what() << '\n'; return 2; }
  return 0;
}
