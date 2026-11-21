#include "CutFemWeight_LS.hpp"
#include <array>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

std::vector<std::vector<long double>> Hex27MomentRows(const CutFemWeight_LS<double,boost::multiprecision::cpp_bin_float_quad>& rule) {
  std::vector<std::vector<long double>> rows(56,std::vector<long double>(rule.GetGaussQuadraturePointNumber()));
  for(unsigned gp = 0; gp < rule.GetGaussQuadraturePointNumber(); ++gp) {
    std::array<std::array<long double,6>,3> powers;
    for(unsigned d = 0; d < 3; ++d) {
      const long double x = (1.0L+rule.GetGaussCoordinatePointer(d)[gp])/2;
      powers[d][0] = 1;
      for(unsigned n = 1; n <= 5; ++n) powers[d][n] = powers[d][n-1]*x;
    }
    unsigned row = 0;
    for(int degree = 0; degree <= 5; ++degree) for(int i = degree; i >= 0; --i) for(int j = degree-i; j >= 0; --j) rows[row++][gp] = rule.GetGaussWeightPointer()[gp]*powers[0][i]*powers[1][j]*powers[2][degree-i-j]/8;
  }
  return rows;
}

}

void TestHex27Elements(const std::string& path) {
  const std::array<unsigned,28> elements = {{7325,7691,23611,62600,62639,62641,62648,62749,62928,62932,63591,63641,63797,63950,63955,63957,64113,64118,67379,67540,67750,67908,68935,69421,69428,69490,69491,69612}};
  const std::array<const char*,5> names = {{"volume","tangential","nx","ny","nz"}};
  std::ifstream input(path);
  if(!input) throw std::runtime_error("Cannot open HEX27 fixture "+path);
  CutFemWeight_LS<double,boost::multiprecision::cpp_bin_float_quad> rule(HEX,5,"legendre");
  const std::vector<std::vector<long double>> rows = Hex27MomentRows(rule);
  for(unsigned element : elements) {
    std::string label;
    unsigned id;
    if(!(input >> label >> id) || label != "HEX27" || id != element) throw std::runtime_error("Missing HEX27 element "+std::to_string(element));
    std::vector<double> psi(27), volume, tangential;
    std::vector<std::vector<double>> coordinates(3,std::vector<double>(27)), normals;
    for(unsigned n = 0; n < 27; ++n) if(!(input >> psi[n] >> coordinates[0][n] >> coordinates[1][n] >> coordinates[2][n])) throw std::runtime_error("Invalid nodes for element "+std::to_string(element));
    if(!(input >> label) || label != "MOMENTS") throw std::runtime_error("Missing moment references for element "+std::to_string(element));
    const unsigned long long warnings = FELevelSetMoments::QuadratureWarningCount();
    rule.GetTangentialAndNormalSurfaceWeightsFromLevelSetData(CutFemLevelSetData(HEX,psi),coordinates,volume,tangential,normals,1e-8,1e-8);
    if(FELevelSetMoments::QuadratureWarningCount() != warnings) throw std::runtime_error("Quadrature accuracy warning on HEX27 element "+std::to_string(element));
    if(normals.size() != 3) throw std::runtime_error("Missing normal weights for element "+std::to_string(element));
    const std::array<const std::vector<double>*,5> weights = {{&volume,&tangential,&normals[0],&normals[1],&normals[2]}};
    for(const std::vector<double>* w : weights) if(w->size() != rule.GetGaussQuadraturePointNumber()) throw std::runtime_error("Incomplete fitted weights for element "+std::to_string(element));
    unsigned row = 0;
    for(int degree = 0; degree <= 5; ++degree) for(int i = degree; i >= 0; --i) for(int j = degree-i; j >= 0; --j, ++row) {
      std::array<long double,5> fitted = {};
      for(unsigned gp = 0; gp < rows[row].size(); ++gp) for(unsigned d = 0; d < 5; ++d) fitted[d] += rows[row][gp]*(*weights[d])[gp];
      for(unsigned d = 0; d < 5; ++d) {
        long double reference;
        if(!(input >> reference) || !std::isfinite(reference)) throw std::runtime_error("Invalid moment reference for element "+std::to_string(element));
        const long double scale = d == 0 ? 1 : (d == 1 ? std::abs(reference) : 1+std::abs(reference));
        if(!std::isfinite(fitted[d]) || std::abs(fitted[d]-reference) > 1e-8L*scale) {
          std::ostringstream message;
          message << std::setprecision(17) << "HEX27 element " << element << ' ' << names[d] << " moment (" << i << ',' << j << ',' << degree-i-j << ") computed " << fitted[d] << ", reference " << reference << ", tolerance " << 1e-8L*scale;
          throw std::runtime_error(message.str());
        }
      }
    }
  }
  std::string extra;
  if(input >> extra) throw std::runtime_error("Unexpected data after the 28 HEX27 elements.");
  std::cout << "28 HEX27 elements passed all 7840 fitted moment checks.\n";
}
