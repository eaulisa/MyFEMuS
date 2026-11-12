#include "CutFemWeight_LS.hpp"
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <type_traits>
#include <utility>
#include <vector>

template <class T> auto HasBuild(int) -> decltype(std::declval<T&>().build(),std::true_type());
template <class T> std::false_type HasBuild(...);
static_assert(!decltype(HasBuild<CutFemWeight_LS<double,long double>>(0))::value,"Removed API Build remains public.");

template <class T> auto HasClear(int) -> decltype(std::declval<T&>().clear(),std::true_type());
template <class T> std::false_type HasClear(...);
static_assert(!decltype(HasClear<CutFemWeight_LS<double,long double>>(0))::value,"Removed API Clear remains public.");

template <class T> auto HasClearMap(int) -> decltype(std::declval<T&>().ClearMap(),std::true_type());
template <class T> std::false_type HasClearMap(...);
static_assert(!decltype(HasClearMap<CutFemWeight_LS<double,long double>>(0))::value,"Removed API ClearMap remains public.");

template <class T> auto HasCounter(int) -> decltype(std::declval<T&>().GetCounter(),std::true_type());
template <class T> std::false_type HasCounter(...);
static_assert(!decltype(HasCounter<CutFemWeight_LS<double,long double>>(0))::value,"Removed API Counter remains public.");

template <class T> auto HasOrderUpdate(int) -> decltype(std::declval<T&>().UpdateQuadratureRule(1U),std::true_type());
template <class T> std::false_type HasOrderUpdate(...);
static_assert(!decltype(HasOrderUpdate<CutFemWeight_LS<double,long double>>(0))::value,"Removed API OrderUpdate remains public.");

template <class T> auto HasPlanarOperator(int) -> decltype(std::declval<T&>()(0,std::declval<std::vector<double>&>(),0.,std::declval<std::vector<double>&>()),std::true_type());
template <class T> std::false_type HasPlanarOperator(...);
static_assert(!decltype(HasPlanarOperator<CutFemWeight_LS<double,long double>>(0))::value,"Removed API PlanarOperator remains public.");

template <class T> auto HasParabolaOperator(int) -> decltype(std::declval<T&>()(0,std::declval<std::vector<long double>&>(),0.L,true,std::declval<std::vector<double>&>()),std::true_type());
template <class T> std::false_type HasParabolaOperator(...);
static_assert(!decltype(HasParabolaOperator<CutFemWeight_LS<double,long double>>(0))::value,"Removed API ParabolaOperator remains public.");

template <class T> auto HasLineMoments(int) -> decltype(std::declval<T&>().GetWeightFromMoments(std::declval<std::vector<double>&>(),std::declval<std::vector<double>&>()),std::true_type());
template <class T> std::false_type HasLineMoments(...);
static_assert(!decltype(HasLineMoments<CutFemWeight_LS<double,long double>>(0))::value,"Removed API LineMoments remains public.");

template <class T> auto HasSurfaceMoments(int) -> decltype(std::declval<T&>().GetWeightFromMoments(std::declval<std::vector<std::vector<double>>&>(),std::declval<std::vector<double>&>()),std::true_type());
template <class T> std::false_type HasSurfaceMoments(...);
static_assert(!decltype(HasSurfaceMoments<CutFemWeight_LS<double,long double>>(0))::value,"Removed API SurfaceMoments remains public.");

template <class T> auto HasVolumeMoments(int) -> decltype(std::declval<T&>().GetWeightFromMoments(std::declval<std::vector<std::vector<std::vector<double>>>&>(),std::declval<std::vector<double>&>()),std::true_type());
template <class T> std::false_type HasVolumeMoments(...);
static_assert(!decltype(HasVolumeMoments<CutFemWeight_LS<double,long double>>(0))::value,"Removed API VolumeMoments remains public.");

template <class T> auto HasFECoefficients(int) -> decltype(std::declval<T&>().GetWeightFromFE(0,std::declval<std::vector<double>&>(),std::declval<std::vector<double>&>()),std::true_type());
template <class T> std::false_type HasFECoefficients(...);
static_assert(!decltype(HasFECoefficients<CutFemWeight_LS<double,long double>>(0))::value,"Removed API FECoefficients remains public.");

template <class T> auto HasFEQuadCoefficients(int) -> decltype(std::declval<T&>().GetWeightFromFE(0,std::declval<std::vector<std::vector<double>>&>(),std::declval<std::vector<double>&>()),std::true_type());
template <class T> std::false_type HasFEQuadCoefficients(...);
static_assert(!decltype(HasFEQuadCoefficients<CutFemWeight_LS<double,long double>>(0))::value,"Removed API FEQuadCoefficients remains public.");

template <class T> auto HasCanonicalCoefficients(int) -> decltype(std::declval<T&>().GetWeightFromCanonicalFE(0,std::declval<std::vector<double>&>(),std::declval<std::vector<double>&>()),std::true_type());
template <class T> std::false_type HasCanonicalCoefficients(...);
static_assert(!decltype(HasCanonicalCoefficients<CutFemWeight_LS<double,long double>>(0))::value,"Removed API CanonicalCoefficients remains public.");

template <class T> auto HasCanonicalQuadCoefficients(int) -> decltype(std::declval<T&>().GetWeightFromCanonicalFE(0,std::declval<std::vector<std::vector<double>>&>(),std::declval<std::vector<double>&>()),std::true_type());
template <class T> std::false_type HasCanonicalQuadCoefficients(...);
static_assert(!decltype(HasCanonicalQuadCoefficients<CutFemWeight_LS<double,long double>>(0))::value,"Removed API CanonicalQuadCoefficients remains public.");

template <class T> auto HasCanonicalHexCoefficients(int) -> decltype(std::declval<T&>().GetWeightFromCanonicalFE(0,std::declval<std::vector<std::vector<std::vector<double>>>&>(),std::declval<std::vector<double>&>()),std::true_type());
template <class T> std::false_type HasCanonicalHexCoefficients(...);
static_assert(!decltype(HasCanonicalHexCoefficients<CutFemWeight_LS<double,long double>>(0))::value,"Removed API CanonicalHexCoefficients remains public.");

template <class T> auto HasSingleNodalWeights(int) -> decltype(std::declval<T&>().GetWeightFromNodalValues(0,std::declval<std::vector<double>&>(),std::declval<std::vector<double>&>()),std::true_type());
template <class T> std::false_type HasSingleNodalWeights(...);
static_assert(!decltype(HasSingleNodalWeights<CutFemWeight_LS<double,long double>>(0))::value,"Removed API SingleNodalWeights remains public.");

template <class T> auto HasLineSliceNodalWeights(int) -> decltype(std::declval<T&>().GetWeightsFromNodalValuesByLineSlices(std::declval<std::vector<double>&>(),std::declval<std::vector<double>&>(),std::declval<std::vector<double>&>()),std::true_type());
template <class T> std::false_type HasLineSliceNodalWeights(...);
static_assert(!decltype(HasLineSliceNodalWeights<CutFemWeight_LS<double,long double>>(0))::value,"Removed API LineSliceNodalWeights remains public.");

template <class T> auto HasAffineSurface(int) -> decltype(std::declval<T&>().GetPhysicalSurfaceWeightsFromNodalValues(std::declval<std::vector<double>&>(),std::declval<std::vector<std::vector<double>>&>(),std::declval<std::vector<double>&>(),std::declval<std::vector<double>&>()),std::true_type());
template <class T> std::false_type HasAffineSurface(...);
static_assert(!decltype(HasAffineSurface<CutFemWeight_LS<double,long double>>(0))::value,"Removed API AffineSurface remains public.");

template <class T> auto HasAffineLineSliceSurface(int) -> decltype(std::declval<T&>().GetPhysicalSurfaceWeightsFromNodalValuesByLineSlices(std::declval<std::vector<double>&>(),std::declval<std::vector<std::vector<double>>&>(),std::declval<std::vector<double>&>(),std::declval<std::vector<double>&>()),std::true_type());
template <class T> std::false_type HasAffineLineSliceSurface(...);
static_assert(!decltype(HasAffineLineSliceSurface<CutFemWeight_LS<double,long double>>(0))::value,"Removed API AffineLineSliceSurface remains public.");

template <class T> auto HasGeometryLineSliceSurface(int) -> decltype(std::declval<T&>().GetPhysicalSurfaceWeightsFromNodalGeometryByLineSlices(std::declval<std::vector<double>&>(),std::declval<std::vector<std::vector<double>>&>(),std::declval<std::vector<double>&>(),std::declval<std::vector<double>&>()),std::true_type());
template <class T> std::false_type HasGeometryLineSliceSurface(...);
static_assert(!decltype(HasGeometryLineSliceSurface<CutFemWeight_LS<double,long double>>(0))::value,"Removed API GeometryLineSliceSurface remains public.");

template <class T> auto HasQuad9Direction(int) -> decltype(std::declval<T&>().SetQuad9IntegrationDirection(0U),std::true_type());
template <class T> std::false_type HasQuad9Direction(...);
static_assert(!decltype(HasQuad9Direction<CutFemWeight_LS<double,long double>>(0))::value,"Removed API Quad9Direction remains public.");

template <class T> auto HasQuad9Validation(int) -> decltype(std::declval<T&>().SetQuad9DirectionValidation(0U),std::true_type());
template <class T> std::false_type HasQuad9Validation(...);
static_assert(!decltype(HasQuad9Validation<CutFemWeight_LS<double,long double>>(0))::value,"Removed API Quad9Validation remains public.");

template <class T> auto HasTri7Direction(int) -> decltype(std::declval<T&>().SetTri7IntegrationDirection(0U),std::true_type());
template <class T> std::false_type HasTri7Direction(...);
static_assert(!decltype(HasTri7Direction<CutFemWeight_LS<double,long double>>(0))::value,"Removed API Tri7Direction remains public.");

template <class T> auto HasTri7Validation(int) -> decltype(std::declval<T&>().SetTri7DirectionValidation(0U),std::true_type());
template <class T> std::false_type HasTri7Validation(...);
static_assert(!decltype(HasTri7Validation<CutFemWeight_LS<double,long double>>(0))::value,"Removed API Tri7Validation remains public.");

template <class T> auto HasHex27Direction(int) -> decltype(std::declval<T&>().SetHex27IntegrationDirection(0U),std::true_type());
template <class T> std::false_type HasHex27Direction(...);
static_assert(!decltype(HasHex27Direction<CutFemWeight_LS<double,long double>>(0))::value,"Removed API Hex27Direction remains public.");

template <class T> auto HasHex27Validation(int) -> decltype(std::declval<T&>().SetHex27DirectionValidation(0U),std::true_type());
template <class T> std::false_type HasHex27Validation(...);
static_assert(!decltype(HasHex27Validation<CutFemWeight_LS<double,long double>>(0))::value,"Removed API Hex27Validation remains public.");

template <class T> auto HasTet15Direction(int) -> decltype(std::declval<T&>().SetTet15IntegrationDirection(0U),std::true_type());
template <class T> std::false_type HasTet15Direction(...);
static_assert(!decltype(HasTet15Direction<CutFemWeight_LS<double,long double>>(0))::value,"Removed API Tet15Direction remains public.");

template <class T> auto HasTet15Validation(int) -> decltype(std::declval<T&>().SetTet15DirectionValidation(0U),std::true_type());
template <class T> std::false_type HasTet15Validation(...);
static_assert(!decltype(HasTet15Validation<CutFemWeight_LS<double,long double>>(0))::value,"Removed API Tet15Validation remains public.");

template <class T> auto HasTetBase(int) -> decltype(std::declval<T&>().SetTetBaseType(1U),std::true_type());
template <class T> std::false_type HasTetBase(...);
static_assert(!decltype(HasTetBase<CutFemWeight_LS<double,long double>>(0))::value,"Removed API TetBase remains public.");

unsigned contractChecks = 0;

void Check(bool condition, const char* message) {
  ++contractChecks;
  if(!condition) throw std::runtime_error(message);
}

void SameWeights(const std::vector<double>& first, const std::vector<double>& second) {
  Check(first.size() == second.size(), "Cache output size changed.");
  for(unsigned i = 0; i < first.size(); ++i) Check(first[i] == second[i], "Cache output changed.");
}

void InvalidDirection(CutFemWeight_LS<double,long double>& rule, unsigned value) {
  bool rejected = false;
  try { rule.SetIntegrationDirection(value); }
  catch(const std::invalid_argument&) { rejected = true; }
  Check(rejected,"Invalid integration direction was accepted.");
}

void InvalidValidation(CutFemWeight_LS<double,long double>& rule, unsigned value) {
  bool rejected = false;
  try { rule.SetDirectionValidation(value); }
  catch(const std::invalid_argument&) { rejected = true; }
  Check(rejected,"Invalid validation mode was accepted.");
}

void InvalidNodalRequest(CutFemWeight_LS<double,long double>& rule, const std::vector<double>& nodal, double absoluteTolerance, double relativeTolerance) {
  std::vector<double> positive;
  std::vector<double> surface;
  bool rejected = false;
  try { rule.GetWeightsFromNodalValues(nodal,positive,surface,absoluteTolerance,relativeTolerance,true); }
  catch(const std::invalid_argument&) { rejected = true; }
  Check(rejected,"Invalid nodal request was accepted.");
}

void InvalidGeometry(CutFemWeight_LS<double,long double>& rule, const std::vector<double>& nodal, const std::vector<std::vector<double>>& coordinates) {
  std::vector<double> positive;
  std::vector<double> surface;
  bool rejected = false;
  try { rule.GetPhysicalSurfaceWeightsFromNodalGeometry(nodal,coordinates,positive,surface); }
  catch(const std::invalid_argument&) { rejected = true; }
  Check(rejected,"Invalid geometry was accepted.");
}

void CheckContract() {
  const GeomElType elements[6] = {LINE,QUAD,TRI,HEX,TET,WEDGE};
  const unsigned nodes[6] = {3,9,7,27,15,21};
  const unsigned maximumDirection[6] = {1,2,3,6,4,6};
  for(unsigned e = 0; e < 6; ++e) {
    CutFemWeight_LS<double,long double> rule(elements[e],0,"legendre");
    const CutFemWeight_LS<double,long double>& constant = rule;
    Check(constant.GetCutFEMQuadratureOrder() == 0,"Constructor degree mismatch.");
    Check(constant.GetGaussQuadratureOrder() == 0,"Gauss degree mismatch.");
    Check(constant.GetDimension() == (e == 0 ? 1U : (e <= 2 ? 2U : 3U)),"Dimension mismatch.");
    const unsigned count = constant.GetGaussQuadraturePointNumber();
    Check(count > 0,"No quadrature points.");
    for(unsigned i = 0; i < count; ++i) Check(std::isfinite(constant.GetGaussWeightPointer()[i]),"Invalid Gauss weight.");
    for(unsigned d = 0; d < constant.GetDimension(); ++d) {
      for(unsigned i = 0; i < count; ++i) Check(std::isfinite(constant.GetGaussCoordinatePointer(d)[i]),"Invalid Gauss coordinate.");
    }
    for(unsigned direction = 0; direction <= maximumDirection[e]; ++direction) rule.SetIntegrationDirection(direction);
    for(unsigned validation = 0; validation <= 2; ++validation) rule.SetDirectionValidation(validation);
    rule.SetIntegrationDirection(0);
    rule.SetDirectionValidation(0);
    const std::vector<double> nodal(nodes[e],1.);
    std::vector<double> positive;
    std::vector<double> surface;
    rule.GetWeightsFromNodalValues(nodal,positive,surface,1.e-8,1.e-8,true);
    InvalidDirection(rule,maximumDirection[e] + 1U);
    InvalidDirection(rule,std::numeric_limits<unsigned>::max());
    InvalidValidation(rule,3);
    InvalidValidation(rule,std::numeric_limits<unsigned>::max());
    std::vector<double> repeatedPositive;
    std::vector<double> repeatedSurface;
    rule.GetWeightsFromNodalValues(nodal,repeatedPositive,repeatedSurface,1.e-8,1.e-8,true);
    SameWeights(positive,repeatedPositive);
    SameWeights(surface,repeatedSurface);
    rule.ClearCache();
    rule.GetWeightsFromNodalValues(nodal,repeatedPositive,repeatedSurface,1.e-8,1.e-8,true);
    SameWeights(positive,repeatedPositive);
    SameWeights(surface,repeatedSurface);
    std::vector<double> invalidNodal = nodal;
    invalidNodal.pop_back();
    InvalidNodalRequest(rule,invalidNodal,1.e-8,1.e-8);
    invalidNodal = nodal;
    invalidNodal[0] = std::numeric_limits<double>::quiet_NaN();
    InvalidNodalRequest(rule,invalidNodal,1.e-8,1.e-8);
    invalidNodal[0] = std::numeric_limits<double>::infinity();
    InvalidNodalRequest(rule,invalidNodal,1.e-8,1.e-8);
    InvalidNodalRequest(rule,nodal,0.,1.e-8);
    InvalidNodalRequest(rule,nodal,-1.,1.e-8);
    InvalidNodalRequest(rule,nodal,1.e-8,-1.);
    InvalidNodalRequest(rule,nodal,std::numeric_limits<double>::infinity(),1.e-8);
    InvalidNodalRequest(rule,nodal,1.e-8,std::numeric_limits<double>::quiet_NaN());
    bool rejected = false;
    try { rule.GetWeightsFromNodalValues(nodal,repeatedPositive,repeatedPositive); }
    catch(const std::invalid_argument&) { rejected = true; }
    Check(rejected,"Aliased output vectors were accepted.");
    std::vector<std::vector<double>> coordinates(constant.GetDimension(),std::vector<double>(nodes[e],0.));
    coordinates[0].pop_back();
    InvalidGeometry(rule,nodal,coordinates);
    coordinates[0].push_back(std::numeric_limits<double>::quiet_NaN());
    InvalidGeometry(rule,nodal,coordinates);
    coordinates.pop_back();
    InvalidGeometry(rule,nodal,coordinates);
  }
  bool rejected = false;
  try { CutFemWeight_LS<double,long double> rule(static_cast<GeomElType>(99),0,"legendre"); }
  catch(const std::invalid_argument&) { rejected = true; }
  Check(rejected,"An invalid geometry was accepted.");
}

int main() {
  try {
    CheckContract();
    std::cout << "Passed " << contractChecks << " input, setter-range, const-getter and cache checks; 29 removed API signatures are inaccessible.\n";
  }
  catch(const std::exception& error) {
    std::cerr << error.what() << '\n';
    return 1;
  }
  return 0;
}
