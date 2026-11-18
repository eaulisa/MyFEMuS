#ifndef FEMUS_EX47_PREPARED_PLANAR_ENDPOINT_INTEGRATION_HPP
#define FEMUS_EX47_PREPARED_PLANAR_ENDPOINT_INTEGRATION_HPP

#include "../../../../../src/quadrature/cutfem_ls/PlanarEndpointIntegration.hpp"
#include "EndpointStripQuadrature.hpp"

namespace FELevelSetMoments {
namespace Ex47Prepared {

class PlanarEndpointIntegration {
public:
  // QUAD and full TRI share prepared LINE graph kernels and exact endpoint quadrature.
  static unsigned GraphOrder(const std::vector<std::vector<long double>>& coefficient, bool triangle, unsigned requested, bool& strict, const PreparedPlane* prepared = nullptr) {
    if(prepared == nullptr) throw std::logic_error("Prepared planar integration requires a level-set record.");
    unsigned selected = 0;
    long double best = -1;
    for(unsigned order = 1; order <= 2; ++order) {
      if(requested != 0 && requested != order) continue;
      const CutFemLevelSetData::GraphCertificate& certificate = prepared->data->GetGraphCertificate(PreparedPlane::InnerAxis(order));
      if(certificate.sign == 0 || (!certificate.strict && certificate.zeroFiber != -1)) continue;
      const long double magnitude = std::max(std::abs(certificate.bounds[0]),std::abs(certificate.bounds[1]));
      const long double lower = certificate.bounds[0] > 0 ? certificate.bounds[0] : (certificate.bounds[1] < 0 ? -certificate.bounds[1] : 0);
      const long double score = certificate.strict ? 1+lower/magnitude : magnitude/(1+magnitude);
      if(score > best) { best = score; selected = order; strict = certificate.strict; }
    }
    return selected == 0 ? ::FELevelSetMoments::PlanarEndpointIntegration::GraphOrder(coefficient,triangle,requested,strict) : selected;
  }

  static void Integrate(const std::vector<std::vector<long double>>& coefficient, unsigned degree, bool triangle, bool computeDelta, long double absoluteH, long double absoluteD, long double relative, unsigned maximumDepth, unsigned maximumEvaluations, std::vector<long double>& values, std::vector<long double>& errors, unsigned& evaluations, unsigned& accepted, unsigned& selected, const std::vector<std::vector<std::vector<long double>>>* gradient = nullptr, const std::array<long double,5>* coordinates = nullptr, const std::vector<unsigned char>* requested = nullptr, unsigned direction = 0, bool fullBoundaryDelta = false, const std::vector<unsigned>* zeroBoundary = nullptr, const OuterQuadratureRules& rules = OuterQuadratureRules(), bool inverseSurface = false, int normalComponent = -1, const PreparedPlane* prepared = nullptr) {
    QuadratureWarningScope warningScope(triangle ? "TRI endpoint integration" : "QUAD endpoint integration");
    rules.Validate();
    if(coefficient.size() != 3 || degree > 64 || maximumEvaluations < 3 || direction > 2 || !(absoluteH > 0) || !(absoluteD > 0) || !(relative >= 0)) throw std::invalid_argument("Invalid endpoint integration input.");
    for(unsigned i = 0; i < 3; ++i) if(coefficient[i].size() != 3) throw std::invalid_argument("Endpoint integration requires the exact 3 by 3 parent slice.");
    bool depends[2] = {false,false};
    for(unsigned i = 0; i < 3; ++i) for(unsigned j = 0; j < 3; ++j) if(coefficient[i][j] != 0) { depends[0] = depends[0] || i != 0; depends[1] = depends[1] || j != 0; }
    if(depends[0] != depends[1]) direction = depends[0] ? 2U : 1U;
    const unsigned count = (degree+1U)*(degree+2U)/2U;
    if(requested != nullptr && requested->size() != 2U*count) throw std::invalid_argument("Invalid selective moment mask.");
    std::vector<unsigned char> needed = requested == nullptr ? std::vector<unsigned char>(2U*count,1) : *requested;
    if(!triangle && coordinates != nullptr) {
      for(unsigned axis = 0; axis < 2; ++axis) if(std::fma(0.5L,(*coordinates)[2U*axis+1],(*coordinates)[2U*axis]) == 0 && prepared != nullptr && prepared->data->GetReflectionParity(axis) == 1 && gradient == nullptr) {
        for(unsigned i = 0; i <= degree; ++i) for(unsigned j = 0; j <= degree-i; ++j) if((axis == 0 ? i : j)%2U != 0) { const unsigned n = ReferenceMomentIndex2D(degree,i,j); needed[n] = 0; needed[count+n] = 0; }
      }
    }
    requested = &needed;
    values.assign(2U*count,0);
    errors.assign(2U*count,0);
    long double minimum;
    long double maximum;
    if(prepared == nullptr) throw std::logic_error("Missing prepared planar data.");
    const std::array<long double,2> bounds = prepared->data->GetBounds();
    minimum = bounds[0]/prepared->normalization;
    maximum = bounds[1]/prepared->normalization;
    const LineIntegrationStrip full = {0.0L,1,0.0L,0,1,triangle ? -1.0L : 0.0L};
    if(minimum == 0 && maximum == 0) {
      if(computeDelta) throw std::domain_error("The delta of an identically zero planar field is undefined.");
      LineIntegration::FullStripMoments(degree,full,full.left,full.right,values,errors,coordinates,requested);
      for(unsigned n = 0; n < count; ++n) { values[n] *= 0.5L; errors[n] *= 0.5L; }
      selected = 0;
      ++accepted;
      return;
    }
    bool zeroTrace = false;
    if(computeDelta && zeroBoundary != nullptr) for(unsigned n = 0; n < zeroBoundary->size(); ++n) zeroTrace = zeroTrace || (*zeroBoundary)[n] != 0;
    const long double padding = 4096*std::numeric_limits<long double>::epsilon()*std::max(std::abs(minimum),std::abs(maximum));
    if(!zeroTrace && (prepared->data->GetPhase() != 0 || minimum > padding || maximum < -padding)) {
      if(prepared->data->GetPhase() > 0 || minimum > 0) LineIntegration::FullStripMoments(degree,full,full.left,full.right,values,errors,coordinates,requested);
      selected = 0;
      ++accepted;
      return;
    }
    bool strict = false;
    selected = GraphOrder(coefficient,triangle,direction,strict,prepared);
    if(selected == 0) {
      if(!EndpointStripQuadrature::Try(coefficient,degree,triangle,computeDelta,absoluteH,absoluteD,relative,maximumDepth,maximumEvaluations,values,errors,evaluations,accepted,selected,gradient,1,nullptr,nullptr,coordinates,requested,direction,fullBoundaryDelta,rules,inverseSurface,normalComponent,prepared)) throw std::runtime_error("Unsupported endpoint topology: "+EndpointStripQuadrature::Statistics().failure);
      return;
    }
    std::vector<std::vector<long double>> working = coefficient;
    std::vector<std::vector<std::vector<long double>>> transposed;
    if(selected == 2) {
      for(unsigned i = 0; i < 3; ++i) for(unsigned j = i+1U; j < 3; ++j) std::swap(working[i][j],working[j][i]);
      if(gradient != nullptr) { SurfaceGradientTranspose(*gradient,transposed); gradient = &transposed; }
    }
    const std::array<long double,5> parent = coordinates == nullptr ? std::array<long double,5>{0,1,0,1,1} : LineIntegration::MomentCoordinates(coordinates,selected);
    const std::vector<unsigned char> active = LineIntegration::MomentRequest(degree,requested,selected,computeDelta);
    LineIntegrationStrip strip = full;
    strip.fullBoundaryDelta = fullBoundaryDelta;
    if(zeroBoundary != nullptr) {
      const unsigned axis = selected == 1 ? 1U : 0U;
      strip.certifiedLowerZero = (*zeroBoundary)[triangle ? axis : 2U*axis] != 0;
      strip.certifiedUpperZero = (*zeroBoundary)[triangle ? 2U : 2U*axis+1U] != 0;
    }
    std::vector<unsigned char> failed = active;
    std::vector<long double> best(2U*count,std::numeric_limits<long double>::quiet_NaN());
    std::vector<long double> bestErrors(2U*count,std::numeric_limits<long double>::infinity());
    long double budgetFactor = 1;
    for(unsigned pass = 0; pass < 5; ++pass) {
      std::vector<long double> passValues(2U*count,0);
      std::vector<long double> passErrors(2U*count,0);
      LineIntegration::IntegrateStrip(working,degree,strip,true,!strict,computeDelta,budgetFactor*absoluteH,budgetFactor*absoluteD,budgetFactor*relative,maximumDepth,maximumEvaluations,evaluations,accepted,passValues,passErrors,gradient,&parent,&failed,rules,inverseSurface,normalComponent,prepared,selected);
      KeepBestQuadrature(passValues,passErrors,failed,best,bestErrors);
      long double worst = 0;
      for(unsigned n = 0; n < failed.size(); ++n) if(failed[n]) {
        const long double tolerance = (n < count ? absoluteH : absoluteD)+relative*std::abs(passValues[n]);
        if(std::isfinite(passValues[n]) && std::isfinite(passErrors[n]) && passErrors[n] <= tolerance) { values[n] = passValues[n]; errors[n] = passErrors[n]; failed[n] = 0; }
        else worst = std::max(worst,passErrors[n]/tolerance);
      }
      if(std::find(failed.begin(),failed.end(),1) == failed.end()) break;
      if(pass == 4 || !std::isfinite(worst) || evaluations >= maximumEvaluations) {
        RetainUncertifiedQuadrature("Planar moments",evaluations >= maximumEvaluations ? "evaluation budget exhausted before accumulated tolerance" : "accumulated error exceeds tolerance after selective retries",best,bestErrors,absoluteH,absoluteD,relative,failed,values,errors);
        break;
      }
      budgetFactor /= std::max(2.0L,2*worst);
    }
    if(selected == 2) for(unsigned i = 0; i <= degree; ++i) for(unsigned j = i+1U; j <= degree-i; ++j) for(unsigned pair = 0; pair < 2; ++pair) {
      std::swap(values[pair*count+ReferenceMomentIndex2D(degree,i,j)],values[pair*count+ReferenceMomentIndex2D(degree,j,i)]);
      std::swap(errors[pair*count+ReferenceMomentIndex2D(degree,i,j)],errors[pair*count+ReferenceMomentIndex2D(degree,j,i)]);
    }
    CheckQuadratureTolerance("Planar moments","accumulated error exceeds requested tolerance",values,errors,absoluteH,absoluteD,relative,requested);
  }
  static unsigned Compare(const std::vector<std::vector<long double>>& coefficient, unsigned degree, bool triangle, bool computeDelta, long double absoluteH, long double absoluteD, long double relative, unsigned maximumDepth, unsigned maximumEvaluations, const std::vector<long double>& values, std::vector<long double>& errors, unsigned& evaluations, unsigned& accepted, unsigned selected, const std::vector<std::vector<std::vector<long double>>>* gradient, const std::array<long double,5>* coordinates, bool fullBoundaryDelta = false, const std::vector<unsigned>* zeroBoundary = nullptr, const OuterQuadratureRules& rules = OuterQuadratureRules(), bool inverseSurface = false, int normalComponent = -1, const PreparedPlane* prepared = nullptr) {
    bool depends[2] = {false,false}, constantMetric = true;
    for(unsigned i = 0; i < 3; ++i) for(unsigned j = 0; j < 3; ++j) if(coefficient[i][j] != 0) { depends[0] = depends[0] || i != 0; depends[1] = depends[1] || j != 0; }
    if(gradient != nullptr) for(const std::vector<std::vector<long double>>& component : *gradient) for(unsigned i = 0; i < component.size(); ++i) for(unsigned j = 0; j < component[i].size(); ++j) if(i+j != 0 && component[i][j] != 0) constantMetric = false;
    if(!triangle && depends[0] != depends[1] && constantMetric) return ::FELevelSetMoments::PlanarEndpointIntegration::Compare(coefficient,degree,triangle,computeDelta,absoluteH,absoluteD,relative,maximumDepth,maximumEvaluations,values,errors,evaluations,accepted,selected,gradient,coordinates,fullBoundaryDelta,zeroBoundary,rules,inverseSurface,normalComponent);

    QuadratureWarningScope warningScope("Independent planar validation");
    std::vector<long double> second;
    std::vector<long double> secondErrors;
    unsigned independent = 0;
    try { Integrate(coefficient,degree,triangle,computeDelta,0.1L*absoluteH,0.1L*absoluteD,0.1L*relative,maximumDepth,maximumEvaluations,second,secondErrors,evaluations,accepted,independent,gradient,coordinates,nullptr,selected == 0 ? 0 : 3U-selected,fullBoundaryDelta,zeroBoundary,rules,inverseSurface,normalComponent,prepared); }
    catch(const std::runtime_error& error) { WarnQuadrature("Independent planar validation",std::string("validation unavailable: ")+error.what(),0,std::numeric_limits<long double>::infinity(),absoluteD); return 0; }
    const unsigned count = values.size()/2U;
    for(unsigned n = 0; n < values.size(); ++n) {
      errors[n] += secondErrors[n]+std::abs(values[n]-second[n]);
      if(errors[n] > (n < count ? absoluteH : absoluteD)+relative*std::abs(values[n])) WarnQuadrature("Independent planar validation","endpoint projections disagree",n,errors[n],(n < count ? absoluteH : absoluteD)+relative*std::abs(values[n]));
    }
    return independent != selected ? 1U : 0U;
  }

};

}

}

#endif
