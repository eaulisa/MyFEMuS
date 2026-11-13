#ifndef FEMUS_PLANAR_ENDPOINT_INTEGRATION_HPP
#define FEMUS_PLANAR_ENDPOINT_INTEGRATION_HPP

#include "EndpointStripQuadrature.hpp"

namespace FELevelSetMoments {

class PlanarEndpointIntegration {
public:
  // QUAD and full TRI share the unchanged LINE graph kernels and exact endpoint quadrature.
  static void Range(const std::vector<std::vector<long double>>& coefficient, bool triangle, long double& minimum, long double& maximum) {
    minimum = std::numeric_limits<long double>::infinity();
    maximum = -minimum;
    if(!triangle) {
      long double bernstein[3][3];
      LineSliceUtilities::CanonicalBernstein(coefficient,bernstein);
      for(unsigned i = 0; i < 3; ++i) for(unsigned j = 0; j < 3; ++j) { minimum = std::min(minimum,bernstein[i][j]); maximum = std::max(maximum,bernstein[i][j]); }
      return;
    }
    unsigned degree = 0;
    for(unsigned i = 0; i < 3; ++i) for(unsigned j = 0; j < 3; ++j) if(coefficient[i][j] != 0) degree = std::max(degree,i+j);
    for(unsigned i = 0; i <= degree; ++i) for(unsigned j = 0; j <= degree-i; ++j) {
      long double value = 0;
      for(unsigned a = 0; a <= i && a < 3; ++a) for(unsigned b = 0; b <= j && b < 3; ++b) {
        long double multiplier = 1;
        for(unsigned n = 0; n < a; ++n) multiplier *= static_cast<long double>(i-n)/(degree-n);
        for(unsigned n = 0; n < b; ++n) multiplier *= static_cast<long double>(j-n)/(degree-a-n);
        value += multiplier*coefficient[a][b];
      }
      minimum = std::min(minimum,value);
      maximum = std::max(maximum,value);
    }
  }

  static unsigned GraphOrder(const std::vector<std::vector<long double>>& coefficient, bool triangle, unsigned requested, bool& strict) {
    unsigned selected = 0;
    long double best = 0;
    for(unsigned order = 1; order <= 2; ++order) {
      if(requested != 0 && order != requested) continue;
      const unsigned axis = order == 1 ? 1U : 0U;
      std::vector<std::vector<long double>> derivative(3,std::vector<long double>(3,0));
      for(unsigned i = 0; i < 3; ++i) for(unsigned j = 0; j < 3; ++j) if((axis == 0 ? i : j) != 0) derivative[i-(axis == 0)][j-(axis == 1)] += (axis == 0 ? i : j)*coefficient[i][j];
      long double minimum;
      long double maximum;
      Range(derivative,triangle,minimum,maximum);
      const long double magnitude = std::max(std::abs(minimum),std::abs(maximum));
      const long double padding = 4096*std::numeric_limits<long double>::epsilon()*magnitude;
      if(!((minimum >= -padding && maximum > padding) || (maximum <= padding && minimum < -padding))) continue;
      const long double lower = minimum > padding ? minimum : (maximum < -padding ? -maximum : 0);
      const long double score = lower > 0 ? 1+lower/magnitude : magnitude/(1+magnitude);
      if(score > best && !LineIntegration::HasVerticalFiber(coefficient,order)) { selected = order; strict = lower > 0; best = score; }
    }
    return selected;
  }

  static void Integrate(const std::vector<std::vector<long double>>& coefficient, unsigned degree, bool triangle, bool computeDelta, long double absoluteH, long double absoluteD, long double relative, unsigned maximumDepth, unsigned maximumEvaluations, std::vector<long double>& values, std::vector<long double>& errors, unsigned& evaluations, unsigned& accepted, unsigned& selected, const std::vector<std::vector<std::vector<long double>>>* gradient = nullptr, const std::array<long double,5>* coordinates = nullptr, const std::vector<unsigned char>* requested = nullptr, unsigned direction = 0, bool fullBoundaryDelta = false, const std::vector<unsigned>* zeroBoundary = nullptr, const OuterQuadratureRules& rules = OuterQuadratureRules(), bool inverseSurface = false) {
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
      std::vector<std::vector<long double>> centered;
      ReferencePolynomialAffine2D(coefficient,0.5L,0.5L,0.5L,0.5L,centered);
      for(unsigned axis = 0; axis < 2; ++axis) if(std::fma(0.5L,(*coordinates)[2U*axis+1U],(*coordinates)[2U*axis]) == 0) {
        bool even = true;
        for(unsigned i = 0; i < 3; ++i) for(unsigned j = 0; j < 3; ++j) if((axis == 0 ? i : j)%2U != 0 && centered[i][j] != 0) even = false;
        if(!even) continue;
        bool metric = true;
        if(gradient != nullptr) for(unsigned component = 0; component < gradient->size(); ++component) {
          std::vector<std::vector<long double>> field;
          ReferencePolynomialAffine2D((*gradient)[component],0.5L,0.5L,0.5L,0.5L,field);
          bool componentEven = true;
          bool componentOdd = true;
          for(unsigned i = 0; i < field.size(); ++i) for(unsigned j = 0; j < field[i].size(); ++j) if(field[i][j] != 0) { if((axis == 0 ? i : j)%2U != 0) componentEven = false; else componentOdd = false; }
          metric = metric && (componentEven || componentOdd);
        }
        for(unsigned i = 0; i <= degree; ++i) for(unsigned j = 0; j <= degree-i; ++j) if((axis == 0 ? i : j)%2U != 0) {
          const unsigned n = ReferenceMomentIndex2D(degree,i,j);
          needed[n] = 0;
          if(metric) needed[count+n] = 0;
        }
      }
    }
    requested = &needed;
    values.assign(2U*count,0);
    errors.assign(2U*count,0);
    long double minimum;
    long double maximum;
    Range(coefficient,triangle,minimum,maximum);
    const LineIntegrationStrip full = {0,1,0,0,1,triangle ? -1.0L : 0.0L};
    if(minimum == 0 && maximum == 0) {
      if(computeDelta) throw std::domain_error("The delta of an identically zero planar field is undefined.");
      LineIntegration::FullStripMoments(degree,full,0,1,values,errors,coordinates,requested);
      for(unsigned n = 0; n < count; ++n) { values[n] *= 0.5L; errors[n] *= 0.5L; }
      selected = 0;
      ++accepted;
      return;
    }
    bool zeroTrace = false;
    if(computeDelta && zeroBoundary != nullptr) for(unsigned n = 0; n < zeroBoundary->size(); ++n) zeroTrace = zeroTrace || (*zeroBoundary)[n] != 0;
    const long double padding = 4096*std::numeric_limits<long double>::epsilon()*std::max(std::abs(minimum),std::abs(maximum));
    if(!zeroTrace && (minimum > padding || maximum < -padding)) {
      if(minimum > 0) LineIntegration::FullStripMoments(degree,full,0,1,values,errors,coordinates,requested);
      selected = 0;
      ++accepted;
      return;
    }
    bool strict = false;
    selected = GraphOrder(coefficient,triangle,direction,strict);
    if(selected == 0) {
      if(!EndpointStripQuadrature::Try(coefficient,degree,triangle,computeDelta,absoluteH,absoluteD,relative,maximumDepth,maximumEvaluations,values,errors,evaluations,accepted,selected,gradient,1,nullptr,nullptr,coordinates,requested,direction,fullBoundaryDelta,rules,inverseSurface)) throw std::runtime_error("Unsupported endpoint topology: "+EndpointStripQuadrature::Statistics().failure);
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
      LineIntegration::IntegrateStrip(working,degree,strip,true,!strict,computeDelta,budgetFactor*absoluteH,budgetFactor*absoluteD,budgetFactor*relative,maximumDepth,maximumEvaluations,evaluations,accepted,passValues,passErrors,gradient,&parent,&failed,rules,inverseSurface);
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
  static unsigned Compare(const std::vector<std::vector<long double>>& coefficient, unsigned degree, bool triangle, bool computeDelta, long double absoluteH, long double absoluteD, long double relative, unsigned maximumDepth, unsigned maximumEvaluations, const std::vector<long double>& values, std::vector<long double>& errors, unsigned& evaluations, unsigned& accepted, unsigned selected, const std::vector<std::vector<std::vector<long double>>>* gradient, const std::array<long double,5>* coordinates, bool fullBoundaryDelta = false, const std::vector<unsigned>* zeroBoundary = nullptr, const OuterQuadratureRules& rules = OuterQuadratureRules(), bool inverseSurface = false) {
    QuadratureWarningScope warningScope("Independent planar validation");
    bool depends[2] = {false,false};
    for(unsigned i = 0; i < 3; ++i) for(unsigned j = 0; j < 3; ++j) if(coefficient[i][j] != 0) { depends[0] = depends[0] || i != 0; depends[1] = depends[1] || j != 0; }
    if(!triangle && depends[0] != depends[1]) {
      bool constantMetric = true;
      if(gradient != nullptr) for(unsigned component = 0; component < gradient->size(); ++component) for(unsigned i = 0; i < (*gradient)[component].size(); ++i) for(unsigned j = 0; j < (*gradient)[component][i].size(); ++j) if(i+j != 0 && (*gradient)[component][i][j] != 0) constantMetric = false;
      if(constantMetric) {
        const unsigned axis = depends[0] ? 0U : 1U;
        const std::array<long double,5> parent = coordinates == nullptr ? std::array<long double,5>{0,1,0,1,1} : *coordinates;
        const long double origin = parent[2U*axis];
        const long double scale = parent[2U*axis+1U];
        const long double a = (axis == 0 ? coefficient[2][0] : coefficient[0][2])/(scale*scale);
        const long double b = (axis == 0 ? coefficient[1][0] : coefficient[0][1])/scale-2*a*origin;
        const long double c = coefficient[0][0]-b*origin-a*origin*origin;
        std::vector<long double> positive;
        std::vector<long double> delta;
        LineIntegration::IntegratePolynomial(a,b,c,degree,std::min(origin,origin+scale),std::max(origin,origin+scale),computeDelta,fullBoundaryDelta,positive,delta);
        const long double metric = gradient == nullptr ? 1 : SurfaceGradientFactor(*gradient,0,0,inverseSurface);
        const unsigned count = values.size()/2U;
        for(unsigned i = 0; i <= degree; ++i) for(unsigned j = 0; j <= degree-i; ++j) {
          const unsigned n = ReferenceMomentIndex2D(degree,i,j);
          const unsigned power = axis == 0 ? i : j;
          const unsigned transverse = axis == 0 ? j : i;
          const long double otherOrigin = parent[2U*(1U-axis)];
          const long double otherScale = parent[2U*(1U-axis)+1U];
          const long double factor = parent[4]*MonomialIntegral(std::min(otherOrigin,otherOrigin+otherScale),std::max(otherOrigin,otherOrigin+otherScale),transverse)/std::abs(scale*otherScale);
          const long double independentH = factor*positive[power];
          const long double independentD = factor*delta[power]*metric;
          errors[n] += std::abs(values[n]-independentH)+64*std::numeric_limits<long double>::epsilon()*std::abs(independentH);
          errors[count+n] += std::abs(values[count+n]-independentD)+64*std::numeric_limits<long double>::epsilon()*std::abs(independentD);
          if(errors[n] > absoluteH+relative*std::abs(values[n])) WarnQuadrature("Independent LINE validation","volume moment disagreement",n,errors[n],absoluteH+relative*std::abs(values[n]));
          if(errors[count+n] > absoluteD+relative*std::abs(values[count+n])) WarnQuadrature("Independent LINE validation","interface moment disagreement",count+n,errors[count+n],absoluteD+relative*std::abs(values[count+n]));
        }
        return 1;
      }
    }
    std::vector<long double> second;
    std::vector<long double> secondErrors;
    unsigned independent = 0;
    try { Integrate(coefficient,degree,triangle,computeDelta,0.1L*absoluteH,0.1L*absoluteD,0.1L*relative,maximumDepth,maximumEvaluations,second,secondErrors,evaluations,accepted,independent,gradient,coordinates,nullptr,selected == 0 ? 0 : 3U-selected,fullBoundaryDelta,zeroBoundary,rules,inverseSurface); }
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

#endif
