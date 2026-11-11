#include "CutFemWeight_LS.hpp"
#include "Basis.hpp"
#include <cmath>
#include <iostream>
#include <sstream>

unsigned ruleSelectionChecks = 0;

void RuleSelectionCheck(bool condition, const char* message) {
  ++ruleSelectionChecks;
  if(!condition) throw std::runtime_error(message);
}

struct RuleSelectionContext {
  unsigned evaluations = 0;
  unsigned mode = 0;
  long double endpoint = 0;
  long double delta = 1;
  long double cutoff = 0;
};

void RuleSelectionSample(RuleSelectionContext& context, const FELevelSetMoments::OuterIntervalPoint& point, const std::vector<unsigned char>& active, std::vector<long double>& values, std::vector<long double>& errors, long double absoluteH, long double absoluteD, long double relative) {
  ++context.evaluations;
  if(context.mode == 3 && point.parameter < context.cutoff) throw FELevelSetMoments::OuterQuadratureResolutionError("An inner slice reached coordinate resolution.");
  if(context.mode == 4 && point.parameter < context.cutoff) throw std::runtime_error("An unrelated integration failure.");
  if(context.mode != 0) RuleSelectionCheck(point.coordinate != context.endpoint,"A collapsed point reached the integrand.");
  values.assign(active.size(),0);
  errors.assign(active.size(),0);
  for(unsigned n = 0; n < active.size(); ++n) if(active[n]) {
    const long double coordinate = context.mode == 2 ? point.offset/context.delta : point.coordinate;
    values[n] = point.jacobian*FELevelSetMoments::ReferencePower(coordinate,n%11U);
  }
}

void RuleSelectionTables() {
  const unsigned orders[15] = {3,5,7,9,11,13,15,17,19,21,23,25,27,29,31};
  for(unsigned index = 0; index < 15; ++index) {
    const FELevelSetMoments::GaussKronrodRule& rule = FELevelSetMoments::GaussKronrodRules::Get(orders[index]);
    RuleSelectionCheck(rule.index == index && rule.kronrod == orders[index] && rule.gauss == (orders[index]-1U)/2U,"Incorrect cached quadrature descriptor.");
    RuleSelectionCheck(&rule == &FELevelSetMoments::GaussKronrodRules::Get(orders[index]),"Quadrature tables were not reused.");
    RuleSelectionCheck(FELevelSetMoments::GaussKronrodRules::Parse(rule.name) == orders[index],"A supported pair was not parsed.");
    const unsigned exactDegree = 3U*rule.gauss+1U+rule.gauss%2U;
    for(unsigned degree = 0; degree <= 10; ++degree) {
      long double high = 0;
      long double low = 0;
      for(unsigned q = 0; q < rule.size; ++q) for(unsigned side = 0; side < (q == 0 ? 1U : 2U); ++side) {
        const long double value = FELevelSetMoments::ReferencePower(side == 0 ? rule.nodes[q] : -rule.nodes[q],degree);
        high += rule.weights[q]*value;
        if(q%2U == (rule.gauss%2U == 0 ? 1U : 0U)) low += rule.gaussWeights[q/2U]*value;
      }
      const long double exact = degree%2U == 0 ? 2.0L/(degree+1U) : 0;
      if(degree <= exactDegree) RuleSelectionCheck(std::abs(high-exact) < 2e-17L,"A cached Kronrod rule failed polynomial exactness.");
      if(degree <= 2U*rule.gauss-1U) RuleSelectionCheck(std::abs(low-exact) < 2e-17L,"An embedded Gauss rule failed polynomial exactness.");
    }
  }
  bool rejected = false;
  try { FELevelSetMoments::GaussKronrodRules::Get(32); } catch(const std::invalid_argument&) { rejected = true; }
  RuleSelectionCheck(rejected,"An unsupported rule order was accepted.");
  rejected = false;
  try { FELevelSetMoments::GaussKronrodRules::Parse("7/31"); } catch(const std::invalid_argument&) { rejected = true; }
  RuleSelectionCheck(rejected,"A mismatched Gauss/Kronrod pair was accepted.");
  FELevelSetMoments::OuterQuadratureRules parsed;
  char program[] = "ex44";
  char unrelated[] = "--dim";
  char two[] = "2";
  char regular[] = "--gk-regular-2d";
  char pair[] = "1/3";
  char tangent[] = "--gk-tangent-3d";
  char largest[] = "15-31";
  char* arguments[] = {program,unrelated,two,regular,pair,tangent,largest};
  parsed = FELevelSetMoments::OuterQuadratureRules::ParseArguments(7,arguments);
  RuleSelectionCheck(parsed.regular2D == 3 && parsed.tangent3D == 31 && parsed.tangent2D == 15 && parsed.regular3D == 15,"Application rule options were not independent.");
}

void RuleSelectionIntervals(unsigned dimension, unsigned selected, bool mapped) {
  FELevelSetMoments::OuterQuadratureRules rules;
  if(dimension == 2) { rules.regular2D = selected; rules.tangent2D = selected; }
  else { rules.regular3D = selected; rules.tangent3D = selected; }
  const FELevelSetMoments::OuterInterval interval = mapped ? FELevelSetMoments::OuterInterval{0,1,0,1,1,2,true} : FELevelSetMoments::OuterInterval{0,1,0,0,0,1,true};
  for(unsigned accuracy = 0; accuracy < 2; ++accuracy) {
    FELevelSetMoments::OuterIntervalQuadrature::Statistics(dimension) = FELevelSetMoments::OuterIntervalStatistics();
    RuleSelectionContext context;
    std::vector<unsigned char> active(22,accuracy == 0 ? 0 : 1);
    active[0] = active[11] = 1;
    std::vector<long double> values(22,0);
    std::vector<long double> errors(22,0);
    unsigned accepted = 0;
    const long double tolerance = accuracy == 0 ? 1e-12L : (selected == 3 ? 1e-5L : (selected == 7 ? 1e-8L : 1e-12L));
    FELevelSetMoments::OuterIntervalQuadrature::Adaptive(context,RuleSelectionSample,interval,0,1,tolerance,tolerance,tolerance,0,32,200000,context.evaluations,accepted,dimension,active,values,errors,rules);
    const FELevelSetMoments::OuterIntervalStatistics& statistics = FELevelSetMoments::OuterIntervalQuadrature::Statistics(dimension);
    const unsigned index = FELevelSetMoments::GaussKronrodRules::Index(selected);
    for(unsigned rule = 0; rule < 15; ++rule) {
      const unsigned long long panels = mapped ? statistics.transformedRulePanels[rule] : statistics.regularRulePanels[rule];
      RuleSelectionCheck(rule == index ? panels > 0 : panels == 0,"An ordinary interval used an unselected rule.");
      RuleSelectionCheck(statistics.emergencyRulePanels[rule] == 0,"A rule was reduced without coordinate collapse.");
    }
    if(accuracy == 0) RuleSelectionCheck(context.evaluations == selected,"The primary panel used the wrong node count.");
    for(unsigned n = 0; n < active.size(); ++n) if(active[n]) {
      const long double exact = 1.0L/(n%11U+1U);
      RuleSelectionCheck(std::abs(values[n]-exact) <= tolerance*(1+std::abs(exact)),"A selected rule failed its adaptive moment tolerance.");
    }
    if(accuracy == 1 && selected <= 11) RuleSelectionCheck(statistics.subdivisions > 0,"A lower rule did not adapt its unresolved moments.");
  }
}

void RuleSelectionFallback(unsigned dimension, unsigned selected, unsigned ulps, unsigned expected, bool inaccurate = false) {
  FELevelSetMoments::OuterQuadratureRules rules;
  rules.tangent2D = rules.tangent3D = selected;
  FELevelSetMoments::OuterIntervalQuadrature::Statistics(dimension) = FELevelSetMoments::OuterIntervalStatistics();
  const long double delta = ulps*std::numeric_limits<long double>::epsilon();
  const FELevelSetMoments::OuterInterval interval = {1,1+delta,1,delta,1,2,true};
  RuleSelectionContext context;
  context.endpoint = 1;
  context.delta = delta;
  context.mode = inaccurate ? 2U : 1U;
  std::vector<unsigned char> active(22,0);
  active[inaccurate ? 10 : 0] = active[inaccurate ? 21 : 11] = 1;
  std::vector<long double> values(22,0);
  std::vector<long double> errors(22,0);
  unsigned accepted = 0;
  bool rejected = false;
  try { FELevelSetMoments::OuterIntervalQuadrature::Adaptive(context,RuleSelectionSample,interval,0,1,delta*1e-12L,delta*1e-12L,1e-12L,0,24,10000,context.evaluations,accepted,dimension,active,values,errors,rules); }
  catch(const std::runtime_error&) { rejected = true; }
  const FELevelSetMoments::OuterIntervalStatistics& statistics = FELevelSetMoments::OuterIntervalQuadrature::Statistics(dimension);
  RuleSelectionCheck(rejected == (expected == 0),"A failsafe result did not match its accuracy or coordinate limit.");
  RuleSelectionCheck(statistics.subdivisions == 0,"A collapsed interval was repeatedly bisected.");
  const unsigned primaryIndex = FELevelSetMoments::GaussKronrodRules::Index(selected);
  for(unsigned rule = primaryIndex; rule < 15; ++rule) RuleSelectionCheck(statistics.emergencyRulePanels[rule] == 0,"The failsafe raised the selected rule.");
  if(expected != 0) {
    RuleSelectionCheck(statistics.emergencyRulePanels[FELevelSetMoments::GaussKronrodRules::Index(expected)] == 1,"The failsafe skipped its required lower rule.");
    RuleSelectionCheck(values[0] > 0 && std::abs(values[0]-delta) < 1e-12L*delta,"A tiny positive contribution was discarded.");
  }
  else RuleSelectionCheck(statistics.resolutionLimits == 1,"An unresolved contribution was not reported explicitly.");
  if(inaccurate && selected == 31) RuleSelectionCheck(statistics.emergency15 > 0 && statistics.emergency7 > 0 && statistics.emergency3 > 0,"The inaccurate failsafe did not descend the complete ladder.");
}

void RuleSelectionConfiguredFallback(unsigned dimension, unsigned selected) {
  const FELevelSetMoments::GaussKronrodRule& primary = FELevelSetMoments::GaussKronrodRules::Get(selected);
  const FELevelSetMoments::GaussKronrodRule& coarse = FELevelSetMoments::GaussKronrodRules::Get(2U*(primary.gauss/2U)+1U);
  for(unsigned ulps = 1; ulps < (1U<<26U); ulps *= 2U) {
    const long double delta = ulps*std::numeric_limits<long double>::epsilon();
    const FELevelSetMoments::OuterInterval interval = {1,1+delta,1,delta,1,2,true};
    if(!FELevelSetMoments::OuterIntervalQuadrature::Resolved(primary,interval,0,1) && FELevelSetMoments::OuterIntervalQuadrature::Resolved(coarse,interval,0,1)) { RuleSelectionFallback(dimension,selected,ulps,coarse.kronrod); return; }
  }
  throw std::runtime_error("A selected rule did not expose a downward precision transition.");
}

std::vector<double> RuleSelectionNodal(unsigned element, bool ball) {
  const unsigned dimension = element < 2 ? 2U : 3U;
  const bool simplex = element%2U != 0;
  femus::QuadBiquadratic quad;
  femus::TriBiquadratic tri;
  femus::HexBiquadratic hex;
  femus::TetBiquadratic tet;
  femus::basis* basis = element == 0 ? static_cast<femus::basis*>(&quad) : (element == 1 ? static_cast<femus::basis*>(&tri) : (element == 2 ? static_cast<femus::basis*>(&hex) : static_cast<femus::basis*>(&tet)));
  const unsigned count[4] = {9,7,27,15};
  std::vector<double> nodal(count[element]);
  for(unsigned node = 0; node < nodal.size(); ++node) {
    long double value = ball ? 1.0L/64 : 0.4375L;
    for(unsigned axis = 0; axis < dimension; ++axis) {
      const long double coordinate = simplex ? basis->GetXcoarse(node)[axis] : (basis->GetXcoarse(node)[axis]+1)/2;
      value -= ball ? (coordinate-0.25L)*(coordinate-0.25L) : coordinate;
    }
    nodal[node] = static_cast<double>(value);
  }
  return nodal;
}

void RuleSelectionIntegrate(unsigned element, bool ball, const FELevelSetMoments::OuterQuadratureRules& rules) {
  const unsigned dimension = element < 2 ? 2U : 3U;
  const std::vector<double> nodal = RuleSelectionNodal(element,ball);
  std::vector<long double> values;
  std::vector<long double> errors;
  const double tolerance = 1e-3;
  for(unsigned dim = 2; dim <= 3; ++dim) FELevelSetMoments::OuterIntervalQuadrature::Statistics(dim) = FELevelSetMoments::OuterIntervalStatistics();
  if(element%2U == 0) FELevelSetMoments::IntegrateUnitTensorNodalMoments(nodal,dimension,0,values,errors,tolerance,tolerance,false,nullptr,0,dimension == 2 ? 2U : 1U,rules);
  else if(element == 1) {
    std::vector<std::vector<double>> h;
    std::vector<std::vector<double>> d;
    FELevelSetMoments::IntegrateTri7NodalSelected(nodal,0,h,d,tolerance,tolerance,true,false,2,0,nullptr,false,nullptr,nullptr,rules);
    values = {h[0][0],d[0][0]};
  }
  else {
    std::vector<std::vector<std::vector<double>>> h;
    std::vector<std::vector<std::vector<double>>> d;
    FELevelSetMoments::IntegrateTet15NodalByTriangleSlices(nodal,0,h,d,tolerance,tolerance,true,24,0,200000,nullptr,1,nullptr,false,nullptr,nullptr,rules);
    values = {h[0][0][0],d[0][0][0]};
  }
  const long double volume = ball ? std::pow(std::acos(-1.0L),0.5L*dimension)*std::pow(0.125L,dimension)/std::tgamma(0.5L*dimension+1) : std::pow(0.4375L,dimension)/std::tgamma(dimension+1.0L);
  const long double delta = ball ? dimension*volume/(2*0.125L*0.125L) : dimension*volume/0.4375L;
  RuleSelectionCheck(std::abs(values[0]-volume) <= tolerance*(1+volume),"A configured element volume failed its analytic reference.");
  RuleSelectionCheck(std::abs(values[1]-delta) <= tolerance*(1+delta),"A configured element delta failed its analytic reference.");
  for(unsigned dim = 2; dim <= dimension; ++dim) {
    const FELevelSetMoments::OuterIntervalStatistics& stat = FELevelSetMoments::OuterIntervalQuadrature::Statistics(dim);
    const unsigned regular = FELevelSetMoments::GaussKronrodRules::Index(rules.Select(dim,false));
    const unsigned tangent = FELevelSetMoments::GaussKronrodRules::Index(rules.Select(dim,true));
    for(unsigned index = 0; index < 15; ++index) {
      if(index != regular) RuleSelectionCheck(stat.regularRulePanels[index] == 0,"A nested regular integral ignored its configured rule.");
      if(index != tangent) RuleSelectionCheck(stat.transformedRulePanels[index] == 0,"A nested tangent integral ignored its configured rule.");
    }
    if(ball) RuleSelectionCheck(stat.transformedRulePanels[tangent] > 0,"The requested tangent rule was not exercised.");
    else RuleSelectionCheck(stat.regularRulePanels[regular] > 0,"The requested regular rule was not exercised.");
  }
}

void RuleSelectionChildResolution(unsigned dimension, unsigned selected, bool mapped, bool unrelated) {
  FELevelSetMoments::OuterQuadratureRules rules;
  rules.regular2D = rules.tangent2D = rules.regular3D = rules.tangent3D = selected;
  const FELevelSetMoments::GaussKronrodRule& primary = FELevelSetMoments::GaussKronrodRules::Get(selected);
  const FELevelSetMoments::GaussKronrodRule& coarse = FELevelSetMoments::GaussKronrodRules::Get(2U*(primary.gauss/2U)+1U);
  const FELevelSetMoments::OuterInterval interval = mapped ? FELevelSetMoments::OuterInterval{0,1,0,1,1,2,true} : FELevelSetMoments::OuterInterval{0,1,0,0,0,1,true};
  RuleSelectionContext context;
  context.mode = unrelated ? 4U : 3U;
  context.cutoff = 0.5L-0.25L*(primary.nodes[primary.size-1U]+coarse.nodes[coarse.size-1U]);
  std::vector<unsigned char> active(2,1);
  std::vector<long double> values(2,0);
  std::vector<long double> errors(2,0);
  unsigned accepted = 0;
  bool rejected = false;
  FELevelSetMoments::OuterIntervalQuadrature::Statistics(dimension) = FELevelSetMoments::OuterIntervalStatistics();
  try { FELevelSetMoments::OuterIntervalQuadrature::Adaptive(context,RuleSelectionSample,interval,0,1,1e-12L,1e-12L,1e-12L,0,24,10000,context.evaluations,accepted,dimension,active,values,errors,rules); }
  catch(const std::runtime_error&) { rejected = true; }
  const FELevelSetMoments::OuterIntervalStatistics& statistics = FELevelSetMoments::OuterIntervalQuadrature::Statistics(dimension);
  RuleSelectionCheck(rejected == unrelated,"A nested precision failure was confused with an unrelated failure.");
  RuleSelectionCheck(statistics.subdivisions == 0,"A failed inner slice caused repeated outer bisection.");
  for(unsigned index = primary.index; index < 15; ++index) RuleSelectionCheck(statistics.emergencyRulePanels[index] == 0,"A failed inner slice raised its parent's selected rule.");
  if(unrelated) RuleSelectionCheck(statistics.childResolutionFailures == 0 && accepted == 0,"A non-precision failure entered the failsafe.");
  else {
    RuleSelectionCheck(statistics.childResolutionFailures == 1 && statistics.emergencyRulePanels[coarse.index] == 1,"A failed inner slice did not reduce its parent rule.");
    RuleSelectionCheck(std::abs(values[0]-1) < 1e-12L && std::abs(values[1]-0.5L) < 1e-12L,"A parent retry retained partial moments from its failed panel.");
  }
}

void RuleSelectionCurvedPrecision() {
  const std::vector<double> nodal = {-0.0624999999999999861222,-0.125000000000000083267,-0.0624999999999999861222,-2.77555756156289135106e-17,-1.20370621524202240816e-35,-0.0625,-1.20370621524202240816e-35,0.0625,-0.0781250000000000277556,-0.0781250000000000277556,-0.0156249999999999930611,-0.0156249999999999930611,-0.015625,-0.015625,0.046875,0.046875,-0.0156250000000000173472,-0.0781250000000000555112,-0.0156250000000000173472,0.0468749999999999930611,-0.0312499999999999861222,-0.03125,0.0312500000000000069389,0.0312500000000000069389,-0.0312499999999999861222,0.03125,0.0156249999999999930611};
  const std::vector<std::vector<double>> coordinates = {{0,0.200000000000000011102,0.200000000000000011102,0,0,0.25,0.25,0,0.100000000000000019429,0.200000000000000038858,0.100000000000000019429,0,0.125,0.25,0.125,0,0,0.225000000000000033307,0.225000000000000033307,0,0.112500000000000002776,0.225000000000000005551,0.112500000000000002776,0,0.100000000000000005551,0.125,0.112500000000000016653},{-0.200000000000000011102,-0.200000000000000011102,0,0,-0.25,-0.25,0,0,-0.200000000000000038858,-0.100000000000000019429,0,-0.100000000000000019429,-0.25,-0.125,0,-0.125,-0.225000000000000033307,-0.225000000000000033307,0,0,-0.224999999999999977796,-0.112500000000000002776,0,-0.112500000000000002776,-0.0999999999999999916733,-0.125,-0.112500000000000016653},{-0.249999999999999972244,-0.275000000000000022204,-0.249999999999999972244,-0.225000000000000033307,-0.0250000000000000048572,-0.0499999999999999888978,-0.0250000000000000048572,0,-0.256249999999999977796,-0.256249999999999977796,-0.231249999999999983347,-0.231249999999999983347,-0.0312499999999999895917,-0.0312499999999999895917,-0.00625000000000000208167,-0.00625000000000000208167,-0.143749999999999961142,-0.168750000000000011102,-0.143749999999999961142,-0.118750000000000008327,-0.149999999999999994449,-0.149999999999999994449,-0.124999999999999972244,-0.124999999999999972244,-0.237499999999999988898,-0.0124999999999999868161,-0.131249999999999977796}};
  const long double exactH = std::acos(-1.0L)/6;
  const long double exactSurface = std::acos(-1.0L)/32*(0.8L*(std::pow(1.04L,1.5L)-std::pow(0.64L,1.5L))+0.5L*(std::pow(1.04L,2.5L)-std::pow(0.64L,2.5L)));
  std::vector<long double> baseline;
  for(unsigned choice = 0; choice < 4; ++choice) {
    FELevelSetMoments::OuterQuadratureRules rules;
    if(choice == 1) rules.regular2D = rules.tangent2D = rules.regular3D = rules.tangent3D = 31;
    if(choice == 2) { rules.regular2D = rules.regular3D = 21; rules.tangent2D = rules.tangent3D = 31; }
    if(choice == 3) rules.tangent2D = 31;
    std::vector<long double> values;
    std::vector<long double> errors;
    FELevelSetMoments::OuterIntervalQuadrature::Statistics(3) = FELevelSetMoments::OuterIntervalStatistics();
    FELevelSetMoments::IntegrateUnitTensorNodalMoments(nodal,3,8,values,errors,1e-5L,1e-5L,true,&coordinates,0,0,rules);
    const unsigned count = values.size()/2U;
    RuleSelectionCheck(std::abs(values[0]-exactH) < 5e-12L && std::abs(values[count]-exactSurface) < 5e-12L,"A refined curved HEX failed its analytic volume or surface reference.");
    if(choice == 0) baseline = values;
    for(unsigned n = 0; n < values.size(); ++n) RuleSelectionCheck(std::abs(values[n]-baseline[n]) < 1e-9L,"A nested precision retry changed a curved HEX moment beyond tolerance.");
    if(choice == 1 || choice == 2) RuleSelectionCheck(FELevelSetMoments::OuterIntervalQuadrature::Statistics(3).childResolutionFailures > 0,"The refined HEX did not exercise the parent precision retry.");
  }
}

void RuleSelectionWeightCache() {
  const std::vector<double> nodal = RuleSelectionNodal(0,false);
  CutFemWeight_LS<long double,long double> first(QUAD,2,"legendre");
  CutFemWeight_LS<long double,long double> second(QUAD,2,"legendre");
  FELevelSetMoments::OuterQuadratureRules selected;
  selected.regular2D = 31;
  std::vector<long double> h;
  std::vector<long double> d;
  first.GetWeightsFromNodalValues(nodal,h,d,1e-6,1e-6,false);
  FELevelSetMoments::OuterIntervalQuadrature::Statistics(2) = FELevelSetMoments::OuterIntervalStatistics();
  first.GetWeightsFromNodalValues(nodal,h,d,1e-6,1e-6,false,selected);
  RuleSelectionCheck(FELevelSetMoments::OuterIntervalQuadrature::Statistics(2).regularRulePanels[FELevelSetMoments::GaussKronrodRules::Index(31)] > 0,"The weight cache ignored a per-call rule change.");
  RuleSelectionCheck(first.GetOuterQuadratureRules().regular2D == 15,"A per-call override changed the object's defaults.");
  FELevelSetMoments::OuterIntervalQuadrature::Statistics(2) = FELevelSetMoments::OuterIntervalStatistics();
  first.GetWeightsFromNodalValues(nodal,h,d,1e-6,1e-6,false,selected);
  RuleSelectionCheck(FELevelSetMoments::OuterIntervalQuadrature::Statistics(2).samples == 0,"Identical configured weights were not cached.");
  first.SetOuterQuadratureRules(selected);
  first.GetWeightsFromNodalValues(nodal,h,d,1e-6,1e-6,false);
  FELevelSetMoments::OuterIntervalQuadrature::Statistics(2) = FELevelSetMoments::OuterIntervalStatistics();
  second.GetWeightsFromNodalValues(nodal,h,d,1e-6,1e-6,false);
  RuleSelectionCheck(FELevelSetMoments::OuterIntervalQuadrature::Statistics(2).regularRulePanels[FELevelSetMoments::GaussKronrodRules::Index(15)] > 0,"Rule settings leaked between weight objects.");
  first.GetWeightsFromNodalValues(nodal,h,d,1e-6,1e-6,false,FELevelSetMoments::OuterQuadratureRules());
  FELevelSetMoments::OuterIntervalQuadrature::Statistics(2) = FELevelSetMoments::OuterIntervalStatistics();
  first.GetWeightsFromNodalValues(nodal,h,d,1e-6,1e-6,false);
  RuleSelectionCheck(FELevelSetMoments::OuterIntervalQuadrature::Statistics(2).regularRulePanels[FELevelSetMoments::GaussKronrodRules::Index(31)] > 0,"A per-call override contaminated the persistent configured cache.");
  selected.tangent3D = 99;
  bool rejected = false;
  try { first.SetOuterQuadratureRules(selected); } catch(const std::invalid_argument&) { rejected = true; }
  RuleSelectionCheck(rejected && first.GetOuterQuadratureRules().tangent3D == 15,"Invalid settings changed a weight object's configuration.");
  FELevelSetMoments::OuterQuadratureRules lineRules;
  lineRules.regular2D = lineRules.tangent2D = lineRules.regular3D = lineRules.tangent3D = 3;
  std::vector<long double> values;
  std::vector<long double> errors;
  FELevelSetMoments::OuterIntervalQuadrature::Statistics(2) = FELevelSetMoments::OuterIntervalStatistics();
  FELevelSetMoments::OuterIntervalQuadrature::Statistics(3) = FELevelSetMoments::OuterIntervalStatistics();
  FELevelSetMoments::IntegrateUnitTensorNodalMoments({-0.375,0.625,0.125},1,10,values,errors,1e-12L,1e-12L,false,nullptr,0,1,lineRules);
  for(unsigned degree = 0; degree <= 10; ++degree) {
    RuleSelectionCheck(std::abs(values[degree]-(1-FELevelSetMoments::ReferencePower(0.375L,degree+1U))/(degree+1U)) < 1e-16L,"Rule selection changed exact LINE3 integration.");
    RuleSelectionCheck(std::abs(values[11U+degree]-FELevelSetMoments::ReferencePower(0.375L,degree)) < 1e-16L,"Rule selection changed the exact LINE3 delta.");
  }
  RuleSelectionCheck(FELevelSetMoments::OuterIntervalQuadrature::Statistics(2).samples == 0 && FELevelSetMoments::OuterIntervalQuadrature::Statistics(3).samples == 0,"Exact LINE3 integration invoked outer quadrature.");
}

int main() {
  try {
    RuleSelectionTables();
    const unsigned orders[4] = {3,7,15,31};
    for(unsigned dimension = 2; dimension <= 3; ++dimension) {
      for(unsigned gauss = 1; gauss <= 15; ++gauss) for(unsigned mapped = 0; mapped < 2; ++mapped) RuleSelectionIntervals(dimension,2U*gauss+1U,mapped != 0);
      for(unsigned gauss = 2; gauss <= 15; ++gauss) RuleSelectionConfiguredFallback(dimension,2U*gauss+1U);
      for(unsigned selected : {21U,31U}) for(unsigned mapped = 0; mapped < 2; ++mapped) for(unsigned unrelated = 0; unrelated < 2; ++unrelated) RuleSelectionChildResolution(dimension,selected,mapped != 0,unrelated != 0);
      RuleSelectionFallback(dimension,31,65536,15);
      RuleSelectionFallback(dimension,31,4096,7);
      RuleSelectionFallback(dimension,31,128,3);
      RuleSelectionFallback(dimension,15,4096,7);
      RuleSelectionFallback(dimension,7,128,3);
      RuleSelectionFallback(dimension,3,1,0);
      RuleSelectionFallback(dimension,31,65536,0,true);
    }
    for(unsigned first = 0; first < 4; ++first) for(unsigned second = 0; second < 4; ++second) {
      FELevelSetMoments::OuterQuadratureRules rules;
      rules.regular2D = orders[first];
      rules.regular3D = orders[second];
      for(unsigned element = 2; element < 4; ++element) RuleSelectionIntegrate(element,false,rules);
      rules.regular2D = rules.regular3D = 31;
      rules.tangent2D = orders[first];
      rules.tangent3D = orders[second];
      for(unsigned element = 2; element < 4; ++element) RuleSelectionIntegrate(element,true,rules);
      if(second == 0) for(unsigned element = 0; element < 2; ++element) { RuleSelectionIntegrate(element,true,rules); rules.regular2D = orders[first]; RuleSelectionIntegrate(element,false,rules); }
    }
    for(unsigned gauss = 1; gauss <= 15; ++gauss) {
      FELevelSetMoments::OuterQuadratureRules rules;
      rules.regular2D = rules.tangent2D = 2U*gauss+1U;
      rules.regular3D = rules.tangent3D = 2U*(16U-gauss)+1U;
      for(unsigned element = 0; element < 4; ++element) for(unsigned ball = 0; ball < 2; ++ball) RuleSelectionIntegrate(element,ball != 0,rules);
    }
    RuleSelectionCurvedPrecision();
    RuleSelectionWeightCache();
    std::cout << ruleSelectionChecks << " configurable-rule, nested-slice, cache, exactness and downward-failsafe checks passed\n";
  }
  catch(const std::exception& error) { std::cerr << error.what() << '\n'; return 1; }
}
