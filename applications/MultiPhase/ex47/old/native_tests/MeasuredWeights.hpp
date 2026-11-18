#ifndef FEMUS_EX47_NATIVE_MEASURED_WEIGHTS_HPP
#define FEMUS_EX47_NATIVE_MEASURED_WEIGHTS_HPP

#include "../include/CutFemWeight_LS.hpp"
#include "CutFemWeightRules_LS.hpp"
#include <chrono>
#include <iomanip>
#include <iostream>
#include <map>

namespace Ex47NativeTests {

struct WorkCounters {
  unsigned long long evaluations = 0;
  unsigned long long panels = 0;
  unsigned long long subdivisions = 0;
  unsigned long long outer_samples = 0;
  unsigned long long integration_passes = 0;
  unsigned long long graph_cells = 0;
  unsigned long long general_surfaces = 0;
  unsigned long long line_restrictions = 0;
  unsigned long long metric_preparations = 0;
  unsigned long long coefficient_expansions = 0;
};

inline WorkCounters Counters() {
  WorkCounters result;
  const FELevelSetMoments::OuterIntervalStatistics& first = FELevelSetMoments::OuterIntervalQuadrature::Statistics(2);
  const FELevelSetMoments::OuterIntervalStatistics& second = FELevelSetMoments::OuterIntervalQuadrature::Statistics(3);
  result.evaluations = first.samples+second.samples;
  result.panels = first.regularPanels+first.transformedPanels+second.regularPanels+second.transformedPanels;
  result.subdivisions = first.subdivisions+second.subdivisions;
#ifdef EX47_NATIVE_PREPARED
  const FELevelSetMoments::Ex47Prepared::SurfaceEndpointStatistics& surface = FELevelSetMoments::Ex47Prepared::SurfaceEndpointIntegration::Statistics();
#else
  const FELevelSetMoments::SurfaceEndpointStatistics& surface = FELevelSetMoments::SurfaceEndpointIntegration::Statistics();
#endif
  result.outer_samples = surface.outerSamples;
  result.integration_passes = surface.integrationPasses;
  result.graph_cells = surface.graphCells;
  result.general_surfaces = surface.generalSurfaces;
  return result;
}

template<class TypeIO, class TypeA> WorkCounters ObjectCounters(const ::CutFemWeight_LS<TypeIO,TypeA>& rule) { return Counters(); }
template<class TypeIO, class TypeA> WorkCounters ObjectCounters(const Ex47::CutFemWeight_LS<TypeIO,TypeA>& rule) { WorkCounters result = Counters(); result.line_restrictions = rule.GetPreparedStatistics().lineRestrictions; result.metric_preparations = rule.GetPreparedStatistics().geometryPreparations; try { result.coefficient_expansions = rule.GetPreparedLevelSetData().GetWorkCounts()[0]; } catch(const std::logic_error&) {} return result; }

struct Timing {
  unsigned long long calls = 0;
  unsigned long long failures = 0;
  long double seconds = 0;
  WorkCounters work;
};

inline std::map<std::string,Timing>& Timings() { static std::map<std::string,Timing> values; return values; }

inline void Record(const char* method, const std::chrono::steady_clock::time_point& start, bool success, const WorkCounters* before = nullptr, const WorkCounters* after = nullptr, const std::chrono::steady_clock::time_point* stop = nullptr) {
  const long double elapsed = std::chrono::duration<long double>((stop == nullptr ? std::chrono::steady_clock::now() : *stop)-start).count();
  Timing& value = Timings()[method];
  ++value.calls;
  value.failures += !success;
  value.seconds += elapsed;
  if(before == nullptr) return;
  const WorkCounters current = after == nullptr ? Counters() : *after;
  if(current.evaluations >= before->evaluations) value.work.evaluations += current.evaluations-before->evaluations;
  if(current.panels >= before->panels) value.work.panels += current.panels-before->panels;
  if(current.subdivisions >= before->subdivisions) value.work.subdivisions += current.subdivisions-before->subdivisions;
  if(current.outer_samples >= before->outer_samples) value.work.outer_samples += current.outer_samples-before->outer_samples;
  if(current.integration_passes >= before->integration_passes) value.work.integration_passes += current.integration_passes-before->integration_passes;
  if(current.graph_cells >= before->graph_cells) value.work.graph_cells += current.graph_cells-before->graph_cells;
  if(current.general_surfaces >= before->general_surfaces) value.work.general_surfaces += current.general_surfaces-before->general_surfaces;
  if(current.line_restrictions >= before->line_restrictions) value.work.line_restrictions += current.line_restrictions-before->line_restrictions;
  if(current.metric_preparations >= before->metric_preparations) value.work.metric_preparations += current.metric_preparations-before->metric_preparations;
  if(current.coefficient_expansions >= before->coefficient_expansions) value.work.coefficient_expansions += current.coefficient_expansions-before->coefficient_expansions;
}

template<class Weight> void RecordWeight(const char* method, const std::chrono::steady_clock::time_point& start, bool success, const WorkCounters& before, const Weight& rule) { const std::chrono::steady_clock::time_point stop = std::chrono::steady_clock::now(); const WorkCounters after = ObjectCounters(rule); Record(method,start,success,&before,&after,&stop); }

inline std::chrono::steady_clock::time_point& ConstructionStart() { static thread_local std::chrono::steady_clock::time_point start; return start; }
inline unsigned BeginConstruction(unsigned degree) { ConstructionStart() = std::chrono::steady_clock::now(); return degree; }

struct TimingOutput {
  TimingOutput() { Timings(); FELevelSetMoments::GaussKronrodRules::All();
#ifdef EX47_NATIVE_PREPARED
    FELevelSetMoments::Ex47Prepared::LineIntegration::FullStripRules();
#endif
  }
  ~TimingOutput() {
    std::cerr << "EX47_NATIVE_TIMING {";
    bool first = true;
    for(const std::pair<const std::string,Timing>& entry : Timings()) {
      if(!first) std::cerr << ',';
      first = false;
      std::cerr << '"' << entry.first << "\":{\"calls\":" << entry.second.calls << ",\"failures\":" << entry.second.failures << ",\"seconds\":" << std::setprecision(18) << entry.second.seconds;
      std::cerr << ",\"evaluations\":" << entry.second.work.evaluations;
      std::cerr << ",\"panels\":" << entry.second.work.panels;
      std::cerr << ",\"subdivisions\":" << entry.second.work.subdivisions;
      std::cerr << ",\"outer_samples\":" << entry.second.work.outer_samples;
      std::cerr << ",\"integration_passes\":" << entry.second.work.integration_passes;
      std::cerr << ",\"graph_cells\":" << entry.second.work.graph_cells;
      std::cerr << ",\"general_surfaces\":" << entry.second.work.general_surfaces;
      std::cerr << ",\"line_restrictions\":" << entry.second.work.line_restrictions;
      std::cerr << ",\"metric_preparations\":" << entry.second.work.metric_preparations;
      std::cerr << ",\"coefficient_expansions\":" << entry.second.work.coefficient_expansions;
      std::cerr << '}';
    }
    std::cerr << "}\n";
  }
};

inline TimingOutput& Printer() { static TimingOutput output; return output; }
static const TimingOutput& timingOutput = Printer();

template<class TypeIO, class TypeA, template<class,class> class Weight> class MeasuredWeight : public Weight<TypeIO,TypeA> {
public:
  MeasuredWeight(GeomElType geometry, unsigned degree, const std::string& family) try : Weight<TypeIO,TypeA>(geometry,BeginConstruction(degree),family) { Record("construction",ConstructionStart(),true); } catch(...) { Record("construction",ConstructionStart(),false); throw; }
  void GetWeightsFromNodalValues(const std::vector<double>& nodal, std::vector<TypeIO>& positive, std::vector<TypeIO>& surface, double absolute = 1e-12, double relative = 1e-12, bool reference = false) { const WorkCounters before = ObjectCounters(static_cast<const Weight<TypeIO,TypeA>&>(*this)); const std::chrono::steady_clock::time_point start = std::chrono::steady_clock::now(); try { Weight<TypeIO,TypeA>::GetWeightsFromNodalValues(nodal,positive,surface,absolute,relative,reference); } catch(...) { RecordWeight("weights",start,false,before,static_cast<const Weight<TypeIO,TypeA>&>(*this)); throw; } RecordWeight("weights",start,true,before,static_cast<const Weight<TypeIO,TypeA>&>(*this)); }
  void GetWeightsFromNodalValues(const std::vector<double>& nodal, std::vector<TypeIO>& positive, std::vector<TypeIO>& surface, double absolute, double relative, bool reference, const FELevelSetMoments::OuterQuadratureRules& rules) { const WorkCounters before = ObjectCounters(static_cast<const Weight<TypeIO,TypeA>&>(*this)); const std::chrono::steady_clock::time_point start = std::chrono::steady_clock::now(); try { Weight<TypeIO,TypeA>::GetWeightsFromNodalValues(nodal,positive,surface,absolute,relative,reference,rules); } catch(...) { RecordWeight("weights",start,false,before,static_cast<const Weight<TypeIO,TypeA>&>(*this)); throw; } RecordWeight("weights",start,true,before,static_cast<const Weight<TypeIO,TypeA>&>(*this)); }
  void GetPhysicalSurfaceWeightsFromNodalGeometry(const std::vector<double>& nodal, const std::vector<std::vector<double>>& coordinates, std::vector<TypeIO>& positive, std::vector<TypeIO>& surface, double absolute = 1e-12, double relative = 1e-12) { const WorkCounters before = ObjectCounters(static_cast<const Weight<TypeIO,TypeA>&>(*this)); const std::chrono::steady_clock::time_point start = std::chrono::steady_clock::now(); try { Weight<TypeIO,TypeA>::GetPhysicalSurfaceWeightsFromNodalGeometry(nodal,coordinates,positive,surface,absolute,relative); } catch(...) { RecordWeight("physical",start,false,before,static_cast<const Weight<TypeIO,TypeA>&>(*this)); throw; } RecordWeight("physical",start,true,before,static_cast<const Weight<TypeIO,TypeA>&>(*this)); }
  void GetPhysicalSurfaceWeightsFromNodalGeometry(const std::vector<double>& nodal, const std::vector<std::vector<double>>& coordinates, std::vector<TypeIO>& positive, std::vector<TypeIO>& surface, double absolute, double relative, const FELevelSetMoments::OuterQuadratureRules& rules) { const WorkCounters before = ObjectCounters(static_cast<const Weight<TypeIO,TypeA>&>(*this)); const std::chrono::steady_clock::time_point start = std::chrono::steady_clock::now(); try { Weight<TypeIO,TypeA>::GetPhysicalSurfaceWeightsFromNodalGeometry(nodal,coordinates,positive,surface,absolute,relative,rules); } catch(...) { RecordWeight("physical",start,false,before,static_cast<const Weight<TypeIO,TypeA>&>(*this)); throw; } RecordWeight("physical",start,true,before,static_cast<const Weight<TypeIO,TypeA>&>(*this)); }
  void GetTangentialSurfaceWeightsFromNodalGeometry(const std::vector<double>& nodal, const std::vector<std::vector<double>>& coordinates, std::vector<TypeIO>& positive, std::vector<TypeIO>& surface, double absolute = 1e-12, double relative = 1e-12) { const WorkCounters before = ObjectCounters(static_cast<const Weight<TypeIO,TypeA>&>(*this)); const std::chrono::steady_clock::time_point start = std::chrono::steady_clock::now(); try { Weight<TypeIO,TypeA>::GetTangentialSurfaceWeightsFromNodalGeometry(nodal,coordinates,positive,surface,absolute,relative); } catch(...) { RecordWeight("tangential",start,false,before,static_cast<const Weight<TypeIO,TypeA>&>(*this)); throw; } RecordWeight("tangential",start,true,before,static_cast<const Weight<TypeIO,TypeA>&>(*this)); }
  void GetTangentialSurfaceWeightsFromNodalGeometry(const std::vector<double>& nodal, const std::vector<std::vector<double>>& coordinates, std::vector<TypeIO>& positive, std::vector<TypeIO>& surface, double absolute, double relative, const FELevelSetMoments::OuterQuadratureRules& rules) { const WorkCounters before = ObjectCounters(static_cast<const Weight<TypeIO,TypeA>&>(*this)); const std::chrono::steady_clock::time_point start = std::chrono::steady_clock::now(); try { Weight<TypeIO,TypeA>::GetTangentialSurfaceWeightsFromNodalGeometry(nodal,coordinates,positive,surface,absolute,relative,rules); } catch(...) { RecordWeight("tangential",start,false,before,static_cast<const Weight<TypeIO,TypeA>&>(*this)); throw; } RecordWeight("tangential",start,true,before,static_cast<const Weight<TypeIO,TypeA>&>(*this)); }
  void GetNormalSurfaceWeightsFromNodalGeometry(const std::vector<double>& nodal, const std::vector<std::vector<double>>& coordinates, std::vector<TypeIO>& positive, std::vector<std::vector<TypeIO>>& surface, double absolute = 1e-12, double relative = 1e-12) { const WorkCounters before = ObjectCounters(static_cast<const Weight<TypeIO,TypeA>&>(*this)); const std::chrono::steady_clock::time_point start = std::chrono::steady_clock::now(); try { Weight<TypeIO,TypeA>::GetNormalSurfaceWeightsFromNodalGeometry(nodal,coordinates,positive,surface,absolute,relative); } catch(...) { RecordWeight("normals",start,false,before,static_cast<const Weight<TypeIO,TypeA>&>(*this)); throw; } RecordWeight("normals",start,true,before,static_cast<const Weight<TypeIO,TypeA>&>(*this)); }
  void GetNormalSurfaceWeightsFromNodalGeometry(const std::vector<double>& nodal, const std::vector<std::vector<double>>& coordinates, std::vector<TypeIO>& positive, std::vector<std::vector<TypeIO>>& surface, double absolute, double relative, const FELevelSetMoments::OuterQuadratureRules& rules) { const WorkCounters before = ObjectCounters(static_cast<const Weight<TypeIO,TypeA>&>(*this)); const std::chrono::steady_clock::time_point start = std::chrono::steady_clock::now(); try { Weight<TypeIO,TypeA>::GetNormalSurfaceWeightsFromNodalGeometry(nodal,coordinates,positive,surface,absolute,relative,rules); } catch(...) { RecordWeight("normals",start,false,before,static_cast<const Weight<TypeIO,TypeA>&>(*this)); throw; } RecordWeight("normals",start,true,before,static_cast<const Weight<TypeIO,TypeA>&>(*this)); }
};

template<class TypeIO, class TypeA> class OriginalWeight : public MeasuredWeight<TypeIO,TypeA,::CutFemWeight_LS> { public: OriginalWeight(GeomElType geometry, unsigned degree, const std::string& family) : MeasuredWeight<TypeIO,TypeA,::CutFemWeight_LS>(geometry,degree,family) {} };
template<class TypeIO, class TypeA> class PreparedWeight : public MeasuredWeight<TypeIO,TypeA,Ex47::CutFemWeight_LS> { public: PreparedWeight(GeomElType geometry, unsigned degree, const std::string& family) : MeasuredWeight<TypeIO,TypeA,Ex47::CutFemWeight_LS>(geometry,degree,family) {} };

template<class TypeIO, class TypeA, template<class,class> class Weight> class MeasuredRules {
public:
  MeasuredRules(unsigned order, const std::vector<GeomElType>& geometries, const FELevelSetMoments::OuterQuadratureRules& rules, const std::string& family = "legendre", unsigned validation = 0) : _order(order), _outerRules(rules), _rules(6) { rules.Validate(); for(GeomElType geometry : geometries) { const unsigned index = static_cast<unsigned>(geometry); if(index >= 6) throw std::invalid_argument("Invalid fixed rule geometry."); if(_rules[index]) continue; _rules[index].reset(new Weight<TypeIO,TypeA>(geometry,order,family)); _rules[index]->SetDirectionValidation(validation); _rules[index]->SetOuterQuadratureRules(rules); } }
  Weight<TypeIO,TypeA>* GetRule(GeomElType geometry) const { const unsigned index = static_cast<unsigned>(geometry); if(index >= 6 || !_rules[index]) throw std::logic_error("Rule was not initialized."); return _rules[index].get(); }
  unsigned GetOrder() const { return _order; }
  const FELevelSetMoments::OuterQuadratureRules& GetOuterRules() const { return _outerRules; }
private:
  unsigned _order;
  FELevelSetMoments::OuterQuadratureRules _outerRules;
  std::vector<std::unique_ptr<Weight<TypeIO,TypeA>>> _rules;
};

template<class TypeIO, class TypeA> class OriginalRules : public MeasuredRules<TypeIO,TypeA,OriginalWeight> { public: OriginalRules(unsigned order, const std::vector<GeomElType>& geometries, const FELevelSetMoments::OuterQuadratureRules& rules, const std::string& family = "legendre", unsigned validation = 0) : MeasuredRules<TypeIO,TypeA,OriginalWeight>(order,geometries,rules,family,validation) {} };
template<class TypeIO, class TypeA> class PreparedRules : public MeasuredRules<TypeIO,TypeA,PreparedWeight> { public: PreparedRules(unsigned order, const std::vector<GeomElType>& geometries, const FELevelSetMoments::OuterQuadratureRules& rules, const std::string& family = "legendre", unsigned validation = 0) : MeasuredRules<TypeIO,TypeA,PreparedWeight>(order,geometries,rules,family,validation) {} };

}

#endif
