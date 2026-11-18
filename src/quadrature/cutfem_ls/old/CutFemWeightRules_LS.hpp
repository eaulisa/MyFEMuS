#ifndef FEMUS_CUT_FEM_WEIGHT_RULES_LS_HPP
#define FEMUS_CUT_FEM_WEIGHT_RULES_LS_HPP

#include "CutFemWeight_LS.hpp"
#include <memory>

template<class TypeIO, class TypeA>
class CutFemWeightRules_LS {
public:
  // Constructs the selected rules once; element assembly only retrieves existing objects.
  CutFemWeightRules_LS(unsigned order, const std::vector<GeomElType>& geometries, const FELevelSetMoments::OuterQuadratureRules& outerRules, const std::string& gaussType = "legendre", unsigned directionValidation = 0) : _order(order), _outerRules(outerRules), _rules(6) {
    _outerRules.Validate();
    for(unsigned n = 0; n < geometries.size(); ++n) {
      const unsigned geometry = static_cast<unsigned>(geometries[n]);
      if(geometry >= _rules.size()) throw std::invalid_argument("Invalid fixed CutFEM rule geometry.");
      if(_rules[geometry]) continue;
      _rules[geometry].reset(new CutFemWeight_LS<TypeIO,TypeA>(geometries[n],_order,gaussType));
      _rules[geometry]->SetIntegrationDirection(0);
      _rules[geometry]->SetDirectionValidation(directionValidation);
      _rules[geometry]->SetOuterQuadratureRules(_outerRules);
    }
  }

  CutFemWeight_LS<TypeIO,TypeA>* GetRule(GeomElType geometry) const {
    const unsigned index = static_cast<unsigned>(geometry);
    if(index >= _rules.size() || !_rules[index]) throw std::logic_error("The CutFEM rule must be initialized before assembly.");
    return _rules[index].get();
  }

  unsigned GetOrder() const {
    return _order;
  }

  const FELevelSetMoments::OuterQuadratureRules& GetOuterRules() const {
    return _outerRules;
  }

private:
  const unsigned _order;
  const FELevelSetMoments::OuterQuadratureRules _outerRules;
  std::vector<std::unique_ptr<CutFemWeight_LS<TypeIO,TypeA>>> _rules;
};

#endif
