#ifndef FEMUS_GAUSS_KRONROD_RULES_HPP
#define FEMUS_GAUSS_KRONROD_RULES_HPP

#include <array>
#include <boost/math/quadrature/gauss_kronrod.hpp>
#include <stdexcept>
#include <string>

namespace FELevelSetMoments {

struct GaussKronrodRule {
  unsigned gauss = 0;
  unsigned kronrod = 0;
  unsigned index = 0;
  unsigned size = 0;
  const long double* nodes = nullptr;
  const long double* weights = nullptr;
  const long double* gaussWeights = nullptr;
  const char* name = nullptr;
};

class GaussKronrodRules {
public:
  template<unsigned Kronrod> static GaussKronrodRule Make(unsigned index, const char* name) {
    return {(Kronrod-1U)/2U,Kronrod,index,static_cast<unsigned>(boost::math::quadrature::gauss_kronrod<long double,Kronrod>::abscissa().size()),boost::math::quadrature::gauss_kronrod<long double,Kronrod>::abscissa().data(),boost::math::quadrature::gauss_kronrod<long double,Kronrod>::weights().data(),boost::math::quadrature::gauss<long double,(Kronrod-1U)/2U>::weights().data(),name};
  }

  static const std::array<GaussKronrodRule,15>& All() {
    static const std::array<GaussKronrodRule,15> rules = {{Make<3>(0,"1/3"),Make<5>(1,"2/5"),Make<7>(2,"3/7"),Make<9>(3,"4/9"),Make<11>(4,"5/11"),Make<13>(5,"6/13"),Make<15>(6,"7/15"),Make<17>(7,"8/17"),Make<19>(8,"9/19"),Make<21>(9,"10/21"),Make<23>(10,"11/23"),Make<25>(11,"12/25"),Make<27>(12,"13/27"),Make<29>(13,"14/29"),Make<31>(14,"15/31")}};
    return rules;
  }

  static unsigned Index(unsigned kronrod) {
    if(kronrod >= 3 && kronrod <= 31 && kronrod%2U != 0) return (kronrod-3U)/2U;
    throw std::invalid_argument("A Gauss-Kronrod rule must be G_n/K_(2n+1) for n=1,...,15 (odd Kronrod order 3,...,31).");
  }

  static const GaussKronrodRule& Get(unsigned kronrod) {
    return All()[Index(kronrod)];
  }

  static unsigned Parse(const std::string& pair) {
    for(unsigned i = 0; i < All().size(); ++i) {
      const std::string slash = All()[i].name;
      std::string dash = slash;
      dash[dash.find('/')] = '-';
      if(pair == slash || pair == dash) return All()[i].kronrod;
    }
    throw std::invalid_argument("Invalid Gauss-Kronrod pair: "+pair+"; expected n/(2n+1) for n=1,...,15.");
  }

  static bool IsOption(const std::string& option) {
    return option == "--gk-regular-2d" || option == "--gk-tangent-2d" || option == "--gk-regular-3d" || option == "--gk-tangent-3d";
  }

  static const char* Help() {
    return "  --gk-regular-2d PAIR         2D regular intervals (default 7/15)\n  --gk-tangent-2d PAIR         2D vertical-tangent intervals (default 7/15)\n  --gk-regular-3d PAIR         3D regular intervals (default 7/15)\n  --gk-tangent-3d PAIR         3D vertical-tangent intervals (default 7/15)\n  PAIR: n/(2n+1), n=1,...,15 (1/3, 2/5, ..., 15/31); dash notation is also accepted.\n";
  }
};

struct OuterQuadratureRules {
  unsigned regular2D = 15;
  unsigned tangent2D = 15;
  unsigned regular3D = 15;
  unsigned tangent3D = 15;

  void Validate() const {
    GaussKronrodRules::Index(regular2D);
    GaussKronrodRules::Index(tangent2D);
    GaussKronrodRules::Index(regular3D);
    GaussKronrodRules::Index(tangent3D);
  }

  bool SameAs(const OuterQuadratureRules& other) const {
    return regular2D == other.regular2D && tangent2D == other.tangent2D && regular3D == other.regular3D && tangent3D == other.tangent3D;
  }

  unsigned Select(unsigned dimension, bool tangent) const {
    if(dimension == 2) return tangent ? tangent2D : regular2D;
    if(dimension == 3) return tangent ? tangent3D : regular3D;
    throw std::invalid_argument("Outer Gauss-Kronrod rules require dimension 2 or 3.");
  }

  void SetOption(const std::string& option, const std::string& pair) {
    const unsigned rule = GaussKronrodRules::Parse(pair);
    if(option == "--gk-regular-2d") regular2D = rule;
    else if(option == "--gk-tangent-2d") tangent2D = rule;
    else if(option == "--gk-regular-3d") regular3D = rule;
    else if(option == "--gk-tangent-3d") tangent3D = rule;
    else throw std::invalid_argument("Unknown outer quadrature option: "+option);
  }

  static OuterQuadratureRules ParseArguments(int argc, char** argv) {
    OuterQuadratureRules result;
    for(int i = 1; i < argc; ++i) {
      const std::string option = argv[i];
      if(option.compare(0,5,"--gk-") != 0) continue;
      if(!GaussKronrodRules::IsOption(option)) throw std::invalid_argument("Unknown outer quadrature option: "+option);
      if(i+1 >= argc) throw std::invalid_argument("Missing value for "+option);
      result.SetOption(option,argv[++i]);
    }
    return result;
  }
};

}

#endif
