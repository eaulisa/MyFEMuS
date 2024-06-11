#pragma once

#include "polyWPar.hpp"

#include "GeomElTypeEnum.hpp"

#include "GaussPoints.hpp"

#include "parabolaIntegration.hpp"

#include <iostream>
#include <boost/multiprecision/cpp_bin_float.hpp>


using namespace femus;

  template <class TypeA>
  polyWParQUAD<TypeA>::polyWParQUAD (const unsigned &qM) {
    _qM = qM;
    _dim = 2;

    this->build();
  }

  template <class TypeA>
  polyWParQUAD<TypeA>::polyWParQUAD (const unsigned &qM, const int &maxDepth, const double &maxRelErr) {
    _qM = qM;
    _dim = 2;
    _maxDepth = maxDepth;
    _maxRelErr = maxRelErr;
    this->build();
  }



  template <class TypeA>
  void polyWParQUAD<TypeA>::build(){
      // TODO this is a generic random initialization of Pweights, needs to be done better
      _s = 0; // TODO in this case we use s = 0 for area integration
      _p1 = { static_cast<cpp_bin_float_oct>(0), static_cast<cpp_bin_float_oct>(0.5) }; // TODO all these variables can be eliminated in the future
      _p2 = { static_cast<cpp_bin_float_oct>(0.5), static_cast<cpp_bin_float_oct>(1) };
      _p3 = { static_cast<cpp_bin_float_oct>((_p1.x + _p2.x) / 2.0), static_cast<cpp_bin_float_oct>(0.125) };

      CutFemWeightParabola <double, cpp_bin_float_oct> Pweights(QUAD, _qM, "legendre");
      Pweights(_s, 0, 1, 0, _p1, _p2, _p3, _weightCF); // TODO here we put a=0, c=1, table=0

      generateAndLoadOctrees<cpp_bin_float_oct>(_maxDepth, _qM, _maxRelErr, Pweights, _loadedRoots);






    }

  template <class TypeA>
  void polyWParQUAD<TypeA>:: GetWeight (const std::vector<std::vector<double>> &xv, const std::vector<double> &A, std::vector<double> &weightCF, bool &twoInt) {
        twoInt = find_Weight_CF(_loadedRoots, xv, A, weightCF);
  }
