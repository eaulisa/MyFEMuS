#pragma once

#include "polyWPar.hpp"

#include "GeomElTypeEnum.hpp"

#include "GaussPoints.hpp"

#include "parabolaIntegration.hpp"

#include <iostream>
#include <boost/multiprecision/cpp_bin_float.hpp>


using namespace femus;

  template <class Type> polyWPar<Type>::polyWPar(const GeomElType &geomElemType, const unsigned &qM, const std::string &gaussType) {
    _geomElemType = geomElemType;
    _qM = qM;
    _gaussType = gaussType;

    if(_geomElemType == HEX || _geomElemType == WEDGE || _geomElemType == TET) _dim = 3;
    if(_geomElemType == QUAD || _geomElemType == TRI) _dim = 2;
    if(_geomElemType == LINE) _dim = 1;
    this->build();
  }


