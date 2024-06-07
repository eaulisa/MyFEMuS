#pragma once

#include "GeomElTypeEnum.hpp"

#include "GaussPoints.hpp"

#include "parabolaIntegration.hpp"

#include <iostream>
#include <boost/multiprecision/cpp_bin_float.hpp>


using namespace femus;

template <class Type>
 class polyWPar {
    public:
      polyWPar(const GeomElType &geomElemType, const unsigned &qM, const std::string &gaussType);
      ~polyWPar();

      void build() {
      // TODO this is a generic random initialization of Pweights, needs to be done better
      _s = 0; // TODO in this case we use s = 0 for area integration
      _p1 = { static_cast<Type>(0), static_cast<Type>(0.5) }; // TODO all these variables can be eliminated in the future
      _p2 = { static_cast<Type>(0.5), static_cast<Type>(1) };
      _p3 = { static_cast<Type>((_p1.x + _p2.x) / 2.0), static_cast<Type>(0.125) };

      CutFemWeightParabola <double, Type> Pweights(_geomElemType, _qM, "legendre");
      Pweights(_s, 0, 1, 0, _p1, _p2, _p3, _weightCF); // TODO here we put a=0, c=1, table=0

      generateAndLoadOctrees<Type>(_maxDepth, _qM, _maxRelErr, Pweights, _loadedRoots);
    };

    void setTableDepth(int depth) {_maxDepth = depth;}
    void setTableRelErr(double relErr) {_maxRelErr = relErr;}

 private:
     unsigned _dim;
     GeomElType _geomElemType;
     unsigned _qM;
     std::string _gaussType;

     int _s;
     PointT <Type> _p1, _p2, _p3;
     std::vector<double>_weightCF;

     int _maxDepth = 5; // max depth of the table (initialized to 5)
     double _maxRelErr = 0.001; // max relative error of the table (initialized to 10^-3)

     std::vector<OctreeNode<Type>>_loadedRoots; // Map

 };

