#pragma once

#include "GeomElTypeEnum.hpp"

#include "GaussPoints.hpp"

#include "parabolaIntegration.hpp"

#include <iostream>
#include <boost/multiprecision/cpp_bin_float.hpp>


using namespace femus;

template <class TypeA>
class polyWPar {
  public:
    // virtual ~polyWPar() {};
    // virtual void GetWeight(const std::vector<double> &A, std::vector<double> &weight) = 0;
};


template <class TypeA>
class polyWParQUAD :
  public polyWPar <TypeA> {
  public:
    polyWParQUAD(const unsigned &qM);
    polyWParQUAD(const unsigned &qM, int &maxDepth, double &maxRelErr);
    ~polyWParQUAD();

    void build();

    void setTableDepth(int depth) {_maxDepth = depth;}
    void setTableRelErr(double relErr) {_maxRelErr = relErr;}

 private:
     unsigned _dim;
     GeomElType _geomElemType;
     unsigned _qM;
     std::string _gaussType;

     int _s;
     PointT <TypeA> _p1, _p2, _p3;
     std::vector<double>_weightCF;

     int _maxDepth = 5; // max depth of the table (initialized to 5)
     double _maxRelErr = 0.001; // max relative error of the table (initialized to 10^-3)

     std::vector<OctreeNode<TypeA>>_loadedRoots; // Map

 };




// template <class Type>
//  class polyWPar {
//     public:
//       polyWPar(const GeomElType &geomElemType, const unsigned &qM, const std::string &gaussType);
//       ~polyWPar();
//
//       void build();
//
//     void setTableDepth(int depth) {_maxDepth = depth;}
//     void setTableRelErr(double relErr) {_maxRelErr = relErr;}
//     void prova();
//
//  private:
//      unsigned _dim;
//      GeomElType _geomElemType;
//      unsigned _qM;
//      std::string _gaussType;
//
//      int _s;
//      PointT <Type> _p1, _p2, _p3;
//      std::vector<double>_weightCF;
//
//      int _maxDepth = 5; // max depth of the table (initialized to 5)
//      double _maxRelErr = 0.001; // max relative error of the table (initialized to 10^-3)
//
//      std::vector<OctreeNode<Type>>_loadedRoots; // Map
//
//  };

