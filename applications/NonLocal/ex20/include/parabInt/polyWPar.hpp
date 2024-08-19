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
    virtual ~polyWPar() {};
    virtual void GetWeight(const std::vector<std::vector<double>> &xv, const std::vector<double> &A, const std::vector < std::vector < std::vector <double > > > &aP, std::vector<double> &weightCF, bool &twoInt) = 0;
};


template <class TypeA>
class polyWParQUAD :
  public polyWPar <TypeA> {
  public:
    polyWParQUAD(const unsigned &qM);
    polyWParQUAD(const unsigned &qM, const int &maxDepth, const double &maxRelErr);
    // ~polyWParQUAD();

    void build();

    void setTableDepth(int depth) {_maxDepth = depth;}
    void setTableRelErr(double relErr) {_maxRelErr = relErr;}

    void GetWeight(const std::vector<std::vector<double>> &xv, const std::vector<double> &A, const std::vector < std::vector < std::vector <double > > > &aP, std::vector<double> &weightCF, bool &twoInt);
    bool find_Weight_CF( const std::vector<std::vector<double>> &xv, const std::vector<double> &A, const std::vector < std::vector < std::vector <double > > > &aP, std::vector<double> &modified_weights);
    void GetCellPointsFromQuadric (const std::vector<std::vector<double>> &xv, const std::vector<double> &Cf, unsigned npt, unsigned & nInt/*, unsigned level*/);

 private:
     unsigned _dim;
     GeomElType _geomElemType;
     unsigned _qM;
     std::string _gaussType;

     int _s;
     PointT <cpp_bin_float_oct> _p1, _p2, _p3;
     std::vector<double>_weightCF;

     int _maxDepth = 6; // max depth of the table (initialized to 5)
     double _maxRelErr = 0.0005; // max relative error of the table (initialized to 10^-3)

     std::vector<OctreeNode<cpp_bin_float_oct>>_loadedRoots; // Map

     std::vector < std::vector <double> > _xe;
     std::vector <double> _ds;

 };





