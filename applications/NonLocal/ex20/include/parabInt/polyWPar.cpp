#pragma once

#include "polyWPar.hpp"

#include "GeomElTypeEnum.hpp"

#include "GaussPoints.hpp"

#include "parabolaIntegration.hpp"

#include <iostream>
#include <boost/multiprecision/cpp_bin_float.hpp>


using namespace femus;

template <class TypeA>
polyWParQUAD<TypeA>::polyWParQUAD(const unsigned &qM) {
  _qM = qM;
  _dim = 2;

  this->build();
}

template <class TypeA>
polyWParQUAD<TypeA>::polyWParQUAD(const unsigned &qM, const int &maxDepth, const double &maxRelErr) {
  _qM = qM;
  _dim = 2;
  _maxDepth = maxDepth;
  _maxRelErr = maxRelErr;
  this->build();
}



template <class TypeA>
void polyWParQUAD<TypeA>::build() {
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
void polyWParQUAD<TypeA>:: GetWeight(const std::vector<std::vector<double>> &xv, const std::vector<double> &A, const std::vector < std::vector < std::vector <double > > > &aP, std::vector<double> &weightCF, bool &twoInt) {
  twoInt = find_Weight_CF(xv, A, aP, weightCF);
}



template <class TypeA>
bool polyWParQUAD<TypeA>:: find_Weight_CF(const std::vector<std::vector<double>> &xv, const std::vector<double> &A, const std::vector < std::vector < std::vector <double > > > &aP, std::vector<double> &modified_weights) {

  unsigned nInt;
  std::vector<std::vector<double>> unitxv = {{0., 1., 1., 0.}, {0., 0., 1., 1.}};
  unsigned nPoints = 3;
  short unsigned ielType = 3; //quad
  unsigned femType = 0; //linear FEM
  std::vector< double > interp_point_weights;
  PointT <TypeA> p1, p2, p3;
  unsigned table_number;
  PointT <double> q1, q2, q3;
  Point3D searchP(0., 0., 0.);

  GetCellPointsFromQuadric(xv, A, nPoints, nInt);     //This finds the intersection points + the middle point in physical space

  if(nInt != 2) {
    return false;
  }
  else {

    std::vector<int> xvsign(4);
    std::vector<int> unitxvsign(4);

    std::vector<std::vector<double>> xi(nPoints, std::vector<double>(2, 0.));
    for(unsigned i = 0; i < nPoints; i++) {
      bool inverseMapping = GetInverseMapping(femType, ielType, aP, _xe[i], xi[i], 100);        //This maps the phsical points to {(-1,-1),(1,1)} box
      xi[i] = {0.5 * (xi[i][0] + 1.), 0.5 * (xi[i][1] + 1.)};                                        // //This maps the points to unit box
    }

    for(unsigned l = 0; l < 4 ; l++) {
      xvsign[l] = ((A[0] * xv[0][l] * xv[0][l] + A[1] * xv[0][l] * xv[1][l] + A[2] * xv[1][l] * xv[1][l] + A[3] * xv[0][l] + A[4] * xv[1][l] + A[5]) >= 0) ? 1 : -1 ;
    }

    // bool vertical = false;
    // if(fabs(xi[0][0] - xi[2][0]) >= fabs(xi[0][1] - xi[2][1])) {
    //   if((xi[0][0] < xi[1][0] && xi[1][0] < xi[2][0]) || (xi[0][0] > xi[1][0] && xi[1][0] > xi[2][0])) vertical = true;
    //   else vertical = false;
    // }
    // else {
    //   if((xi[0][1] < xi[1][1] && xi[1][1] < xi[2][1]) || (xi[0][1] > xi[1][1] && xi[1][1] > xi[2][1])) vertical = false;
    //   else vertical = true;
    // }

    bool vertical = false;
    bool xSpan = false;
    bool ySpan = false;
    if((xi[0][0] < xi[1][0] && xi[1][0] < xi[2][0]) || (xi[0][0] > xi[1][0] && xi[1][0] > xi[2][0])) xSpan = true ;
    if((xi[0][1] < xi[1][1] && xi[1][1] < xi[2][1]) || (xi[0][1] > xi[1][1] && xi[1][1] > xi[2][1])) ySpan = true ;

//     if(xSpan && ySpan){
//       if(fabs(xi[0][0] - xi[2][0]) >= fabs(xi[0][1] - xi[2][1])) vertical = true;
//       else vertical = false ;
//     }
//     else if(xSpan && !ySpan) vertical = true;
//     else if(!xSpan && ySpan) vertical = false;
//     else std::cout << " The parabola formed by ths three points is not a function. Use line cut " << std::endl;


    if(xSpan) {
      if(ySpan) {
        if(fabs(xi[0][0] - xi[2][0]) >= fabs(xi[0][1] - xi[2][1])) vertical = true;
        else vertical = false;
      }
      else {
        vertical = true;
      }
    }
    else {
      if(ySpan) vertical = false;
      else {
        std::cout << " The parabola formed by this three points is not a function. Use line cut " << std::endl;
        return false;
      }
    }

    if(vertical) {

      q1 = { xi[0][0], xi[0][1] };
      q2 = { xi[2][0], xi[2][1] };
      q3 = { xi[1][0], xi[1][1] };

      Parabola <double> parabola = get_parabola_equation(q1, q2, q3);
      int normal;

      for(unsigned l = 0; l < 4 ; l++) {
        unitxvsign[l] = ((parabola.k * unitxv[0][l] * unitxv[0][l] + parabola.b * unitxv[0][l] + parabola.d + unitxv[1][l]) > 0) ? 1 : -1;
      }
      normal = checkVectorRelation(xvsign, unitxvsign);

      q3.x = (q1.x + q2.x) / 2;
      q3.y = -parabola.k * q3.x * q3.x - parabola.b * q3.x - parabola.d ;
      find_search_table(q1, q2, q3, table_number, searchP);
      OctreeNode<cpp_bin_float_oct>* result = _loadedRoots[table_number].search(searchP);

      if(result) {

        std::vector<double>interp_point = {searchP.x, searchP.y, searchP.z};
        trilinier_interpolation_vector(result->corners, result->cornerWeights, interp_point, interp_point_weights);
        modified_weights.resize(interp_point_weights.size());

        if(normal == -1) {
          for(unsigned aq = 0; aq < interp_point_weights.size(); aq++) {
            modified_weights[aq] = 1 - interp_point_weights[aq];
          }
        }
        else modified_weights = interp_point_weights;
      }
      else {
        std::cout << "Search point not found in the Octree." << std::endl;
      }
    }

    else {//horizontal

      q1 = { xi[0][1], xi[0][0] };
      q2 = { xi[2][1], xi[2][0] };
      q3 = { xi[1][1], xi[1][0] };

      Parabola <double> parabola = get_parabola_equation(q1, q2, q3);
      int normal;

      //use horizotal parabola for the normal
      for(unsigned l = 0; l < 4 ; l++) {
        unitxvsign[l] = ((static_cast<double>(parabola.k) * unitxv[1][l] * unitxv[1][l] + static_cast<double>(parabola.b) * unitxv[1][l] + static_cast<double>(parabola.d) + unitxv[0][l]) > 0) ? 1 : -1;
      }


      normal = checkVectorRelation(xvsign, unitxvsign);
      q3.x = (q1.x + q2.x) / 2.;
      q3.y = -parabola.k * q3.x * q3.x - parabola.b * q3.x - parabola.d ;
      find_search_table(q1, q2, q3, table_number, searchP);

      OctreeNode<cpp_bin_float_oct>* result = _loadedRoots[table_number].search(searchP);
      if(result) {

        std::vector<double>interp_point = {searchP.x, searchP.y, searchP.z};
        trilinier_interpolation_vector(result->corners, result->cornerWeights, interp_point, interp_point_weights);
        modified_weights.resize(interp_point_weights.size());

        if(normal == -1) {
          int sqrt_size = sqrt(interp_point_weights.size());
          for(unsigned ai = 0; ai < sqrt_size; ai++) {
            for(unsigned aj = 0; aj < sqrt_size; aj++) {
              modified_weights[ai * sqrt_size + aj] = 1 - interp_point_weights[aj * sqrt_size + ai];
            }
          }
        }
        else {
          int sqrt_size = sqrt(interp_point_weights.size());
          for(unsigned ai = 0; ai < sqrt_size; ai++) {
            for(unsigned aj = 0; aj < sqrt_size; aj++) {
              modified_weights[ai * sqrt_size + aj] = interp_point_weights[aj * sqrt_size + ai];
            }
          }
        }


      }
      else {
        std::cout << "Search point not found in the Octree." << std::endl;
      }

    }
  }

  return true ;

}


template <class TypeA> void polyWParQUAD<TypeA>:: GetCellPointsFromQuadric(const std::vector<std::vector<double>> &xv, const std::vector<double> &Cf, unsigned npt, unsigned & nInt/*, unsigned level*/) {

  //typedef cpp_bin_float_oct myType;
  typedef double myType;

  unsigned cnt = 0;
  _xe.resize(((8 < npt) ? npt : 8), std::vector<double>(_dim));       //xe is a matrix with size (min 8 by 2) TODO is this the other way around?
  _ds.resize(npt);     //ds has size number of marker.

  const unsigned nve = xv[0].size();     //nve = 2 (x,y)
  std::vector<double> v(_dim, 0.); // zero vector with size 4

  double pm[2] = {1, -1};

  for(unsigned i = 0; i < nve; i++) {
    unsigned ip1 = (i + 1) % nve;
    for(unsigned k = 0; k < _dim; k++) v[k] = xv[k][ip1] - xv[k][i];   // (xi-yi)

    const double &x0 = xv[0][i];     //isn't it more like x0 and x1
    const double &y0 = xv[1][i];

    myType a = Cf[0] * v[0] * v[0] + Cf[1] * v[0] * v[1] + Cf[2] * v[1] * v[1];
    myType b = 2 * Cf[0] * v[0] * x0 + Cf[1] * v[1] * x0 + Cf[1] * v[0] * y0 + 2 * Cf[2] * v[1] * y0 + Cf[3] * v[0] + Cf[4] * v[1];
    myType c = Cf[0] * x0 * x0 + Cf[1] * x0 * y0 + Cf[2] * y0 * y0 + Cf[3] * x0 + Cf[4] * y0 + Cf[5];

    // myType norm = sqrt(a * a + b * b + c * c);
    // a /= norm;
    // b /= norm;
    // c /= norm;
    //if(fabs(a) > 1.e-8) {

    myType norm2 = a * a + b * b + c * c;
    if(a * a / norm2 > 1.e-16) {
      myType delta = b * b - 4 * a * c;
      if(delta > 0) {
        myType sqrtDelta = sqrt(delta);
        for(unsigned j = 0; j < 2; j++) {
          double t = static_cast<double>((- b + pm[j] * sqrtDelta) / (2 * a));
          if(t >= 0 && t <= 1) {
            for(unsigned  k = 0; k < _dim; k++) {
              _xe[cnt][k] = xv[k][i]  + t * v[k];
            }
            cnt++;
          }
        }
      }
    }
    else if(b != 0) {
      double t = static_cast<double>(-c / b);
      if(t >= 0 && t <= 1) {
        for(unsigned  k = 0; k < _dim; k++) {
          _xe[cnt][k] = xv[k][i]  + t * v[k];
        }
        cnt++;
      }
    }
  }

  nInt = cnt;

  if(cnt == 2) {

    std::vector<double> Xg(2, 0);


    for(unsigned i = 0; i < nve; i++) {
      Xg[0] += xv[0][i];
      Xg[1] += xv[1][i];
    }
    Xg = {Xg[0] / nve, Xg[1] / nve};

    std::vector<double> P1 = _xe[0];
    std::vector<double> P2 = _xe[1];

    BuildMarkersOnConicArc(2 * M_PI, npt, Cf, Xg, P1, P2, _xe);

    // npt = _xe.size();

    // _ds.assign(npt, 0);
    // for(unsigned i = 0; i < _xe.size() - 1; i++) {
    //   double DS = 0.5 * sqrt((_xe[i][0] - _xe[i + 1][0]) * (_xe[i][0] - _xe[i + 1][0]) + (_xe[i][1] - _xe[i + 1][1]) * (_xe[i][1] - _xe[i + 1][1]));
    //   _ds[i] += DS;
    //   _ds[i + 1] += DS;
    // }

  }
}
