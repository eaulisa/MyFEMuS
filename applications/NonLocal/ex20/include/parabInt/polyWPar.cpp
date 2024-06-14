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
        twoInt = find_Weight_CF(xv, A, weightCF);
  }



  template <class TypeA>
  bool polyWParQUAD<TypeA>:: find_Weight_CF( const std::vector<std::vector<double>> &xv, const std::vector<double> &A, std::vector<double> &modified_weights){

    unsigned nInt;
    std::vector<std::vector<double>> unitxv = {{0., 1., 1., 0.}, {0., 0., 1., 1.}};
    unsigned nPoints = 3;
    short unsigned ielType = 3; //quad
    unsigned femType = 0; //linear FEM
    std::vector< double > interp_point_weights;
    PointT <TypeA> p1, p2, p3;
    bool twoInt = true;

    Fem fem = Fem(3 * 2, _dim);
    unsigned quad = 3;
    unsigned linear = 0;
    const elem_type *femQuad = fem.GetFiniteElement(quad, linear);

    GetCellPointsFromQuadric(xv, A, nPoints, nInt);     //This fins the points in physical space

    if(nInt != 2){return false;}
    else{

    std::vector < std::vector < std::vector <double > > > aP(1);
    ProjectNodalToPolynomialCoefficients(aP[femType], xv, ielType, femType);

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


    if((xi[0][0] < xi[1][0] && xi[1][0] < xi[2][0]) || (xi[0][0] > xi[1][0] && xi[1][0] > xi[2][0])) {  //vertical


      p1 = { static_cast<TypeA>(xi[0][0]), static_cast<TypeA>(xi[0][1]) };
      p2 = { static_cast<TypeA>(xi[2][0]), static_cast<TypeA>(xi[2][1]) };
      p3 = { static_cast<TypeA>(xi[1][0]), static_cast<TypeA>(xi[1][1]) };

      Parabola <TypeA> parabola = get_parabola_equation(p1, p2, p3);
      int normal;

      for(unsigned l = 0; l < 4 ; l++) {
        unitxvsign[l] = ((static_cast<double>(parabola.k) * unitxv[0][l] * unitxv[0][l] + static_cast<double>(parabola.b) * unitxv[0][l] + static_cast<double>(parabola.d) + unitxv[1][l]) > 0) ? 1 : -1;
      }

      normal = checkVectorRelation(xvsign, unitxvsign);

      int intersect_number;
      unsigned table_number;
      std::vector <TypeA> intersection;
      std::vector <TypeA> interp_point;
      CheckIntersection<TypeA>(intersect_number, table_number, intersection, interp_point, parabola);
      p3.x = (p1.x + p2.x) / 2;
      p3.y = -parabola.k * p3.x * p3.x - parabola.b * p3.x - parabola.d ;

      // if(interp_point.size() == 2) {
        Point3D searchP(static_cast<double>(interp_point[0]), static_cast<double>(interp_point[1]), static_cast<double>(p3.y));
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

          std::vector<double> phi, gradPhi;
          std::vector<double> Xg(femQuad->GetGaussPointNumber(),0);
          std::vector<double> Yg(femQuad->GetGaussPointNumber(),0);
          std::vector<double> Jg(femQuad->GetGaussPointNumber(),0);
          for(unsigned ig = 0; ig < femQuad->GetGaussPointNumber(); ig++) {
            // *** get gauss point weight, test function and test function partial derivatives ***
            femQuad->Jacobian(xv, ig, Jg[ig], phi, gradPhi);
            for(unsigned i =0;i<phi.size();i++){
              Xg[ig] += phi[i]*xv[0][i];
              Yg[ig] += phi[i]*xv[1][i];
            }
          }


        }
        else {
          std::cout << "Search point not found in the Octree." << std::endl;
        }
      // }
      // else{
      //   twoInt = false;
      //
      // }
    }



    else if((xi[0][1] < xi[1][1] && xi[1][1] < xi[2][1]) || (xi[0][1] > xi[1][1] && xi[1][1] > xi[2][1])) { //horizontal

      p1 = { static_cast<TypeA>(xi[0][1]), static_cast<TypeA>(xi[0][0]) };
      p2 = { static_cast<TypeA>(xi[2][1]), static_cast<TypeA>(xi[2][0]) };
      p3 = { static_cast<TypeA>(xi[1][1]), static_cast<TypeA>(xi[1][0]) };

      Parabola <TypeA> parabola = get_parabola_equation(p1, p2, p3);
      int normal;

      //use horizotal parabola for the normal
      for(unsigned l = 0; l < 4 ; l++) {
        unitxvsign[l] = ((static_cast<double>(parabola.k) * unitxv[1][l] * unitxv[1][l] + static_cast<double>(parabola.b) * unitxv[1][l] + static_cast<double>(parabola.d) + unitxv[0][l]) > 0) ? 1 : -1;
      }


      normal = checkVectorRelation(xvsign, unitxvsign);


      int intersect_number;
      unsigned table_number;
      std::vector <TypeA> intersection;
      std::vector <TypeA> interp_point;
      CheckIntersection<TypeA>(intersect_number, table_number, intersection, interp_point, parabola);

      p3.x = (p1.x + p2.x) / 2;
      p3.y = -parabola.k * p3.x * p3.x - parabola.b * p3.x - parabola.d ;

      // if(interp_point.size() == 2) {

        Point3D searchP(static_cast<double>(interp_point[0]), static_cast<double>(interp_point[1]), static_cast<double>(p3.y));
        OctreeNode<cpp_bin_float_oct>* result = _loadedRoots[table_number].search(searchP);
        if(result) {
          std::vector<double>interp_point = {searchP.x, searchP.y, searchP.z};

          trilinier_interpolation_vector(result->corners, result->cornerWeights, interp_point, interp_point_weights);

          modified_weights.resize(interp_point_weights.size());


          if (table_number == 2 || table_number == 4){
            if(normal == -1) {
              for(unsigned aq = 0; aq < interp_point_weights.size(); aq++) {
                modified_weights[aq] = 1 - interp_point_weights[interp_point_weights.size()-1-aq];
              }
            }
            else{
              for(unsigned aq = 0; aq < interp_point_weights.size(); aq++) {
                modified_weights[aq] = interp_point_weights[interp_point_weights.size()-1-aq];
              }
            }
          }

          else{
            if(normal == -1) {
              for(unsigned aq = 0; aq < interp_point_weights.size(); aq++) {
                modified_weights[aq] = 1 - interp_point_weights[aq];
              }
            }
            else{
              modified_weights = interp_point_weights;
            }
          }

          std::vector<double> phi, gradPhi;
          std::vector<double> Xg(femQuad->GetGaussPointNumber(),0);
          std::vector<double> Yg(femQuad->GetGaussPointNumber(),0);
          std::vector<double> Jg(femQuad->GetGaussPointNumber(),0);
          for(unsigned ig = 0; ig < femQuad->GetGaussPointNumber(); ig++) {
            // *** get gauss point weight, test function and test function partial derivatives ***
            femQuad->Jacobian(xv, ig, Jg[ig], phi, gradPhi);
            for(unsigned i =0;i<phi.size();i++){
              Xg[ig] += phi[i]*xv[0][i];
              Yg[ig] += phi[i]*xv[1][i];
            }
          }
        }
        else {
          std::cout << "Search point not found in the Octree." << std::endl;
        }
      // }
      // else{
      //   twoInt = false;
      // }

    }
    }

    return twoInt ;

  }


  template <class TypeA> void polyWParQUAD<TypeA>:: GetCellPointsFromQuadric (const std::vector<std::vector<double>> &xv, const std::vector<double> &Cf, unsigned npt, unsigned & nInt/*, unsigned level*/) {

  typedef cpp_bin_float_oct oct;

  unsigned cnt = 0;
  _xe.resize(((8 < npt) ? npt : 8), std::vector<double>(_dim));       //xe is a matrix with size (min 8 by 2) TODO is this the other way around?
  _ds.resize(npt);     //ds has size number of marker.

  const unsigned nve = xv[0].size();     //nve = 2 (x,y)
  std::vector<double> v(_dim, 0.); // zero vector with size 4

  for(unsigned i = 0; i < nve; i++) {
    unsigned ip1 = (i + 1) % nve;
    for(unsigned k = 0; k < _dim; k++) v[k] = xv[k][ip1] - xv[k][i];   // (xi-yi)

    const double &x0 = xv[0][i];     //isn't it more like x0 and x1
    const double &y0 = xv[1][i];

    oct a = Cf[0] * v[0] * v[0] + Cf[1] * v[0] * v[1] + Cf[2] * v[1] * v[1];
    oct b = 2 * Cf[0] * v[0] * x0 + Cf[1] * v[1] * x0 + Cf[1] * v[0] * y0 + 2 * Cf[2] * v[1] * y0 + Cf[3] * v[0] + Cf[4] * v[1];
    oct c = Cf[0] * x0 * x0 + Cf[1] * x0 * y0 + Cf[2] * y0 * y0 + Cf[3] * x0 + Cf[4] * y0 + Cf[5];

    oct norm = sqrt(a * a + b * b + c * c);
    a /= norm;
    b /= norm;
    c /= norm;

    if(fabs(a) > 1.e-8) {
      oct delta = b * b - 4 * a * c;
      if(delta > 0) {
        for(unsigned j = 0; j < 2; j++) {
          double t = static_cast<double>((- b + pow(-1, j) * sqrt(delta)) / (2 * a));
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

    BuildMarkersOnConicArc(2*M_PI, npt, Cf, Xg, P1, P2, _xe);

    // npt = _xe.size();

    // _ds.assign(npt, 0);
    // for(unsigned i = 0; i < _xe.size() - 1; i++) {
    //   double DS = 0.5 * sqrt((_xe[i][0] - _xe[i + 1][0]) * (_xe[i][0] - _xe[i + 1][0]) + (_xe[i][1] - _xe[i + 1][1]) * (_xe[i][1] - _xe[i + 1][1]));
    //   _ds[i] += DS;
    //   _ds[i + 1] += DS;
    // }

  }
}
