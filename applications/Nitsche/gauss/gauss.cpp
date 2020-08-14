/** tutorial/Ex1
 * This example shows how to:
 * initialize a femus application;
 * define the multilevel-mesh object mlMsh;
 * read from the file ./input/square.neu the coarse-level mesh and associate it to mlMsh;
 * add in mlMsh uniform refined level-meshes;
 * define the multilevel-solution object mlSol associated to mlMsh;
 * add in mlSol different types of finite element solution variables;
 * initialize the solution varables;
 * define vtk and gmv writer objects associated to mlSol;
 * print vtk and gmv binary-format files in ./output directory.
 **/

#include "GaussPoints.hpp"
#include "ElemType.hpp"
#include <cmath>
#include <iostream>
#include <boost/math/special_functions/gamma.hpp>

#include "gauss.hpp"

using namespace femus;

enum ORDER {
  ZERO_ORDER = 0, FIRST_ORDER, SECOND_ORDER, THIRD_ORDER, FOURTH_ORDER, FIFTH_ORDER,
  SIXTH_ORDER, SEVENTH_ORDER, EIGHTH_ORDER, NINTH_ORDER, TENTH_ORDER,
  ELEVENTH_ORDER, TWELFTH_ORDER, THIRTEENTH_ORDER, FOURTEENTH_ORDER, FIFTEENTH_ORDER,
  SIXTEENTH_ORDER, SEVENTEENTH_ORDER, EIGHTEENTH_ORDER, NINETEENTH_ORDER, TWENTIETH_ORDER,
  TWENTY_FIRST_ORDER, TWENTY_SECOND_ORDER, TWENTY_THIRD_ORDER, TWENTY_FOURTH_ORDER, TWENTY_FIFTH_ORDER,
  TWENTY_SIXTH_ORDER, TWENTY_SEVENTH_ORDER, TWENTY_EIGHTH_ORDER, TWENTY_NINTH_ORDER, THIRTIETH_ORDER,
  THIRTY_FIRST_ORDER, THIRTY_SECOND_ORDER, THIRTY_THIRD_ORDER, THIRTY_FOURTH_ORDER, THIRTY_FIFTH_ORDER,
  THIRTY_SIXTH_ORDER, THIRTY_SEVENTH_ORDER, THIRTY_EIGHTH_ORDER, THIRTY_NINTH_ORDER
};

void CppPrint(GeomElType &geomElType, const std::vector < std::vector < std::vector<double> > > & gaussQ);
void GnuPrint(GeomElType &geomElType, const std::vector < std::vector < std::vector<double> > > & gaussQ);
double GetIntegral(const std::vector<unsigned> &n, const std::vector < std::vector<double> > & gaussQ);

void TestLineIntegral(const unsigned &m, const std::vector < std::vector < std::vector<double> > > & gaussQ);
void TestQuadIntegral(const unsigned &m, const std::vector < std::vector < std::vector<double> > > & gaussQ);
void TestTriIntegral(const unsigned &m, const std::vector < std::vector < std::vector<double> > > & gaussQ);
void TestHexIntegral(const unsigned &m, const std::vector < std::vector < std::vector<double> > > & gaussQ);
void TestWedgeIntegral(const unsigned &m, const std::vector < std::vector < std::vector<double> > > & gaussQ);
void TestTetIntegral(const unsigned &m, const std::vector < std::vector < std::vector<double> > > & gaussQ);

void BuildGaussPoints(GeomElType &geomElType, const unsigned &m, const elem_type *femAll[3], std::vector < std::vector < double> > & gaussQ);

int main(int argc, char** args) {

  GeomElType geomElType[6] = {HEX, TET, WEDGE, QUAD, TRI, LINE};

  const elem_type *femAll[3];
  femAll[0] = new const elem_type_1D("line", "linear", "zero");
  femAll[1] = new const elem_type_2D("quad", "linear", "zero");
  femAll[2] = new const elem_type_3D("hex", "linear", "zero");

  unsigned maxOrder = THIRTIETH_ORDER;
  
  for(unsigned k = 0; k < 6; k++) {

    std::vector < std::vector < std::vector < double> > > gaussQ(maxOrder + 1);
    
    for(unsigned m = 0; m <= maxOrder; m++) {
      BuildGaussPoints(geomElType[k], m, femAll,  gaussQ[m]);
    }

    CppPrint(geomElType[k], gaussQ);
    //GnuPrint(geomElType[k], gaussQ);

    for(unsigned m = 0; m <= maxOrder + 5; m++) {
      if(HEX == geomElType[k]) {
        TestHexIntegral(m, gaussQ);
      }
      else if(TET == geomElType[k])  {
        TestTetIntegral(m, gaussQ);
      }
      else if(WEDGE == geomElType[k]) {
        TestWedgeIntegral(m, gaussQ);
      }
      else if(QUAD == geomElType[k])  {
        TestQuadIntegral(m, gaussQ);
      }
      else if(TRI == geomElType[k]) {
        TestTriIntegral(m, gaussQ);
      }
      else if(LINE == geomElType[k]) {
        TestLineIntegral(m, gaussQ);
      }
    }
    std::cout << std::endl;
  }

  for(unsigned k = 0; k < 3; k++)
    delete femAll[k];
}

void CppPrint(GeomElType &geomElType, const std::vector < std::vector < std::vector<double> > > & gaussQ) {

  std::string name[6] = {"hex", "tet", "wedge", "quad", "tri", "line"};

  std::ostringstream ofilename;
  ofilename << name[geomElType] << "GaussPoints.cpp";

  std::string filename(ofilename.str());

  std::ofstream fout;
  fout.open(filename);

  unsigned gQsize = gaussQ.size();
  unsigned dim = gaussQ[0].size() - 1;


  fout << "/*=========================================================================" << std::endl << std::endl;

  fout << "Program: FEMUS" << std::endl;
  fout << "Module: Gauss" << std::endl;
  fout << "Authors: Eugenio Aulisa, Giorgio Bornia, Erdi Kara" << std::endl << std::endl;

  fout << "Copyright (c) FEMTTU" << std::endl;
  fout << "All rights reserved." << std::endl << std::endl;

  fout << "This software is distributed WITHOUT ANY WARRANTY; without even" << std::endl;
  fout << "the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR" << std::endl;
  fout << "PURPOSE.  See the above copyright notice for more information." << std::endl << std::endl;

  fout << "=========================================================================*/" << std::endl << std::endl;

  fout << "#include \"GaussPoints.hpp\"" << std::endl;
  fout << "#include <iostream>" << std::endl;
  fout << "#include <stdlib.h>" << std::endl;
  fout << "#include <string.h>" << std::endl << std::endl;
  fout << "namespace femus {" << std::endl;

  fout << "const unsigned " << name[geomElType] << "_gauss::GaussPoints[" << gQsize << "] = {";
  for(unsigned m = 0; m < gQsize; m++) {
    fout << gaussQ[m][0].size();
    if(m < gQsize - 1) fout << ", ";
  }
  fout << "};\n\n";

  fout << "const double * " << name[geomElType] << "_gauss::Gauss[" << gQsize << "] = {";
  for(unsigned m = 0; m < gQsize; m++) {
    fout << "Gauss" << m << "[0]";
    if(m < gQsize - 1) fout << ", ";
  }
  fout << "};\n\n";

  for(unsigned m = 0; m < gQsize; m++) {
    unsigned size = gaussQ[m][0].size();
    fout.precision(14);
    fout << "const double " << name[geomElType] << "_gauss::Gauss" << m << "[" << dim + 1 << "][" <<  size << "] = {{";
    for(unsigned i = 0; i < size; i++) {
      fout << gaussQ[m][0][i];
      if(i < size - 1) fout << ", ";
    }
    for(unsigned k = 1; k <= dim; k++) {
      fout << "},\n{";

      for(unsigned i = 0; i < size; i++) {
        fout << gaussQ[m][k][i];
        if(i < size - 1) fout << ", ";
      }
    }
    fout << "}};" << std::endl << std::endl;
  }
  fout << "}" << std::endl;
  fout.close();
}

void GnuPrint(GeomElType &geomElType, const std::vector < std::vector < std::vector<double> > > & gaussQ) {

  std::string name[6] = {"hex", "tet", "wedge", "quad", "tri", "line"};

  unsigned gQsize = gaussQ.size();
  unsigned dim = gaussQ[0].size() - 1;

  for(unsigned m = 0; m < gQsize; m++) {

    std::ostringstream ofilename;
    ofilename << name[geomElType] << m << ".txt";

    std::string filename(ofilename.str());

    std::ofstream fout;
    fout.open(filename);


    for(unsigned i = 0; i < gaussQ[m][0].size(); i++) {
      for(unsigned k = 1; k <= dim; k++) {
        fout << gaussQ[m][k][i] << " ";
      }
      fout << gaussQ[m][0][i] << std::endl;
    }
    fout.close();

  }
}

double GetIntegral(const std::vector<unsigned> &n, const std::vector < std::vector<double> > & gaussQ) {

  unsigned dim = gaussQ.size() - 1;
  double I = 0;

  for(unsigned i = 0; i < gaussQ[0].size(); i++) {
    double a = 1.;
    for(unsigned k = 1; k <= dim; k++) {
      a *= pow(gaussQ[k][i], n[k - 1]);
    }
    I += a * gaussQ[0][i];
  }
  return I;
}

void TestLineIntegral(const unsigned & order, const std::vector < std::vector < std::vector<double> > > & gaussQ) {

  std::cout << "LINE ";

  unsigned m = (order % 2 == 0) ? order : order - 1;
  std::vector < unsigned > n(1, m);

  double I = (1. + pow(-1., n[0])) / (n[0] + 1.);

  m = (order > gaussQ.size() - 1) ? gaussQ.size() - 1 : order;
  double I1 = GetIntegral(n, gaussQ[m]);

  std::cout << n[0] << " " << m << " " << I << " " << I1 << " " << (I - I1) / I << std::endl;

}


void TestTriIntegral(const unsigned & order, const std::vector < std::vector < std::vector<double> > > & gaussQ) {

  std::cout << "TRI ";

  std::vector < unsigned > n(2);
  n[0] = order / 2;
  n[1] = order - n[0];

  double I = boost::math::tgamma(1 + n[0]) * boost::math::tgamma(1 + n[1]) / boost::math::tgamma(3 + n [0] + n[1]);

  unsigned m = (order > gaussQ.size() - 1) ? gaussQ.size() - 1 : order;
  double I1 = GetIntegral(n, gaussQ[m]);
  std::cout << n[0] << " " << n[1] << " " << n[0] + n[1] << " " << m << " " << I << " " << I1 << " " << (I - I1) / I << std::endl;

}


void TestQuadIntegral(const unsigned & order, const std::vector < std::vector < std::vector<double> > > & gaussQ) {

  std::cout << "QUAD ";

  unsigned m = (order % 2 == 0) ? order : order - 1;
  std::vector < unsigned > n(2, m);

  double I = (1. + pow(-1., n[0])) / (n[0] + 1.) * (1. + pow(-1., n[1])) / (n[1] + 1.);

  m = (order > gaussQ.size() - 1) ? gaussQ.size() - 1 : order;
  double I1 = GetIntegral(n, gaussQ[m]);

  std::cout << n[0] << " " << n[1] << " " << n[0] + n[1] << " " << m << " " << I << " " << I1 << " " << (I - I1) / I << std::endl;

}


void TestHexIntegral(const unsigned & order, const std::vector < std::vector < std::vector<double> > > & gaussQ) {

  std::cout << "HEX ";

  unsigned m = (order % 2 == 0) ? order : order - 1;

  std::vector < unsigned > n(3, m);

  double I = (1. + pow(-1., n[0])) / (n[0] + 1.) * (1. + pow(-1., n[1])) / (n[1] + 1) * (1. + pow(-1., n[2])) / (n[2] + 1) ;

  m = (order > gaussQ.size() - 1) ? gaussQ.size() - 1 : order;
  double I1 = GetIntegral(n, gaussQ[m]);

  std::cout << n[0] << " " << n[1] <<  " " << n[2] << " " << n[0] + n[1] + n[2] << " " << m << " " << I << " " << I1 << " " << (I - I1) / I << std::endl;

}

void TestWedgeIntegral(const unsigned & order, const std::vector < std::vector < std::vector<double> > > & gaussQ) {

  std::cout << "WEDGE ";

  std::vector < unsigned > n(3);

  n[0] = order / 2;
  n[1] = order - n[0];
  n[2] = (order % 2 == 0) ? order : order - 1;

  double I = boost::math::tgamma(1 + n[0]) * boost::math::tgamma(1 + n[1]) / boost::math::tgamma(3 + n [0] + n[1]) * (1. + pow(-1., n[2])) / (n[2] + 1) ;

  unsigned m = (order > gaussQ.size() - 1) ? gaussQ.size() - 1 : order;
  double I1 = GetIntegral(n, gaussQ[m]);

  std::cout << n[0] << " " << n[1] <<  " " << n[2] << " " << n[0] + n[1] << " " << m << " " << I << " " << I1 << " " << (I - I1) / I << std::endl;

}

void TestTetIntegral(const unsigned & order, const std::vector < std::vector < std::vector<double> > > & gaussQ) {

  std::cout << "TET ";

  std::vector < unsigned > n(3);
  n[0] = order / 3;
  n[1] = (order - n[0]) / 2;
  n[2] = order - n[0] - n[1];

  double I = boost::math::tgamma(1 + n[0]) * boost::math::tgamma(1 + n[1]) * boost::math::tgamma(1 + n[2]) / boost::math::tgamma(4 + n [0] + n[1] + n[2]);

  unsigned m = (order > gaussQ.size() - 1) ? gaussQ.size() - 1 : order;
  double I1 = GetIntegral(n, gaussQ[m]);

  std::cout << n[0] << " " << n[1] << " " << n[2] << " " << n[0] + n[1] + n[2] << " " << m << " " << I << " " << I1 << " " << (I - I1) / I << std::endl;

}

void BuildGaussPoints(GeomElType &geomElType, const unsigned &m, const elem_type *femAll[3],
                      std::vector < std::vector < double> > &gaussQ) {

  std::vector < std::vector < double > > xv;
  unsigned nv;
  unsigned dim;
  std::vector < unsigned > dm;
  unsigned m1 = 0;

  if(HEX == geomElType) {
    nv = 8;
    dim = 3;
    dm.assign(dim, 0);
    xv = {{ -1., 1., 1., -1., -1., 1., 1., -1.},
      { -1., -1., 1., 1., -1., -1., 1., 1.},
      { -1., -1., -1., -1., 1., 1., 1., 1.}
    }; //hex to hex mapping
  }
  else if(TET == geomElType) {
    nv = 8;
    dim = 3;
    dm.assign(dim, 2);
    m1 = 9;

    xv = {{0, 1, 0.5, 0., 0., 0.5, 1. / 3., 0},
      {0, 0, 0.5, 1., 0., 0., 1. / 3., 0.5},
      {0, 0, 0, 0, 1, 0.5, 1. / 3., 0.5}
    }; //hex to tet mapping
  }
  else if(WEDGE == geomElType) {
    nv = 8;
    dim = 3;
    dm.assign(dim, 1);
    dm[2] = 0;
    m1 = 14;
    xv = {{0., 1., 0.5, 0., 0., 1., 0.5, 0.},
      {0., 0., 0.5, 1., 0., 0., 0.5, 1.},
      { -1., -1., -1., -1., 1., 1., 1., 1.}
    }; //hex to wedge mapping
  }
  else if(QUAD == geomElType) {
    nv = 4;
    dim = 2;
    dm.assign(dim, 0);
    xv = {{ -1., 1., 1., -1.}, { -1., -1., 1., 1.}}; //quad to quad mapping
  }
  else if(TRI == geomElType) {
    nv = 4;
    dim = 2;
    dm.assign(dim, 1);
    m1 = 14;
    xv = {{0., 1., 0.5, 0.}, {0., 0., 0.5, 1.}}; //quad to tri mapping
  }
  else if(LINE == geomElType) {
    nv = 2;
    dim = 1;
    dm.assign(dim, 0);
    m1 = 40;
    xv = {{ -1., 1.}}; //line to line mapping
  }
  else {
    abort();
  }

  gaussQ.resize(dim + 1);

  if(m < m1) {
    if(geomElType == TET) {
      gaussQ = tetGauss[m];
    }
    else if(geomElType == WEDGE) {
      unsigned sizeT = triGauss[m][0].size();
      unsigned sizeL = lineGauss[m][0].size();
      unsigned size = sizeT * sizeL;

      for(unsigned k = 0; k <= dim; k++) {
        gaussQ[k].resize(size);
      }

      unsigned cnt = 0;
      for(unsigned i = 0; i < sizeT ; i++) {
        for(unsigned j = 0; j < sizeL ; j++) {
          gaussQ[0][cnt] = triGauss[m][0][i] * lineGauss[m][0][j];
          gaussQ[1][cnt] = triGauss[m][1][i];
          gaussQ[2][cnt] = triGauss[m][2][i];
          gaussQ[3][cnt] = lineGauss[m][1][j];
          cnt++;
        }
      }
    }
    else if(geomElType == TRI) {
      gaussQ = triGauss[m];
    }
    else if(geomElType == LINE) {
      gaussQ = lineGauss[m];
    }
  }
  else {
    std::vector <double> phi;  // local test function
    std::vector <double> phi_x; // local test function first order partial derivatives
    double jac;

    const elem_type *fem = femAll[dim - 1];

    unsigned size = 1;
    std::vector < unsigned > ng(dim);
    for(unsigned k = 0; k < dim; k++) {
      ng[k] = (m + dm[k]) / 2 + 1;
      size *= ng[k];
    }

    for(unsigned k = 0; k <= dim; k++) {
      gaussQ[k].resize(size);
    }
    std::vector < double > xi(dim);
    std::vector <unsigned> I(dim);
    std::vector <unsigned> NG(dim);
    NG[dim - 1] = 1;
    for(unsigned k = dim - 1 ; k > 0;  k--) {
      NG[k - 1] = NG[k] * ng[k];
    }

    for(unsigned cnt = 0; cnt < size ; cnt++) {
      I[0] = cnt / NG[0];
      for(unsigned k = 1; k < dim ; k++) {
        unsigned pk = cnt % NG[k - 1];
        I[k] = pk / NG[k];
      }
      gaussQ[0][cnt] = 1.;
      for(unsigned k = 0; k < dim; k++) {
        gaussQ[0][cnt] *= lineGauss[m + dm[k]][0][I[k]];
        xi[k] = lineGauss[m + dm[k]][1][I[k]];
      }
      fem->Jacobian(xv, xi, jac, phi, phi_x);
      gaussQ[0][cnt] *= jac;
      for(unsigned k = 0; k < dim; k++) {
        gaussQ[k + 1][cnt] = 0.;
        for(unsigned ii = 0; ii < nv; ii++) {
          gaussQ[k + 1][cnt] += phi[ii] * xv[k][ii];
        }
      }
    }
  }
  //delete fem;
}


