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

void CppPrint(GeomElType &geomElType, const std::vector < std::vector<double> > & weight, const std::vector < std::vector < std::vector<double> > > & x);
void GnuPrint(GeomElType &geomElType, const std::vector < std::vector<double> > & weight, const std::vector < std::vector < std::vector<double> > > & x);
double GetIntegral(const unsigned &m1, const std::vector<unsigned> &n, const std::vector < std::vector<double> > & weight, const std::vector < std::vector < std::vector<double> > > & x);
void TestQuadIntegral(const unsigned &m, const std::vector < std::vector<double> > & weight, const std::vector < std::vector < std::vector<double> > > & x);
void TestTriIntegral(const unsigned &m, const std::vector < std::vector<double> > & weight, const std::vector < std::vector < std::vector<double> > > & x);
void TestHexIntegral(const unsigned &m, const std::vector < std::vector<double> > & weight, const std::vector < std::vector < std::vector<double> > > & x);
void TestWedgeIntegral(const unsigned &m, const std::vector < std::vector<double> > & weight, const std::vector < std::vector < std::vector<double> > > & x);
void TestTetIntegral(const unsigned &m, const std::vector < std::vector<double> > & weight, const std::vector < std::vector < std::vector<double> > > & x);
void BuildGaussPoints(GeomElType &geomElType, const unsigned &maxNG,  std::vector < std::vector < double> > &weight,  std::vector < std::vector < std::vector < double> > > &x);

int main(int argc, char** args) {

  std::cout << tetGauss.size() << " "<< tetGauss[3].size() <<" " << tetGauss[3][0].size() <<std::endl;  
    
  GeomElType geomElType[6] = {HEX, TET, WEDGE, QUAD, TRI, LINE };

  ORDER order = ZERO_ORDER;
  std::cout << order / 2 << std::endl;

  order = FIRST_ORDER;
  std::cout << order / 2 << std::endl;

  order = SECOND_ORDER;
  std::cout << order / 2 << std::endl;


  order = TWENTY_SEVENTH_ORDER;
  std::cout << order / 2 << std::endl;

  order = TWENTY_EIGHTH_ORDER;
  std::cout << order / 2 << std::endl;

  order = TWENTY_NINTH_ORDER;
  std::cout << order / 2 << std::endl;

  //return 1;

  for(unsigned k = 0; k < 5; k++) {

    unsigned maxNG = 13;
    std::vector < std::vector < double> > weight(maxNG);
    std::vector < std::vector < std::vector < double> > > x(maxNG);

    BuildGaussPoints(geomElType[k], maxNG,  weight,  x);

    CppPrint(geomElType[k], weight, x);
    GnuPrint(geomElType[k], weight, x);

    for(unsigned m = 0; m < maxNG + 5; m++) {
      if(HEX == geomElType[k]) {
        TestHexIntegral(m, weight, x);
      }
      else if(TET == geomElType[k])  {
        TestTetIntegral(m, weight, x);
      }
      else if(WEDGE == geomElType[k]) {
        TestWedgeIntegral(m, weight, x);
      }
      else if(QUAD == geomElType[k])  {
        TestQuadIntegral(m, weight, x);
      }
      else if(TRI == geomElType[k]) {
        TestTriIntegral(m, weight, x);
      }
    }
    std::cout << std::endl;
  }
}

void CppPrint(GeomElType &geomElType, const std::vector < std::vector<double> > & weight, const std::vector < std::vector < std::vector<double> > > & x) {

  std::string name[6] = {"hex", "tet", "wedge", "quad", "tri", "line"};

  std::ostringstream ofilename;
  ofilename << name[geomElType] << "GaussPoints.cpp";

  std::string filename(ofilename.str());

  std::ofstream fout;
  fout.open(filename);

  unsigned maxNG = weight.size();
  unsigned dim = x[0].size();


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

  fout << "const unsigned " << name[geomElType] << "_gauss::GaussPoints[" << maxNG << "] = {";
  for(unsigned m = 0; m < maxNG; m++) {
    fout << weight[m].size();
    if(m < maxNG - 1) fout << ", ";
  }
  fout << "};\n\n";

  fout << "const double * " << name[geomElType] << "_gauss::Gauss[" << maxNG << "] = {";
  for(unsigned m = 0; m < maxNG; m++) {
    fout << "Gauss" << m << "[0]";
    if(m < maxNG - 1) fout << ", ";
  }
  fout << "};\n\n";

  for(unsigned m = 0; m < maxNG; m++) {
    fout.precision(14);
    fout << "const double " << name[geomElType] << "_gauss::Gauss" << m << "[" << dim + 1 << "][" <<  weight[m].size() << "] = {{";
    for(unsigned i = 0; i < weight[m].size(); i++) {
      fout << weight[m][i];
      if(i < weight[m].size() - 1) fout << ", ";
    }
    for(unsigned k = 0; k < dim; k++) {
      fout << "},\n{";

      for(unsigned i = 0; i < weight[m].size(); i++) {
        fout << x[m][k][i];
        if(i < weight[m].size() - 1) fout << ", ";
      }
    }
    fout << "}};" << std::endl << std::endl;
  }
  fout << "}" << std::endl;
  fout.close();
}

void GnuPrint(GeomElType &geomElType, const std::vector < std::vector<double> > & weight, const std::vector < std::vector < std::vector<double> > > & x) {

  std::string name[6] = {"hex", "tet", "wedge", "quad", "tri", "line"};

  unsigned maxNG = weight.size();
  unsigned dim = x[0].size();

  for(unsigned m = 0; m < maxNG; m++) {

    std::ostringstream ofilename;
    ofilename << name[geomElType] << m << ".txt";

    std::string filename(ofilename.str());

    std::ofstream fout;
    fout.open(filename);


    for(unsigned i = 0; i < weight[m].size(); i++) {
      for(unsigned k = 0; k < dim; k++) {
        fout << x[m][k][i] << " ";
      }
      fout << weight[m][i] << std::endl;
    }
    fout.close();

  }
}

double GetIntegral(const std::vector<unsigned> &n, const std::vector<double> & weight, const std::vector < std::vector<double> > & x) {

  unsigned dim = x.size();

  double I = 0;

  for(unsigned i = 0; i < weight.size(); i++) {
    double a = 1.;
    for(unsigned k = 0; k < dim; k++) {
      a *= pow(x[k][i], n[k]);
    }
    I += a * weight[i];
  }
  return I;
}


void TestTriIntegral(const unsigned & order, const std::vector < std::vector<double> > & weight, const std::vector < std::vector < std::vector<double> > > & x) {

  std::cout << "TRI ";

  std::vector < unsigned > n(2);
  n[0] = order / 2;
  n[1] = order - n[0];

  double I = boost::math::tgamma(1 + n[0]) * boost::math::tgamma(1 + n[1]) / boost::math::tgamma(3 + n [0] + n[1]);

  unsigned m = (order > weight.size() - 1) ? weight.size() - 1 : order;
  double I1 = GetIntegral(n, weight[m], x[m]);
  std::cout << n[0] << " " << n[1] << " " << n[0] + n[1] << " " << m << " " << I << " " << I1 << " " << (I - I1) / I << std::endl;

}


void TestQuadIntegral(const unsigned & order, const std::vector < std::vector<double> > & weight, const std::vector < std::vector < std::vector<double> > > & x) {

  std::cout << "QUAD ";

  unsigned m = (order % 2 == 0) ? order : order - 1;
  std::vector < unsigned > n(2, m);

  double I = (1. + pow(-1., n[0])) / (n[0] + 1.) * (1. + pow(-1., n[1])) / (n[1] + 1.);

  m = (order > weight.size() - 1) ? weight.size() - 1 : order;

  double I1 = GetIntegral(n, weight[m], x[m]);
  std::cout << n[0] << " " << n[1] << " " << n[0] + n[1] << " " << m << " " << I << " " << I1 << " " << (I - I1) / I << std::endl;

}


void TestHexIntegral(const unsigned & order, const std::vector < std::vector<double> > & weight, const std::vector < std::vector < std::vector<double> > > & x) {

  std::cout << "HEX ";

  unsigned m = (order % 2 == 0) ? order : order - 1;

  std::vector < unsigned > n(3, m);

  double I = (1. + pow(-1., n[0])) / (n[0] + 1.) * (1. + pow(-1., n[1])) / (n[1] + 1) * (1. + pow(-1., n[2])) / (n[2] + 1) ;

  m = (order > weight.size() - 1) ? weight.size() - 1 : order;
  double I1 = GetIntegral(n, weight[m], x[m]);

  std::cout << n[0] << " " << n[1] <<  " " << n[2] << " " << n[0] + n[1] + n[2] << " " << m << " " << I << " " << I1 << " " << (I - I1) / I << std::endl;

}

void TestWedgeIntegral(const unsigned & order, const std::vector < std::vector<double> > & weight, const std::vector < std::vector < std::vector<double> > > & x) {

  std::cout << "WEDGE ";

  std::vector < unsigned > n(3);

  n[0] = order / 2;
  n[1] = order - n[0];
  n[2] = (order % 2 == 0) ? order : order - 1;

  double I = boost::math::tgamma(1 + n[0]) * boost::math::tgamma(1 + n[1]) / boost::math::tgamma(3 + n [0] + n[1]) * (1. + pow(-1., n[2])) / (n[2] + 1) ;

  unsigned m = (order > weight.size() - 1) ? weight.size() - 1 : order;
  double I1 = GetIntegral(n, weight[m], x[m]);

  std::cout << n[0] << " " << n[1] <<  " " << n[2] << " " << n[0] + n[1] << " " << m << " " << I << " " << I1 << " " << (I - I1) / I << std::endl;

}

void TestTetIntegral(const unsigned & order, const std::vector < std::vector<double> > & weight, const std::vector < std::vector < std::vector<double> > > & x) {

  std::cout << "TET ";

  std::vector < unsigned > n(3);
  n[0] = order / 3;
  n[1] = (order - n[0]) / 2;
  n[2] = order - n[0] - n[1];

  double I = boost::math::tgamma(1 + n[0]) * boost::math::tgamma(1 + n[1]) * boost::math::tgamma(1 + n[2]) / boost::math::tgamma(4 + n [0] + n[1] + n[2]);

  unsigned m = (order > weight.size() - 1) ? weight.size() - 1 : order;
  double I1 = GetIntegral(n, weight[m], x[m]);

  std::cout << n[0] << " " << n[1] << " " << n[2] << " " << n[0] + n[1] + n[2] << " " << m << " " << I << " " << I1 << " " << (I - I1) / I << std::endl;

}

void BuildGaussPoints(GeomElType &geomElType, const unsigned &order,
                      std::vector < std::vector < double> > &weight,  std::vector < std::vector < std::vector < double> > > &x) {

  const unsigned GaussSize[6][40] = {{
    }, {
      1, 1, 5, 5, 15, 15, 31, 31, 45, 45
    }, {
      1, 1, 8, 8, 21, 21, 52, 52, 95, 95,168, 168, 259, 259
    }, {

    }, {
      1, 1, 4, 4, 7, 7, 13, 13, 19, 19, 28, 28, 37, 37
    }, {
      1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7  
    }
  };
  const double *Gauss[6][40] = {{
    }, {
      tet_gauss::Gauss0[0], tet_gauss::Gauss0[0],
      tet_gauss::Gauss1[0], tet_gauss::Gauss1[0],
      tet_gauss::Gauss2[0], tet_gauss::Gauss2[0],
      tet_gauss::Gauss3[0], tet_gauss::Gauss3[0],
      tet_gauss::Gauss4[0], tet_gauss::Gauss4[0]
    }, {
      wedge_gauss::Gauss0[0], wedge_gauss::Gauss0[0],
      wedge_gauss::Gauss1[0], wedge_gauss::Gauss1[0],
      wedge_gauss::Gauss2[0], wedge_gauss::Gauss2[0],
      wedge_gauss::Gauss3[0], wedge_gauss::Gauss3[0],
      wedge_gauss::Gauss4[0], wedge_gauss::Gauss4[0],
      wedge_gauss::Gauss5[0], wedge_gauss::Gauss5[0],
      wedge_gauss::Gauss6[0], wedge_gauss::Gauss6[0],
    }, {
    }, {
      tri_gauss::Gauss0[0], tri_gauss::Gauss0[0],
      tri_gauss::Gauss1[0], tri_gauss::Gauss1[0],
      tri_gauss::Gauss2[0], tri_gauss::Gauss2[0],
      tri_gauss::Gauss3[0], tri_gauss::Gauss3[0],
      tri_gauss::Gauss4[0], tri_gauss::Gauss4[0],
      tri_gauss::Gauss5[0], tri_gauss::Gauss5[0],
      tri_gauss::Gauss6[0], tri_gauss::Gauss6[0]
    }, {
      line_gauss::Gauss0[0], line_gauss::Gauss0[0],
      line_gauss::Gauss1[0], line_gauss::Gauss1[0],
      line_gauss::Gauss2[0], line_gauss::Gauss2[0],
      line_gauss::Gauss3[0], line_gauss::Gauss3[0],
      line_gauss::Gauss4[0], line_gauss::Gauss4[0],
      line_gauss::Gauss5[0], line_gauss::Gauss5[0],
      line_gauss::Gauss6[0], line_gauss::Gauss6[0],
      line_gauss::Gauss7[0], line_gauss::Gauss7[0],
      line_gauss::Gauss8[0], line_gauss::Gauss8[0],
      line_gauss::Gauss9[0], line_gauss::Gauss9[0],
      line_gauss::Gauss10[0], line_gauss::Gauss10[0],
      line_gauss::Gauss11[0], line_gauss::Gauss11[0],
      line_gauss::Gauss12[0], line_gauss::Gauss12[0],
      line_gauss::Gauss13[0], line_gauss::Gauss13[0],
      line_gauss::Gauss14[0], line_gauss::Gauss14[0],
      line_gauss::Gauss15[0], line_gauss::Gauss15[0],
      line_gauss::Gauss16[0], line_gauss::Gauss16[0],
      line_gauss::Gauss17[0], line_gauss::Gauss17[0],
      line_gauss::Gauss18[0], line_gauss::Gauss18[0],
      line_gauss::Gauss19[0], line_gauss::Gauss19[0]
    }
  };

  const elem_type *fem;
  std::vector < std::vector < double > > xv;
  unsigned nv;
  unsigned dim;
  std::vector < unsigned > dm;
  unsigned m1 = 0;

  if(HEX == geomElType) {
    nv = 8;
    dim = 3;
    dm.assign(dim, 0);
    fem = new const elem_type_3D("hex", "linear", "seventh");
    xv = {{ -1., 1., 1., -1., -1., 1., 1., -1.},
      { -1., -1., 1., 1., -1., -1., 1., 1.},
      { -1., -1., -1., -1., 1., 1., 1., 1.}
    }; //hex

  }
  else if(TET == geomElType) {

    nv = 8;
    dim = 3;
    dm.assign(dim, 2);
    m1 = 9;

    fem = new const elem_type_3D("hex", "linear", "seventh");
    xv = {{0, 1, 0.5, 0., 0., 0.5, 1. / 3., 0},
      {0, 0, 0.5, 1., 0., 0., 1. / 3., 0.5},
      {0, 0, 0, 0, 1, 0.5, 1. / 3., 0.5}
    }; //tet1
  }
  else if(WEDGE == geomElType) {

    nv = 8;
    dim = 3;
    dm.assign(dim, 1);
    dm[2] = 0;
    m1 = 14;

    fem = new const elem_type_3D("hex", "linear", "seventh");
    xv = {{0., 1., 0.5, 0., 0., 1., 0.5, 0.},
      {0., 0., 0.5, 1., 0., 0., 0.5, 1.},
      { -1., -1., -1., -1., 1., 1., 1., 1.}
    }; //wedge1
  }
  else if(QUAD == geomElType) {
    nv = 4;
    dim = 2;
    dm.assign(dim, 0);

    fem = new const elem_type_2D("quad", "linear", "seventh");
    xv = {{ -1., 1., 1., -1.}, { -1., -1., 1., 1.}}; //quad
  }
  else if(TRI == geomElType) {
    nv = 4;
    dim = 2;
    dm.assign(dim, 1);
    m1 = 14;

    fem = new const elem_type_2D("quad", "linear", "seventh");
    xv = {{0., 1., 0.5, 0.}, {0., 0., 0.5, 1.}}; // tri

  }
  else {
    abort();
  }


  std::vector <double> phi;  // local test function
  std::vector <double> phi_x; // local test function first order partial derivatives
  double jac;

  weight.resize(order + 1);
  x.resize(order + 1);
  for(unsigned m = 0; m < order + 1; m++) {
    x[m].resize(dim);
  }

  for(unsigned m = 0; m < m1; m++) {
    //if(geomElType != WEDGE) {
      unsigned size = GaussSize[geomElType][m];
      weight[m].resize(size);
      for(unsigned k = 0; k < dim; k++) {
        x[m][k].resize(size);
      }

      for(unsigned cnt = 0; cnt < size ; cnt++) {
        weight[m][cnt] = Gauss[geomElType][m][cnt];
        for(unsigned k = 0; k < dim; k++) {
          x[m][k][cnt] = Gauss[geomElType][m][size * (k + 1) + cnt];
        }
      }
//     //}
//     else{
//       unsigned sizeT = GaussSize[TRI][m];
//       unsigned sizeL = GaussSize[LINE][m];
//       unsigned size = sizeT * sizeL;
//       std::cout << "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA" <<sizeT << " " << sizeL <<" "<< size << std::endl;
//       weight[m].resize(size);
//       for(unsigned k = 0; k < dim; k++) {
//         x[m][k].resize(size);
//       }

//       unsigned cnt = 0;
//       for(unsigned i = 0; i < sizeT ; i++) {
//         for(unsigned j = 0; j < sizeL ; j++) {  
//           weight[m][cnt] = Gauss[TRI][m][i] * Gauss[LINE][m][j];
//           x[m][0][cnt] = Gauss[TRI][m][sizeT * 1 + i];
//           x[m][1][cnt] = Gauss[TRI][m][sizeT * 2 + i];
//           x[m][2][cnt] = Gauss[LINE][m][sizeL * 1 + j];
//           cnt++;
//         }
//       }  
//     }

  }

  for(unsigned m = m1; m < order + 1; m++) {
    unsigned size = 1;
    std::vector < unsigned > ng(dim);
    for(unsigned k = 0; k < dim; k++) {
      ng[k] = (m + dm[k]) / 2 + 1;
      size *= ng[k];
    }

    weight[m].resize(size);
    for(unsigned k = 0; k < dim; k++) {
      x[m][k].resize(size);
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
      weight[m][cnt] = 1.;
      for(unsigned k = 0; k < dim; k++) {
        weight[m][cnt] *= Gauss[LINE][m + dm[k] ][I[k]];
        xi[k] = Gauss[LINE][m + dm[k] ][ng[k] + I[k]];
      }
      fem->Jacobian(xv, xi, jac, phi, phi_x);
      weight[m][cnt] *= jac;
      for(unsigned k = 0; k < dim; k++) {
        x[m][k][cnt] = 0.;
        for(unsigned ii = 0; ii < nv; ii++) {
          x[m][k][cnt] += phi[ii] * xv[k][ii];
        }
      }
    }
  }

  delete fem;
}

