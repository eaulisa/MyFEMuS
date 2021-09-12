
#ifndef __femus_cut_fem_HCItest_hpp__
#define __femus_cut_fem_HCItest_hpp__

#include "HyperCube.hpp"

template <class myType>
void TestQuad(myType &eps) {

  std::cout << "testing the Quadrilateral \n";

  //typedef double myTypeB;
  //typedef cpp_bin_float_quad myTypeB;
  typedef cpp_bin_float_oct myTypeB;

  std::vector<myType> a = {-1.1, -1};
  myType norm = sqrt(a[0] * a[0] + a[1] * a[1]);
  a[0] /= norm;
  a[1] /= norm;

  std::vector<unsigned> m = {0, 0};
  myType d = 0. / norm;

  std::cout.precision(16);

  std::cout << HyperCube<myType, myType>(0, m, a, d) << std::endl;


  for(unsigned i = 0; i < a.size(); i++) {
    d -= a[i];
    a[i] *= 2;
  }

  std::vector<myType> ma(2);
  ma[0] = -a[0];
  ma[1] = -a[1];

  std::cout << 4 * HyperCubeA1(1, 0, m, a, ma, d, -d) << std::endl;
  std::cout << 4 * HyperCubeC(1, 0, m, a, ma, d, -d) << std::endl << std::endl;

  ///////////////////////////////////////////////////////////////////

  m = {13, 5};

  std::cout.precision(14);

  unsigned s = 5;

  myType eps1 = 1.0e-12;
  std::vector<std::vector<myType>> epsCut{{0, 0, 0},
    {eps1, 0, 0}, {-eps1, 0, 0}, {0, eps1, 0}, {0, -eps1, 0}, {0, 0, eps1}, {0, 0, -eps1},
    {eps1, eps1, 0}, {eps1, -eps1, 0}, {-eps1, eps1, 0}, {-eps1, -eps1, 0},
    {eps1, -eps1, eps1}, {-eps1, -eps1, eps1}, {eps1, -eps1, -eps1}, {-eps1, -eps1, -eps1}
  };

  std::vector<std::vector<myType>> smallCut{{0, 0, 0}, {0, 0, 1}, {0, 0, -1},
    {-1, -1, -2}, {-1, 0, -1}, {-1, 1, -2}, {0, 1, -1},
    {1, 1, -2}, {1, 0, -1}, {1, -1, -2}, {0, -1, -1},
    {1, 1, 2}, {1, 0, 1}, {1, -1, 2}, {0, -1, 1},
    {-1, -1, 2}, {-1, 0, 1}, {-1, 1, 2}, {0, 1, 1}
  };

  for(unsigned i = 0; i < smallCut.size(); i++) {
    for(unsigned j = 0; j < epsCut.size(); j++) {

      a.resize(2);
      a[0] = myType(smallCut[i][0] + epsCut[j][0]);
      a[1] = myType(smallCut[i][1] + epsCut[j][1]);
      d = myType(smallCut[i][2] + epsCut[j][2]);

      myType I1 = HyperCube<myType, myType>(s,  m, a, d);
      myType I2 = HyperCube<myType, myTypeB>(s,  m, a, d);

      if((I2 != 0. &&  fabs((I1 - I2) / I2) < eps) || I1 == 0.) {
        //std::cout << "passed "<< i << " " << j<<" ";
        //std::cout << "s = " << s << " a = " << a[0] << " b = " << a[1] << " d = " << d << " I = " << I1 << std::endl;
      }
      else {
        std::cout << "Warning failed " << i << " " << j << " ";
        std::cout << "s = " << s << " a = " << a[0] << " b = " << a[1] << " d = " << d << " I1 =" << I1 << " I2 = " << I2 << std::endl;
      }
    }
  }
}

template <class myType>
void TestHex(myType &eps) {

  std::cout << "testing the Hexahedron \n";

  //typedef long double myTypeB;
  //typedef cpp_bin_float_quad myTypeB;
  typedef cpp_bin_float_oct myTypeB;

  std::vector<myType> a = {-1.1, 1, -0.98};
  myType norm = sqrt(a[0] * a[0] + a[1] * a[1] + a[2] * a[2]);

  a[0] /= norm;
  a[1] /= norm;
  a[2] /= norm;

  std::vector<unsigned> m = {0, 0, 0};
  myType d = 0 / norm;

  std::cout << HyperCube<myType, myType>(0, m, a, d) << std::endl;

  for(unsigned i = 0; i < a.size(); i++) {
    d -= a[i];
    a[i] *= 2;
  }

  std::vector<myType> ma(3);
  ma[0] = -a[0];
  ma[1] = -a[1];
  ma[2] = -a[2];

  std::cout << 8 * HyperCubeB(2, 0, m, a, ma, d, -d) << std::endl;
  std::cout << 8 * HyperCubeC(2, 0, m, a, ma, d, -d) << std::endl;

  ////////////////////////////////////////

  m = {5, 6, 7};

  std::cout.precision(14);

  double eps1 = 1.0e-12;
  std::vector<std::vector<double>> epsCut = { {0, 0, 0, 0},
    {eps1, 0, 0, 0}, {-eps1, 0, 0, 0},
    {0, eps1, 0, 0}, {0, -eps1, 0, 0},
    {0, 0, eps1, 0}, {0, 0, -eps1, 0},
    {0, 0, 0, eps1}, {0, 0, 0, -eps1},
    {eps1, eps1, eps1, 0}, {eps1, -eps1, eps1, 0},
    {eps1, eps1, -eps1, 0}, {eps1, -eps1, -eps1, 0},
    {-eps1, eps1, eps1, 0}, {-eps1, -eps1, eps1, 0},
    {-eps1, eps1, -eps1, 0}, {-eps1, -eps1, -eps1, eps1},
    {eps1, eps1, eps1, eps1}, {eps1, -eps1, eps1, eps1},
    {eps1, eps1, -eps1, eps1}, {eps1, -eps1, -eps1, eps1},
    {-eps1, eps1, eps1, eps1}, {-eps1, -eps1, eps1, eps1},
    {-eps1, eps1, -eps1, eps1}, {-eps1, -eps1, -eps1, eps1},
    {eps1, eps1, eps1, -eps1}, {eps1, -eps1, eps1, -eps1},
    {eps1, eps1, -eps1, -eps1}, {eps1, -eps1, -eps1, -eps1},
    {-eps1, eps1, eps1, -eps1}, {-eps1, -eps1, eps1, -eps1},
    {-eps1, eps1, -eps1, -eps1}, {-eps1, -eps1, -eps1, -eps1},
  };

  std::vector<std::vector<double>> smallCut = {
    {0, 0, 0, 0}, {0, 0, 0, 1}, {0, 0, 0, -1},
    {-1, -1, -1, -3}, {-1, 1, -1, -3}, {1, 1, -1, -3}, {1, -1, -1, -3},
    {-1, -1, 1, -3}, {-1, 1, 1, -3}, {1, 1, 1, -3}, {1, -1, 1, -3},
    {1, 1, 1, 3}, {1, -1, 1, 3}, {-1, -1, 1, 3}, {-1, 1, 1, 3},
    {1, 1, -1, 3}, {1, -1, -1, 3}, {-1, -1, -1, 3}, {-1, 1, -1, 3},
    {-1, 0, 0, -1}, {1, 0, 0, - 1}, {0, -1, 0, -1}, {0, 1, 0, -1}, {0, 0, -1, -1}, {0, 0, 1, -1},
    {1, 0, 0, 1}, {-1, 0, 0, 1}, {0, 1, 0, 1}, {0, -1, 0, 1}, {0, 0, 1, 1}, {0, 0, -1, 1},
  };


  int s = 0;

  for(unsigned i = 0; i < smallCut.size(); i++) {
    for(unsigned j = 0; j < epsCut.size(); j++) {
      a.resize(3);
      a[0] = myType(smallCut[i][0] + epsCut[j][0]);
      a[1] = myType(smallCut[i][1] + epsCut[j][1]);
      a[2] = myType(smallCut[i][2] + epsCut[j][2]);
      d = myType(smallCut[i][3] + epsCut[j][3]);

      myType I1 = HyperCube<myType, myType>(s,  m, a, d);
      myType I2 = HyperCube<myType, myTypeB>(s,  m, a, d);

      if((I2 != 0. &&  fabs((I1 - I2) / I2) < eps) || I1 == 0.) {
        //std::cout << "passed " << i << " " << j << " ";
        //std::cout << "s = " << s << " a = " << a[0] << " b = " << a[1] << " c = " << a[2] << " d = " << d << " I = " << I1 << std::endl;
      }
      else {
        std::cout << "Warning failed " << i << " " << j << " ";
        std::cout << "s = " << s << " a = " << a[0] << " b = " << a[1] << " c = " << a[2] << " d = " << d << " I1 =" << I1 << " I2 = " << I2 << std::endl;
      }

    }
  }

}

#endif
