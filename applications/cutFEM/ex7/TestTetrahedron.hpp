
#ifndef __femus_cut_fem_TETtest_hpp__
#define __femus_cut_fem_TETtest_hpp__

#include "Tetrahedron.hpp"

template <class myType>
void TestTetrahedron(myType &eps) {



  std::cout << "testing the Tetrahedron \n";
  //typedef long double myTypeB;
  //typedef cpp_bin_float_quad myTypeB;
  typedef cpp_bin_float_oct myTypeB;

  std::vector<myType> a = {0.98, 1, 1.1};
  myType norm = sqrt(a[0] * a[0] + a[1] * a[1] + a[2] * a[2]);

  a[0] /= norm;
  a[1] /= norm;
  a[2] /= norm;

  std::vector<unsigned> m = {3, 4, 5};
  myType d = -0.25 / norm;

  std::cout.precision(14);

  int s = 0;
  
  std::cout << TetrahedronB(s, m, a, d) << std::endl;
  std::cout << TetrahedronC(s, m, a, d) << std::endl;

  std::cout << Tetrahedron<myType, myType>(s, m, a, d) << std::endl;
  std::cout << Tetrahedron<myType, myTypeB>(s, m, a, d) << std::endl;

// return;

  double eps1 = 1.e-12;
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
    {-1, -1, -1, 0}, {1, 0, 0, -1}, {0, 1, 0, -1}, {0, 0, 1, -1},
    {0, 0, -1, 0}, {0, -1, 0, 0}, {-1, 0, 0, 0}, {1, 1, 1, -1},
    {1, 1, 1, 0}, {-1, 0, 0, 1}, {0, -1, 0, 1}, {0, 0, -1, 1},
    {0, 0, 1, 0}, {0, 1, 0, 0}, {1, 0, 0, 0}, {-1, -1, -1, 1}
  };

  s = 0;

  for(unsigned i = 0; i < smallCut.size(); i++) {
    for(unsigned j = 0; j < epsCut.size(); j++) {
      a.resize(3);


      a[0] = myType(smallCut[i][0] + epsCut[j][0]);
      a[1] = myType(smallCut[i][1] + epsCut[j][1]);
      a[2] = myType(smallCut[i][2] + epsCut[j][2]);
      d = myType(smallCut[i][3] + epsCut[j][3]);

      myType I1 = Tetrahedron<myType, myType>(s, m, a, d);
      myType I2 = Tetrahedron<myType, myTypeB>(s, m, a, d);

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
