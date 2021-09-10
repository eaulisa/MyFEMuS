
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

  std::vector<unsigned> m = {1,3, 2};
  myType d = 0 / norm;

  std::cout.precision(14);

  //std::cout << Tetrahedron(1, m, a, d) << std::endl;
  //std::cout << TetrahedronA(1, m, a, d) << std::endl;
  std::cout << TetrahedronB(0, m, a, d) << std::endl;
  std::cout << TetrahedronC(0, m, a, d) << std::endl;

  return;
  
  std::vector<myTypeB> af(3);
  af[0] = a[0];
  af[1] = a[1];
  af[2] = a[2];

  myTypeB df = d;

//   std::cout << TetrahedronA(1, m, af, df) << std::endl;
//   std::cout << TetrahedronB(1, m, af, df) << std::endl;
  //std::cout << TetrahedronB1(1, m, af, df) << std::endl;

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

  int s = 0;

  std::vector<myType> a1(3);
  myType d1;

//   for(unsigned i = 8; i <9 + 0*smallCut.size(); i++) {
//     for(unsigned j = 7; j <8 + 0* epsCut.size(); j++) {

  for(unsigned i = 0; i < smallCut.size(); i++) {
    for(unsigned j = 0; j < epsCut.size(); j++) {
      a.resize(3);
      a1.resize(3);
      af.resize(3);

      a1[0] = myType(smallCut[i][0] + epsCut[j][0]);
      a1[1] = myType(smallCut[i][1] + epsCut[j][1]);
      a1[2] = myType(smallCut[i][2] + epsCut[j][2]);
      d1 = myType(smallCut[i][3] + epsCut[j][3]);

      if(fabs(a1[0]) <= fabs(a1[1]) && fabs(a1[0]) <= fabs(a1[2])) {
        if(fabs(a1[1]) <= fabs(a1[2])) {
          a = {-a1[0], a1[1], a1[2]};
          d = d1 + a1[0];
        }
        else {
          a = {-a1[0], a1[2], a1[1]};
          d = d1 + a1[0];
        }
      }
      else if(fabs(a1[1]) <= fabs(a1[2])) {
        if(fabs(a1[0]) <= fabs(a1[2])) {
          a = {-a1[1], a1[0], a1[2]};
          d = d1 + a1[1];
        }
        else {
          a = {-a1[1], a1[2], a1[0]};
          d = d1 + a1[1];
        }
      }
      else {
        if(fabs(a1[0]) <= fabs(a1[1])) {
          a = {-a1[2], a1[0], a1[1]};
          d = d1 + a1[2] ;
        }
        else {
          a = {-a1[2], a1[1], a1[0]};
          d = d1 + a1[2];
        }
      }

//       a=a1;
//       d=d1;

      af[0] = myTypeB(a[0]);
      af[1] = myTypeB(a[1]);
      af[2] = myTypeB(a[2]);
      df = myTypeB(d);

      myType I1 = TetrahedronCast(s,  m, a, d);
      myTypeB I2 = TetrahedronCast(s,  m, af, df);

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

  std::cout <<" cast = "<< cast <<" uncast = "<< uncast <<std::endl;


}
#endif
