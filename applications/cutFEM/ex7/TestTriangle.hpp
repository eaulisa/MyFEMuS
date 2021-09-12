
#ifndef __femus_cut_fem_TRItest_hpp__
#define __femus_cut_fem_TRItest_hpp__

#include "Triangle.hpp"

template <class Float1>
void TestTriangle(const Float1 &eps) {

  typedef typename boost::math::tools::promote_args<Float1>::type myType;
  typedef boost::multiprecision::cpp_bin_float_oct myTypeB;

  myType eps1 = 1.0e-12;

  std::cout << "testing the Triangle \n";
  std::vector<unsigned>m = {9, 9};
  //std::vector<myType>a = { -1e-05, -1e-5};
  std::vector<myType>a = { 0., 0. };
//   std::vector<myTypeB> af(2);
//   af[0] = a[0];
//   af[1] = a[1];

  myType d = 1;
//   myTypeB df = d;

  int s = 4;

  std::cout.precision(16);

  std::vector<std::vector<myType>> epsCut{{0, 0, 0},
    {eps1, 0, 0}, {-eps1, 0, 0}, {0, eps1, 0}, {0, -eps1, 0}, {0, 0, eps1}, {0, 0, -eps1},
    {eps1, eps1, 0}, {eps1, -eps1, 0}, {-eps1, eps1, 0}, {-eps1, -eps1, 0},
    {eps1, -eps1, eps1},  {-eps1, -eps1, eps1}, {eps1, -eps1, -eps1}, {-eps1, -eps1, -eps1}
  };
//   std::vector<std::vector<myType>> smallCut{{0, 0, 0}, {0, 0, 1}, {0, 0, -1},
//     {-1, -1, 0}, {-1, 1, 0}, {1, 1, -2}, {1, 0, -1}, {1, -1, -2}, {0, -1, 0},
//     {1, 1, 0}, {1, -1, 0}, {-1, -1, 2}, {-1, 0, 1}, {-1, 1, 2}, {0, 1, 0}
//   };

  std::vector<std::vector<myType>> smallCut{{0, 0, 0}, {0, 0, 1}, {0, 0, -1},
    {-1, -1, 0}, {-1, 0, 0}, {0, 1, -1}, {1, 1, -1}, {1, 0, -1}, {0, -1, 0},
    {1, 1, 0}, {1, 0, 0}, {0, -1, 1}, {-1, -1, 1}, {-1, 0, 1}, {0, 1, 0}
  };


  for(unsigned i = 0; i < smallCut.size(); i++) {
    for(unsigned j = 0; j < epsCut.size(); j++) {

      a[0] = smallCut[i][0] + epsCut[j][0];
      a[1] = smallCut[i][1] + epsCut[j][1];
      d = smallCut[i][2] + epsCut[j][2];

//       af[0] = static_cast<myType>(a[0]);
//       af[1] = static_cast<myType>(a[1]);
//       df = static_cast<myType>(d);

//       myType I1 = TriangleA(s,  m, a, d);
//       myTypeB I2 = TriangleA(s,  m, af, df);

      myType I1 = Triangle<myType, myType>(s, m, a, d);
      myType I2 = Triangle<myType, myTypeB>(s, m, a, d);

      if((I2 != 0. &&  fabs((I1 - I2) / I2) < eps) || I1 == 0.) {
        //std::cout << "passed, i = " << i << " j = " << j << " ";
        //std::cout << "s = " << s << " a = " << a[0] << " b = " << a[1] << " d = " << d << " I = " << I1 << std::endl;
      }
      else {
        std::cout << "Warning failed, i = " << i << " j = " << j << " ";
        std::cout << "s = " << s << " a = " << a[0] << " b = " << a[1] << " d = " << d << " I1 =" << I1 << " I2 = " << I2 << std::endl;
      }
    }
    //std::cout<<std::endl;
  }


  std::vector<myTypeB> af(2);
  myTypeB df = d;
  

  for(unsigned i = 0; i < 1000; i++) {
    a[0] = 2. * (0.5 * RAND_MAX - static_cast <myType>(rand())) / RAND_MAX;
    a[1] = sqrt(1 - a[0] * a[0]);
    d = 10. * (0.5 * RAND_MAX - static_cast <myType>(rand())) / RAND_MAX;

    if(a[0] + a[1] + d <= 0) {
      a[0] = -a[0];
      a[1] = -a[1];
      d = -d;
    }

    af[0] = a[0];
    af[1] = a[1];
    df = d;

    if(fabs(TriangleA(-1,  m, a, d) - TriangleBFull(-1,  m, af, df)) > eps) {
      std::cout << "surface test failed" << std::endl;
      std::cout << a[0] << " " << a[1] << " " << " " << d << "\n";
      std::cout << TriangleBFull(-1,  m, a, d) << " " << TriangleBFull(-1,  m, af, df) << std::endl;
    }

    if(fabs(TriangleA(0,  m, a, d) - TriangleBFull(0,  m, af, df)) > eps) {
      std::cout << "volume test failed" << std::endl;
      std::cout << a[0] << " " << a[1] << " " << d << "\n";
      std::cout << TriangleBFull(0,  m, a, d) << " " << TriangleBFull(0,  m, af, df) << std::endl;
    }
  }

  for(unsigned i = 0; i < 1000; i++) {

    a[0] = -0.05;
    a[1] = -a[0];
    d = 500. * (0.5 * RAND_MAX - static_cast <myType>(rand())) / RAND_MAX;

    if(a[0] + a[1] + d <= 0) {
      a[0] = -a[0];
      a[1] = -a[1];
      d = -d;
    }

    af[0] = a[0];
    af[1] = a[1];
    df = d;

    if(fabs(TriangleA(-1,  m, a, d) - TriangleBFull(-1,  m, af, df)) > eps) {
      std::cout << "surface test failed" << std::endl;
      std::cout << a[0] << " " << a[1] << " " << " " << d << "\n";
      std::cout << TriangleA(-1,  m, a, d) << " " << TriangleBFull(-1,  m, af, df) << std::endl;
    }

    if(fabs(TriangleA(0,  m, a, d) - TriangleBFull(0,  m, af, df)) > eps) {
      std::cout << "volume test failed" << std::endl;
      std::cout << a[0] << " " << a[1] << " " << d << "\n";
      std::cout << TriangleA(0,  m, a, d) << " " << TriangleBFull(0,  m, af, df) << std::endl;
    }
  }

  for(unsigned i = 0; i < 1000; i++) {

    a[0] = 0.;
    a[1] = 1.;
    d = 5. * (0.5 * RAND_MAX - static_cast <myType>(rand())) / RAND_MAX;

    if(a[0] + a[1] + d <= 0) {
      a[0] = -a[0];
      a[1] = -a[1];
      d = -d;
    }

    af[0] = a[0];
    af[1] = a[1];
    df = d;

    if(fabs(TriangleA(-1,  m, a, d) - TriangleBFull(-1,  m, af, df)) > eps) {
      std::cout << "surface test failed" << std::endl;
      std::cout << a[0] << " " << a[1] << " " << " " << d << "\n";
      std::cout << TriangleA(-1,  m, a, d) << " " << TriangleBFull(-1,  m, af, df) << std::endl;
    }

    if(fabs(TriangleA(0,  m, a, d) - TriangleBFull(0,  m, af, df)) > eps) {
      std::cout << "volume test failed" << std::endl;
      std::cout << a[0] << " " << a[1] << " " << d << "\n";
      std::cout << TriangleA(0,  m, a, d) << " " << TriangleBFull(0,  m, af, df) << std::endl;
    }
  }

  for(unsigned i = 0; i < 1000; i++) {

    a[0] = 1.;
    a[1] = 0.;
    d = 5. * (0.5 * RAND_MAX - static_cast <myType>(rand())) / RAND_MAX;

    if(a[0] + a[1] + d <= 0) {
      a[0] = -a[0];
      a[1] = -a[1];
      d = -d;
    }

    af[0] = a[0];
    af[1] = a[1];
    df = d;

    if(fabs(TriangleA(-1,  m, a, d) - TriangleBFull(-1,  m, af, df)) > eps) {
      std::cout << "surface test failed" << std::endl;
      std::cout << a[0] << " " << a[1] << " " << " " << d << "\n";
      std::cout << TriangleA(-1,  m, a, d) << " " << TriangleBFull(-1,  m, af, df) << std::endl;
    }

    if(fabs(TriangleA(0,  m, a, d) - TriangleBFull(0,  m, af, df)) > eps) {
      std::cout << "volume test failed" << std::endl;
      std::cout << a[0] << " " << a[1] << " " << d << "\n";
      std::cout << TriangleA(0,  m, a, d) << " " << TriangleBFull(0,  m, af, df) << std::endl;
    }
  }
}
#endif
