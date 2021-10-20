
#include "Line.hpp"
#include "Square.hpp"
#include "Cube.hpp"
#include "HyperCube.hpp"
#include "Triangle.hpp"
#include "Tetrahedron.hpp"


#include <iostream>
#include <boost/multiprecision/cpp_bin_float.hpp>


namespace boost {
  namespace multiprecision {
    typedef number < backends::cpp_bin_float < 24, backends::digit_base_2, void, boost::int16_t, -126, 127 >, et_off >         cpp_bin_float_single;
    typedef number < backends::cpp_bin_float < 53, backends::digit_base_2, void, boost::int16_t, -1022, 1023 >, et_off >       cpp_bin_float_double;
    typedef number < backends::cpp_bin_float < 64, backends::digit_base_2, void, boost::int16_t, -16382, 16383 >, et_off >     cpp_bin_float_double_extended;
    typedef number < backends::cpp_bin_float < 113, backends::digit_base_2, void, boost::int16_t, -16382, 16383 >, et_off >    cpp_bin_float_quad;
    typedef number < backends::cpp_bin_float < 237, backends::digit_base_2, void, boost::int32_t, -262142, 262143 >, et_off >  cpp_bin_float_oct;
  }
} // namespaces

using boost::multiprecision::cpp_bin_float_oct;
using boost::multiprecision::cpp_bin_float_quad;


//#include "./old/TestHyperCubeOld.hpp"
//#include "./old/TetrahedronOld.hpp"
//#include "./old/HyperCubeOld.hpp"

#include "TestHyperCube.hpp"
// #include "TestTriangle.hpp"
// #include "TestPrism.hpp"



int main(int, char**) {

  typedef double Type;

  bool quad = false;
  bool hex = true;
//     bool triangle = false;
  bool tetrahedron = false; //true;
//     bool prism = false;
//
  Type eps = 1.0e-11;

  int s = 0;
  unsigned qMax = 8;
  double a = 0.8;
  double b = 0.9;
  double c = 1;
  double d = 0;

  if(quad) {
    clock_t time = clock();
    //TestQuad(eps);
    //TestQuadOld(eps);
    std::cout << "Time = " << static_cast<double>((clock() - time)) / CLOCKS_PER_SEC << std::endl;

    SQImap <double, double> sqr(qMax, 0);

    for(unsigned tt = 0; tt < 10000; tt++) {
      sqr.clear();
      a = -1. + 2. * rand() / RAND_MAX;
      b = -1. + 2. * rand() / RAND_MAX;
      d = -1. + 2. * rand() / RAND_MAX;
      for(unsigned q = 0; q <= qMax; q++) {
        for(unsigned j = 0; j <= q; j++) {
          unsigned i = q - j;
          sqr(-1, {i, j}, {a, b}, d);
          sqr(0, {i, j}, {a, b}, d);
        }
      }
    }
    sqr.printCounter();
    std::cout << "Map Time = " << static_cast<double>((clock() - time)) / CLOCKS_PER_SEC << std::endl;

    HCImap <double, double> sqr2(2, qMax, 0);
    time = clock();
    for(unsigned tt = 0; tt < 10000; tt++) {
      sqr2.clear();
      a = -1. + 2. * rand() / RAND_MAX;
      b = -1. + 2. * rand() / RAND_MAX;
      d = -1. + 2. * rand() / RAND_MAX;
      for(unsigned q = 0; q <= qMax; q++) {
        for(unsigned j = 0; j <= q; j++) {
          unsigned i = q - j;
          sqr2(-1, {i, j}, {a, b}, d);
          sqr2(0, {i, j}, {a, b}, d);
        }
      }
    }
    sqr2.printCounter();
    std::cout << "Old Time = " << static_cast<double>((clock() - time)) / CLOCKS_PER_SEC << std::endl;
  }



  if(hex) {

    TestHex(eps);
    //TestHexOld(eps);


    CBImap <double, double> cube(qMax, 0);
    srand(0);
    clock_t time = clock();
    for(unsigned tt = 0; tt < 1000; tt++) {
      cube.clear();
      a = -1. + 2. * rand() / RAND_MAX;
      b = -1. + 2. * rand() / RAND_MAX;
      c = -1. + 2. * rand() / RAND_MAX;
      d = -1. + 2. * rand() / RAND_MAX;
      for(unsigned q = 0; q <= qMax; q++) {
        for(int ii = q; ii >= 0; ii--) {
          for(int jj = q - ii; jj >= 0; jj--) {
            unsigned i = static_cast<unsigned>(ii);
            unsigned j = static_cast<unsigned>(jj);
            unsigned k = q - i - j;
            cube(-1, {i, j, k}, {a, b, c}, d);
            cube(0, {i, j, k}, {a, b, c}, d);
          }
        }
      }
      //cube.printCounter();
    }
    std::cout << "Map Time = " << static_cast<double>((clock() - time)) / CLOCKS_PER_SEC << std::endl;
    cube.printCounter();


    HCImap <double, double> hci3(3, qMax, 0);
    srand(0);
    time = clock();
    for(unsigned tt = 0; tt < 1000; tt++) {
      hci3.clear();
      a = -1. + 2. * rand() / RAND_MAX;
      b = -1. + 2. * rand() / RAND_MAX;
      c = -1. + 2. * rand() / RAND_MAX;
      d = -1. + 2. * rand() / RAND_MAX;
      for(unsigned q = 0; q <= qMax; q++) {
        for(int ii = q; ii >= 0; ii--) {
          for(int jj = q - ii; jj >= 0; jj--) {
            unsigned i = static_cast<unsigned>(ii);
            unsigned j = static_cast<unsigned>(jj);
            unsigned k = q - i - j;
            hci3(-1, {i, j, k}, {a, b, c}, d);
            hci3(0, {i, j, k}, {a, b, c}, d);
          }
        }
      }
    }
    std::cout << "Old Time = " << static_cast<double>((clock() - time)) / CLOCKS_PER_SEC << std::endl;
    hci3.printCounter();


  }

  return 1;

//
//
//     eps = 5.0e-11;
//     if(triangle) TestTriangle(eps);
//
  if(tetrahedron) {
    //TestTetrahedron(eps);
    clock_t time = clock();
    TTImap <double, double> tet(3, qMax, 0);

    for(unsigned tt = 0; tt < 1000; tt++) {
      tet.clear();
      for(unsigned q = 0; q <= qMax; q++) {
        for(int ii = q; ii >= 0; ii--) {
          for(int jj = q - ii; jj >= 0; jj--) {
            unsigned i = static_cast<unsigned>(ii);
            unsigned j = static_cast<unsigned>(jj);
            unsigned k = q - i - j;
            tet(-1, {i, j, k}, {a, b, c}, d);
            tet(0, {i, j, k}, {a, b, c}, d);
          }
        }
      }
    }
    std::cout << "Map Time = " << static_cast<double>((clock() - time)) / CLOCKS_PER_SEC << std::endl;

    time = clock();
    for(unsigned tt = 0; tt < 1000; tt++) {
      for(unsigned q = 0; q <= qMax; q++) {
        for(int ii = q; ii >= 0; ii--) {
          for(int jj = q - ii; jj >= 0; jj--) {
            unsigned i = static_cast<unsigned>(ii);
            unsigned j = static_cast<unsigned>(jj);
            unsigned k = q - i - j;
            tet(-1, {i, j, k}, {a, b, c}, d);
            tet(0, {i, j, k}, {a, b, c}, d);
          }
        }
      }
    }
    std::cout << "Old Time = " << static_cast<double>((clock() - time)) / CLOCKS_PER_SEC << std::endl;

  }



//
//     if(prism) TestPrism(eps);




//   LSImap <double> lsi(qMax);
//   for(unsigned i = 0; i <= qMax ; i++) {
//     std::cout << lsi(-1, i, a, d) << " ";
//     std::cout << lsi(s, i, a, d) << std::endl;
//   }
//   for(unsigned i = 0; i <= qMax ; i++) {
//     std::cout << lsi(-1, i, a, d) << " ";
//     std::cout << lsi(s, i, a, d) << std::endl;
//   }
//   lsi.clear();
//   for(unsigned i = 0; i <= qMax ; i++) {
//     std::cout << lsi(-1, i, a, d) << " ";
//     std::cout << lsi(s, i, a, d) << std::endl;
//   }
//
//   std::cout << "Triangle\n";
//
//   TRImap <double, double> tri(qMax); // the first template argument is for the IO the second is for the algebra
//   tri.clear();
//   for(unsigned q = 0; q <= qMax; q++) {
//     for(unsigned j = 0; j <= q; j++) {
//       unsigned i = q - j;
//       std::cout << q << " " << i << " " << j << " ";
//       std::cout << tri(-1, {i, j}, {a, b}, d) << " ";
//       std::cout << tri(-1, {i, j}, {-a, -b}, -d) << " ";
//       std::cout << tri(s, {i, j}, {a, b}, d) << " ";
//       std::cout << tri(s, {i, j}, {-a, -b}, -d) << " ";
//       std::cout << tri(s, {i, j}, {a, b}, d) + tri(s, {i, j}, {-a, -b}, -d) - 1. / ((j + 1) * (i + j + 2)) << std::endl;
//     }
//   }
//
//
//
//   TRImap <double, cpp_bin_float_oct> tri2(qMax, 3);
//   tri2.clear();
//   for(unsigned q = 0; q <= qMax; q++) {
//     for(unsigned j = 0; j <= q; j++) {
//       unsigned i = q - j;
//       std::cout << q << " " << i << " " << j << " ";
//       std::cout << tri2(-1, {i, j}, {a, b}, d) << " ";
//       std::cout << tri2(-1, {i, j}, {-a, -b}, -d) << " ";
//       std::cout << tri2(3, {i, j}, {a, b}, d) << " ";
//       std::cout << tri2(3, {i, j}, {-a, -b}, -d) << " ";
//       std::cout << tri2(s, {i, j}, {a, b}, d) + tri2(s, {i, j}, {-a, -b}, -d) - 1. / ((j + 1) * (i + j + 2)) << std::endl;
//     }
//   }

//   std::cout << "Tetrahedron\n";
//
//   s = 0;
//   TTImap <double, double> tti(qMax, s);
//   tti.clear();
//
//   unsigned cnt = 0;
//   for(unsigned q = 0; q <= qMax; q++) {
//     for(int ii = q; ii >= 0; ii--) {
//       for(int jj = q - ii; jj >= 0; jj--) {
//         unsigned i = static_cast<unsigned>(ii);
//         unsigned j = static_cast<unsigned>(jj);
//         unsigned k = q - i - j;
//         cnt++;
//         std::cout << cnt << " " << q << " " << i << " " << j << " " << k << " ";
//         std::cout << tti(-1, {i, j, k}, {a, b, c}, d) - Tetrahedron<double, double>(-1, {i, j, k}, {a, b, c}, d) << " ";
//         std::cout << tti(-1, {i, j, k}, {-a, -b, -c}, -d) - Tetrahedron<double, double>(-1, {i, j, k}, {-a, -b, -c}, -d) << " ";
//         std::cout << tti(s, {i, j, k}, {-a, -b, -c}, -d) - Tetrahedron<double, double>(s, {i, j, k}, {-a, -b, -c}, -d) << " ";
//         std::cout << tti(s, {i, j, k}, {a, b, c}, d) - Tetrahedron<double, double>(s, {i, j, k}, {a, b, c}, d) << " ";
//         std::cout << tti(0, {i, j, k}, {c, b, a}, d) - tti(0, {i, j, k}, {-c, -b, -a}, -d) - 1. / ((1 + k) * (2 + j + k) * (3 + i + j + k)) << " ";
//         std::cout << std::endl;
//       }
//     }
//   }
//   tti.printCounter();
//
//   for(unsigned q = 0; q <= qMax; q++) {
//     for(int ii = q; ii >= 0; ii--) {
//       for(int jj = q - ii; jj >= 0; jj--) {
//         unsigned i = static_cast<unsigned>(ii);
//         unsigned j = static_cast<unsigned>(jj);
//         unsigned k = q - i - j;
//         cnt++;
//         std::cout << cnt << " " << q << " " << i << " " << j << " " << k << " ";
//         std::cout << tti(-1, {i, j, k}, {a, b, c}, d) - Tetrahedron<double, double>(-1, {i, j, k}, {a, b, c}, d) << " ";
//         std::cout << tti(-1, {i, j, k}, {-a, -b, -c}, -d) - Tetrahedron<double, double>(-1, {i, j, k}, {-a, -b, -c}, -d) << " ";
//         std::cout << tti(s, {i, j, k}, {-a, -b, -c}, -d) - Tetrahedron<double, double>(s, {i, j, k}, {-a, -b, -c}, -d) << " ";
//         std::cout << tti(s, {i, j, k}, {a, b, c}, d) - Tetrahedron<double, double>(s, {i, j, k}, {a, b, c}, d) << " ";
//         std::cout << tti(0, {i, j, k}, {c, b, a}, d) - tti(0, {i, j, k}, {-c, -b, -a}, -d) - 1. / ((1 + k) * (2 + j + k) * (3 + i + j + k)) << " ";
//         std::cout << std::endl;
//       }
//     }
//   }
//   tti.printCounter();
//
//
//   std::cout << "Line \n";
//   qMax = 5;
//   s = 4;
//   HCImap <double, double> hci1(1, qMax, s);
//
//   hci1.clear();
//
//   cnt = 0;
//   for(unsigned q = 0; q <= qMax; q++) {
//     for(unsigned i = 0; i <= q; i++) {
//       std::cout << ++cnt << " " << q << " " << i << " ";
//       std::cout << hci1(-1, {i}, {a}, d) - HyperCube<double, double>(-1, {i}, {a}, d) << " ";
//       std::cout << hci1(-1, {i}, {-a}, -d) - HyperCube<double, double>(-1, {i}, {-a}, -d) << " ";
//       std::cout << hci1(s, {i}, {a}, d) - HyperCube<double, double>(s, {i}, {a}, d) << " ";
//       std::cout << hci1(s, {i}, {-a}, -d) - HyperCube<double, double>(s, {i}, {-a}, -d) << " ";
//       std::cout << hci1(0, {i}, {a}, d) + hci1(0, {i}, {-a}, -d) - 2. / (i + 1);
//       std::cout << std::endl;
//     }
//   }
//   hci1.printCounter();
//
//   cnt = 0;
//   for(unsigned q = 0; q <= qMax; q++) {
//     for(unsigned i = 0; i <= q; i++) {
//       std::cout << ++cnt << " " << q << " " << i << " ";
//       std::cout << hci1(-1, {i}, {a}, d) - HyperCube<double, double>(-1, {i}, {a}, d) << " ";
//       std::cout << hci1(-1, {i}, {-a}, -d) - HyperCube<double, double>(-1, {i}, {-a}, -d) << " ";
//       std::cout << hci1(s, {i}, {a}, d) - HyperCube<double, double>(s, {i}, {a}, d) << " ";
//       std::cout << hci1(s, {i}, {-a}, -d) - HyperCube<double, double>(s, {i}, {-a}, -d) << " ";
//       std::cout << hci1(0, {i}, {a}, d) + hci1(0, {i}, {-a}, -d) - 2. / (i + 1);
//       std::cout << std::endl;
//     }
//   }
//   hci1.printCounter();
//
//
//
//
  std::cout << "Square \n";
  qMax = 5;
  s = 4;
  SQImap <double, double> sqi(qMax, s);
  HCImap <double, double> sqr2(2, qMax, 0);
  sqi.clear();

  unsigned cnt = 0;
  for(unsigned q = 0; q <= qMax; q++) {
    for(unsigned j = 0; j <= q; j++) {
      unsigned i = q - j;
      std::cout << ++cnt << " " << q << " " << i << " " << j << " ";
      std::cout << sqi(-1, {i, j}, {a, b}, d) - sqr2(-1, {i, j}, {a, b}, d) << " ";
      std::cout << sqi(-1, {i, j}, {-a, -b}, -d) - sqr2(-1, {i, j}, {-a, -b}, -d) << " ";
      std::cout << sqi(s, {i, j}, {a, b}, d) - sqr2(s, {i, j}, {a, b}, d) << " ";
      std::cout << sqi(s, {i, j}, {-a, -b}, -d) - sqr2(s, {i, j}, {-a, -b}, -d) << " ";
      std::cout << sqi(0, {i, j}, {a, b}, d) + sqi(0, {i, j}, {-a, -b}, -d) - 4. / ((i + 1) * (j + 1));
      std::cout << std::endl;
    }
  }
  sqi.printCounter();

  cnt = 0;
  for(unsigned q = 0; q <= qMax; q++) {
    for(unsigned j = 0; j <= q; j++) {
      unsigned i = q - j;
      std::cout << ++cnt << " " << q << " " << i << " " << j << " ";
      std::cout << sqi(-1, {i, j}, {a, b}, d) - sqr2(-1, {i, j}, {a, b}, d) << " ";
      std::cout << sqi(-1, {i, j}, {-a, -b}, -d) - sqr2(-1, {i, j}, {-a, -b}, -d) << " ";
      std::cout << sqi(s, {i, j}, {a, b}, d) - sqr2(s, {i, j}, {a, b}, d) << " ";
      std::cout << sqi(s, {i, j}, {-a, -b}, -d) - sqr2(s, {i, j}, {-a, -b}, -d) << " ";
      std::cout << sqi(0, {i, j}, {a, b}, d) + sqi(0, {i, j}, {-a, -b}, -d) - 4. / ((i + 1) * (j + 1));
      std::cout << std::endl;
    }
  }


  sqi.printCounter();


//
//   std::cout << "Cube \n";
//   qMax = 3;
//   s = 4;
//   HCImap <double, double> hci3(3, qMax, s);
//
//   hci3.clear();
//
//   cnt = 0;
//   for(unsigned q = 0; q <= qMax; q++) {
//     for(int ii = q; ii >= 0; ii--) {
//       for(int jj = q - ii; jj >= 0; jj--) {
//         unsigned i = static_cast<unsigned>(ii);
//         unsigned j = static_cast<unsigned>(jj);
//         unsigned k = q - i - j;
//         std::cout << ++cnt << " " << q << " " << i << " " << j << " ";
//         std::cout << hci3(-1, {i, j, k}, {a, b, c}, d) - HyperCube<double, double>(-1, {i, j, k}, {a, b, c}, d) << " ";
//         std::cout << hci3(-1, {i, j, k}, {-a, -b, -c}, -d) - HyperCube<double, double>(-1, {i, j, k}, {-a, -b, -c}, -d) << " ";
//         std::cout << hci3(s, {i, j, k}, {a, b, c}, d) - HyperCube<double, double>(s, {i, j, k}, {a, b, c}, d) << " ";
//         std::cout << hci3(s, {i, j, k}, {-a, -b, -c}, -d) - HyperCube<double, double>(s, {i, j, k}, {-a, -b, -c}, -d) << " ";
//         std::cout << hci3(0, {i, j, k}, {a, b, c}, d) + hci3(0, {i, j, k}, {-a, -b, -c}, -d) - 8. / ((i + 1) * (j + 1) *(k + 1));
//         std::cout << std::endl;
//       }
//     }
//   }
//   hci3.printCounter();
//
//   cnt = 0;
//   for(unsigned q = 0; q <= qMax; q++) {
//     for(int ii = q; ii >= 0; ii--) {
//       for(int jj = q - ii; jj >= 0; jj--) {
//         unsigned i = static_cast<unsigned>(ii);
//         unsigned j = static_cast<unsigned>(jj);
//         unsigned k = q - i - j;
//         std::cout << ++cnt << " " << q << " " << i << " " << j << " ";
//         std::cout << hci3(-1, {i, j, k}, {a, b, c}, d) - HyperCube<double, double>(-1, {i, j, k}, {a, b, c}, d) << " ";
//         std::cout << hci3(-1, {i, j, k}, {-a, -b, -c}, -d) - HyperCube<double, double>(-1, {i, j, k}, {-a, -b, -c}, -d) << " ";
//         std::cout << hci3(s, {i, j, k}, {a, b, c}, d) - HyperCube<double, double>(s, {i, j, k}, {a, b, c}, d) << " ";
//         std::cout << hci3(s, {i, j, k}, {-a, -b, -c}, -d) - HyperCube<double, double>(s, {i, j, k}, {-a, -b, -c}, -d) << " ";
//         std::cout << hci3(0, {i, j, k}, {a, b, c}, d) + hci3(0, {i, j, k}, {-a, -b, -c}, -d) - 8. / ((i + 1) * (j + 1) *(k + 1));
//         std::cout << std::endl;
//       }
//     }
//   }
//
//   hci3.printCounter();


  return 1;






}
























