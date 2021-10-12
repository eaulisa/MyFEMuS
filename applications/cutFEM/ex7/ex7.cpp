
#include "Line.hpp"
//#include "HyperCube.hpp"

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

// #include "TestHyperCube.hpp"
// #include "TestTriangle.hpp"
#include "Tetrahedron.hpp"
// #include "TestPrism.hpp"

#include "Triangle.hpp"

int main(int, char**) {

  typedef double Type;

//     bool quad = false;
//     bool hex = false;
//     bool triangle = false;
//     bool tetrahedron = false;
//     bool prism = false;
//
//     Type eps = 1.0e-11;
//
//     if(quad) {
//         clock_t time = clock();
//         TestQuad(eps);
//         std::cout << "Time = " << static_cast<double>((clock() - time)) / CLOCKS_PER_SEC << std::endl;
//     }
//
//     if(hex) {
//         clock_t time = clock();
//         TestHex(eps);
//         std::cout << "Time = " << static_cast<double>((clock() - time)) / CLOCKS_PER_SEC << std::endl;
//     }
//
//
//     eps = 5.0e-11;
//     if(triangle) TestTriangle(eps);
//
//     if(tetrahedron) TestTetrahedron(eps);
//
//     if(prism) TestPrism(eps);


  int s = 0;
  unsigned qMax = 8;
  double a = 0.8;
  double b = 0.9;
  double c = 1;
  double d = -0.5;

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
//       std::cout << tri2(s, {i, j}, {a, b}, d) + tri(s, {i, j}, {-a, -b}, -d) - 1. / ((j + 1) * (i + j + 2)) << std::endl;
//     }
//   }

  std::cout << "Tetrahedron\n";

  TTImap <double, double> tti(qMax, 4);
  tti.clear();
  for(unsigned q = 0; q <= qMax; q++) {
    for(unsigned k = 0; k <= q; k++) {
      for(unsigned j = 0; j <= q - k ; j++) {
        unsigned i = q - j - k;
        //std::cout << q << " " << i << " " << j << " " << k << " ";
        //std::cout << tti(-1, {i, j, k}, {a, b, c}, d) << " ";
        //std::cout << tti(-1, {i, j, k}, {-a, -b, -c}, -d) << " ";
        //std::cout << tti(0, {i, j, k}, {-a, -b, -c}, -d) << " ";
        //std::cout << tti(0, {i, j, k}, {a, b, c}, d) <<" ";
        //std::cout << tti(4, {i, j, k}, {-a, -b, -c}, -d) << " ";
        //std::cout << 
        tti(0, {i, j, k}, {a, b, c}, d);// << " ";
        tti(0, {i, j, k}, {-a, -b, -c}, -d);
        tti(-1, {i, j, k}, {-a, -b, -c}, -d);
        //std::cout << std::endl;

      }
    }
  }
  tti.printCounter();

  return 1;






}























