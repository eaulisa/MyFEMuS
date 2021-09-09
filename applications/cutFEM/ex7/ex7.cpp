
#include "Line.hpp"
#include "HyperCube.hpp"

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

#include "TestHyperCube.hpp"
#include "TestTriangle.hpp"
#include "TestTetrahedron.hpp"

int main(int, char**) {

  typedef double Type;

  bool quad = false;//false
  bool hex = false;//false
  bool triangle = false;//false
  bool tetrahedron = true;//false

  Type eps = 1.0e-11;

  if(quad) {
    clock_t time = clock();
    TestQuad(eps);
    std::cout << "Time = " << static_cast<double>((clock() - time)) / CLOCKS_PER_SEC << std::endl;
  }

  if(hex) {
    clock_t time = clock();
    TestHex(eps);
    std::cout << "Time = " << static_cast<double>((clock() - time)) / CLOCKS_PER_SEC << std::endl;
  }
  
  
  eps = 5.0e-11;
  if(triangle) TestTriangle(eps);

  
  //cpp_bin_float_quad deps = 5.0e-10;
  double deps = 5.0e-10;
  if(tetrahedron) TestTetrahedron(deps);
    
  
  return 1;






}























