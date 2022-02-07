
#include "GeomElTypeEnum.hpp"

#include "Line.hpp"
#include "Square.hpp"
#include "Cube.hpp"
#include "HyperCube.hpp"
#include "Triangle.hpp"
#include "Tetrahedron.hpp"
#include "Prism.hpp" 


#include <iostream>
#include <boost/multiprecision/cpp_bin_float.hpp>

#include "CutFemIntegration.hpp"


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

#include "GramSchmidt.hpp"

#include "TestHyperCube.hpp"
// #include "TestTriangle.hpp"
// #include "TestPrism.hpp"



int main(int, char**) {

  typedef double Type;

  typedef cpp_bin_float_oct Type2;

  bool quad = true;
  bool hex = false;
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


    SQImap <double, double> sqr(qMax, 0);

    for(unsigned tt = 0; tt < 1; tt++) {
      sqr.clear();
      a = -1. + 2. * rand() / RAND_MAX;
      b = -1. + 2. * rand() / RAND_MAX;
      d = -1. + 2. * rand() / RAND_MAX;
      for(unsigned q = 0; q <= qMax; q++) {
        for(unsigned j = 0; j <= q; j++) {
          unsigned i = q - j;
          std::cout << i << " " << j << " " << sqr(-1, {i, j}, {a, b}, d) << " ";
          std::cout << sqr(0, {i, j}, {a, b}, d) << std::endl;
        }
      }
    }
    sqr.printCounter();
    std::cout << "Map Time = " << static_cast<double>((clock() - time)) / CLOCKS_PER_SEC << std::endl;

//     HCImap <double, double> sqr2(2, qMax, 0);
//     time = clock();
//     for(unsigned tt = 0; tt < 10000; tt++) {
//       sqr2.clear();
//       a = -1. + 2. * rand() / RAND_MAX;
//       b = -1. + 2. * rand() / RAND_MAX;
//       d = -1. + 2. * rand() / RAND_MAX;
//       for(unsigned q = 0; q <= qMax; q++) {
//         for(unsigned j = 0; j <= q; j++) {
//           unsigned i = q - j;
//           sqr2(-1, {i, j}, {a, b}, d);
//           sqr2(0, {i, j}, {a, b}, d);
//         }
//       }
//     }
//     sqr2.printCounter();
//     std::cout << "Old Time = " << static_cast<double>((clock() - time)) / CLOCKS_PER_SEC << std::endl;
  }



  if(hex) {
    clock_t time;
    TestHex(eps);
    //TestHexOld(eps);


    CBImap <double, double> cube(qMax, 0);
    srand(0);
    time = clock();

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
    cube.clear();


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
    hci3.clear();




  }

  /*
    std::vector<std::vector<Type2>> ATA;
    std::cout << "Hex" << std::endl;
    Get_GS_ATA_Matrix(HEX, 2, ATA, false);
    std::cout << std::endl;
    for(unsigned i = 0; i < ATA.size(); i++) {
      for(unsigned j = 0; j < ATA[i].size(); j++) {
        std::cout << ATA[i][j] << " ";
      }
      std::cout << std::endl;
    }


    std::cout << "Tet" << std::endl;
    Get_GS_ATA_Matrix(TET, 2, ATA, false);
    std::cout << std::endl;
    for(unsigned i = 0; i < ATA.size(); i++) {
      for(unsigned j = 0; j < ATA[i].size(); j++) {
        std::cout << ATA[i][j] << " ";
      }
      std::cout << std::endl;
    }


    std::cout << "Wedge" << std::endl;
    Get_GS_ATA_Matrix(WEDGE, 2, ATA, false);
    std::cout << std::endl;
    for(unsigned i = 0; i < ATA.size(); i++) {
      for(unsigned j = 0; j < ATA[i].size(); j++) {
        std::cout << ATA[i][j] << " ";
      }
      std::cout << std::endl;
    }

    std::cout << "Quad" << std::endl;
    Get_GS_ATA_Matrix(QUAD, 2, ATA, false);
    std::cout << std::endl;
    for(unsigned i = 0; i < ATA.size(); i++) {
      for(unsigned j = 0; j < ATA[i].size(); j++) {
        std::cout << ATA[i][j] << " ";
      }
      std::cout << std::endl;
    }

    std::cout << "Tri" << std::endl;
    Get_GS_ATA_Matrix(TRI, 2, ATA, false);
    std::cout << std::endl;
    for(unsigned i = 0; i < ATA.size(); i++) {
      for(unsigned j = 0; j < ATA[i].size(); j++) {
        std::cout << ATA[i][j] << " ";
      }
      std::cout << std::endl;
    }*/

  {

    unsigned qM = 6;
    HCImap <Type2, Type2> hci1(1, qM, 0);
    std::vector< std::vector<Type2> > f(1, std::vector<Type2>(qM + 1));

    std::cout << "line" << std::endl;

    for(unsigned i = 0; i <= qM; i++) {
      f[0][i] = hci1(0, {i}, {1.}, 0);
    }

    std::vector<std::vector<Type2>> ATA;
    Get_GS_ATA_Matrix(LINE, qM, ATA, false);
    std::vector< std::vector<Type2> > Co = MatrixMatrixMultiply(f, ATA);


    //print
    std::cout.precision(20);
    std::cout << "fo = ";
    for(unsigned i = 0; i < f[0].size(); i++) {
      std::cout << f[0][i] << " ";
    }
    std::cout << std::endl;

    std::cout << "Co = ";
    for(unsigned i = 0; i < Co.size(); i++) {
      for(unsigned j = 0; j < Co[i].size(); j++) {
        std::cout << Co[i][j] << " ";
      }
      std::cout << std::endl;
    }
    std::cout << std::endl;

  }
  {
    unsigned qM = 3;
    HCImap <Type2, Type2> hci2(2, qM, 0);
    std::vector< std::vector<Type2> > f(1, std::vector<Type2>((qM + 1) * (qM + 2) / 2));

    std::cout << "quad" << std::endl;

    unsigned count = 0;
    for(unsigned q = 0; q <= qM; q++) {
      for(unsigned j = 0; j <= q; j++) {
        unsigned i = q - j;
        f[0][count] = hci2(0, {i, j}, {1., 1.}, 0);
        count++;
      }
    }

    std::vector<std::vector<Type2>> ATA;
    Get_GS_ATA_Matrix(QUAD, qM, ATA, false);

    std::vector< std::vector<Type2> > Co = MatrixMatrixMultiply(f, ATA);

    //print
    std::cout.precision(20);
    std::cout << "fo = ";
    for(unsigned i = 0; i < f[0].size(); i++) {
      std::cout << f[0][i] << " ";
    }
    std::cout << std::endl;

    std::cout << "Co = ";
    for(unsigned i = 0; i < Co.size(); i++) {
      for(unsigned j = 0; j < Co[i].size(); j++) {
        std::cout << Co[i][j] << " ";
      }
      std::cout << std::endl;
    }
    std::cout << std::endl;
  }


  {
    unsigned qM = 3;
    TRImap <Type2, Type2> tri(2, qM, 0);
    std::vector< std::vector<Type2> > f(1, std::vector<Type2>((qM + 1) * (qM + 2) / 2));

    std::cout << "tri" << std::endl;

    unsigned count = 0;
    for(unsigned q = 0; q <= qM; q++) {
      for(unsigned j = 0; j <= q; j++) {
        unsigned i = q - j;
        f[0][count] = tri(-1, {i, j}, {-1 / sqrt(2), 1 / sqrt(2)}, 0);
        count++;
      }
    }

    std::vector<std::vector<Type2>> ATA;
    Get_GS_ATA_Matrix(TRI, qM, ATA, false);

    std::vector< std::vector<Type2> > Co = MatrixMatrixMultiply(f, ATA);

    //print
    std::cout.precision(20);
    std::cout << "fo = ";
    for(unsigned i = 0; i < f[0].size(); i++) {
      std::cout << f[0][i] << " ";
    }
    std::cout << std::endl;

    std::cout << "Co = ";
    for(unsigned i = 0; i < Co.size(); i++) {
      for(unsigned j = 0; j < Co[i].size(); j++) {
        std::cout << Co[i][j] << " ";
      }
      std::cout << std::endl;
    }
    std::cout << std::endl;
  }

  {
    unsigned qM = 3;
    HCImap <Type2, Type2> hci3(3, qM, 0);
    std::vector< std::vector<Type2> > f(1, std::vector<Type2>((qM + 1) * (qM + 2) * (qM + 3) / 6));

    std::cout << "hex" << std::endl;

    unsigned count = 0;

    for(unsigned q = 0; q <= qM; q++) {
      for(int ii = q; ii >= 0; ii--) {
        for(int jj = q - ii; jj >= 0; jj--) {
          unsigned i = static_cast<unsigned>(ii);
          unsigned j = static_cast<unsigned>(jj);
          unsigned k = q - i - j;
          f[0][count] = hci3(-1, {i, j, k}, {0.1/sqrt(1.02), 0.1/sqrt(1.02), -1./sqrt(1.02)}, 0);
//           std::cout << count << " " << i << " " << j << " " << k << std::endl << std::flush;
          count++;

        }
      }
    }


    std::vector<std::vector<Type2>> ATA;
    Get_GS_ATA_Matrix(HEX, qM, ATA, false);

    std::vector< std::vector<Type2> > Co = MatrixMatrixMultiply(f, ATA);

    //print
    std::cout.precision(20);
    std::cout << "fo = ";
    for(unsigned i = 0; i < f[0].size(); i++) {
      std::cout << f[0][i] << " ";
    }
    std::cout << std::endl;

    std::cout << "Co = ";
    for(unsigned i = 0; i < Co.size(); i++) {
      for(unsigned j = 0; j < Co[i].size(); j++) {
        std::cout << "C"<<j<<" = "<<((fabs(Co[i][j]) < 1.0e-60) ? 0. : Co[i][j]) << "; ";
      }
      std::cout << std::endl;
    }
    std::cout << std::endl;
  }
  
  {
    unsigned qM = 3;
    HCImap <Type2, Type2> hci3(3, qM, 0);
    std::vector< std::vector<Type2> > f(1, std::vector<Type2>((qM + 1) * (qM + 2) * (qM + 3) / 6));

    std::cout << "pri" << std::endl;

    unsigned count = 0;

    for(unsigned q = 0; q <= qM; q++) {
      for(int ii = q; ii >= 0; ii--) {
        for(int jj = q - ii; jj >= 0; jj--) {
          unsigned i = static_cast<unsigned>(ii);
          unsigned j = static_cast<unsigned>(jj);
          unsigned k = q - i - j;
          f[0][count] = Prism<Type2, Type2>(-1, {i, j, k}, {-0.1/sqrt(1.02), 0.1/sqrt(1.02), 1./sqrt(1.02)}, 0);
//           std::cout << count << " " << i << " " << j << " " << k << std::endl << std::flush;
          count++;
        }
      }
    }

    std::vector<std::vector<Type2>> ATA;
    Get_GS_ATA_Matrix(WEDGE, qM, ATA, false);

    std::vector< std::vector<Type2> > Co = MatrixMatrixMultiply(f, ATA);

    //print
    std::cout.precision(20);
    std::cout << "fo = ";
    for(unsigned i = 0; i < f[0].size(); i++) {
      std::cout << f[0][i] << " ";
    }
    std::cout << std::endl;

    std::cout << "Co = ";
    for(unsigned i = 0; i < Co.size(); i++) {
      for(unsigned j = 0; j < Co[i].size(); j++) {
        std::cout << "C"<<j<<" = "<<((fabs(Co[i][j]) < 1.0e-60) ? 0. : Co[i][j]) << "; ";
      }
      std::cout << std::endl;
    }
    std::cout << std::endl;
  }
  
  {
    unsigned qM = 3;
    CutFEMmap <Type2> *tet3p  = new TTImap <Type2, Type2> (3, qM, 0);
//     TTImap <Type2, Type2> tet3(3, qM, 0);
    std::vector< std::vector<Type2> > f(1, std::vector<Type2>((qM + 1) * (qM + 2) * (qM + 3) / 6));

    std::cout << "tet" << std::endl;

    unsigned count = 0;

    for(unsigned q = 0; q <= qM; q++) {
      for(int ii = q; ii >= 0; ii--) {
        for(int jj = q - ii; jj >= 0; jj--) {
          unsigned i = static_cast<unsigned>(ii);
          unsigned j = static_cast<unsigned>(jj);
          unsigned k = q - i - j;
//           f[0][count] = tet3(-1, {i, j, k}, {0/sqrt(5), 2/sqrt(5), -1/sqrt(5)}, 0);
//           f[0][count] = (*tet3p)(-1, {i, j, k}, {0/sqrt(5), 2/sqrt(5), -1/sqrt(5)}, 0);
          f[0][count] = (*tet3p)(0, {i, j, k}, {-1,-1,-2}, 1);
//           f[0][count] = tet3(0, {i, j, k}, {0.1, 1, -0.5}, 0.05);
          count++;
        }
      }
    }

    std::vector<std::vector<Type2>> ATA;
    Get_GS_ATA_Matrix(TET, qM, ATA, false);

    std::vector< std::vector<Type2> > Co = MatrixMatrixMultiply(f, ATA);

    //print
    std::cout.precision(20);
    std::cout << "fo = ";
    for(unsigned i = 0; i < f[0].size(); i++) {
      std::cout << f[0][i] << " ";
    }
    std::cout << std::endl;

    std::cout << "Co = ";
    for(unsigned i = 0; i < Co.size(); i++) {
      for(unsigned j = 0; j < Co[i].size(); j++) {
        std::cout << "C"<<j<<" = "<<((fabs(Co[i][j]) < 1.0e-60) ? 0. : Co[i][j]) << "; ";
      }
      std::cout << std::endl;
    }
    std::cout << std::endl;
    
    delete tet3p;
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






























