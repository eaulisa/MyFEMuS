template <class Float1>
void TestTriangle( const Float1 &eps) {
    
  typedef typename boost::math::tools::promote_args<Float1>::type myType;    
  typedef boost::multiprecision::cpp_bin_float_oct myTypeB;

  std::cout << "testing the Triangle \n";
  std::vector<unsigned>m = {6, 6};
  std::vector<myType>a = {0., 0.};
  std::vector<myTypeB> af = {0., 0.};

  myType d;
  myTypeB df;

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

    if(fabs(TriangleFull(-1,  m, a, d) - Triangle(-1,  m, af, df)) > eps) {
      std::cout << "surface test failed" << std::endl;
      std::cout << a[0] << " " << a[1] << " " << " " << d << "\n";
      std::cout << TriangleFull(-1,  m, a, d) << " " << Triangle(-1,  m, af, df) << std::endl;
    }

    if(fabs(TriangleFull(0,  m, a, d) - Triangle(0,  m, af, df)) > eps) {
      std::cout << "volume test failed" << std::endl;
      std::cout << a[0] << " " << a[1] << " " << d << "\n";
      std::cout << TriangleFull(0,  m, a, d) << " " << Triangle(0,  m, af, df) << std::endl;
    }
  }

//   for(unsigned i = 0; i < 1000; i++) {
// 
//     a[0] = -0.05;
//     a[1] = -a[0];
//     d = 500. * (0.5 * RAND_MAX - static_cast <myType>(rand())) / RAND_MAX;
// 
//     if(a[0] + a[1] + d <= 0) {
//       a[0] = -a[0];
//       a[1] = -a[1];
//       d = -d;
//     }
// 
//     af[0] = a[0];
//     af[1] = a[1];
//     df = d;
// 
//     if(fabs(Triangle(-1,  m, a, d) - TriangleFull(-1,  m, af, df)) > eps) {
//       std::cout << "surface test failed" << std::endl;
//       std::cout << a[0] << " " << a[1] << " " << " " << d << "\n";
//       std::cout << Triangle(-1,  m, a, d) << " " << TriangleFull(-1,  m, af, df) << std::endl;
//     }
// 
//     if(fabs(Triangle(0,  m, a, d) - TriangleFull(0,  m, af, df)) > eps) {
//       std::cout << "volume test failed" << std::endl;
//       std::cout << a[0] << " " << a[1] << " " << d << "\n";
//       std::cout << Triangle(0,  m, a, d) << " " << TriangleFull(0,  m, af, df) << std::endl;
//     }
//   }
// 
//   for(unsigned i = 0; i < 1000; i++) {
// 
//     a[0] = 0.;
//     a[1] = 1.;
//     d = 5. * (0.5 * RAND_MAX - static_cast <myType>(rand())) / RAND_MAX;
// 
//     if(a[0] + a[1] + d <= 0) {
//       a[0] = -a[0];
//       a[1] = -a[1];
//       d = -d;
//     }
// 
//     af[0] = a[0];
//     af[1] = a[1];
//     df = d;
// 
//     if(fabs(Triangle(-1,  m, a, d) - TriangleFull(-1,  m, af, df)) > eps) {
//       std::cout << "surface test failed" << std::endl;
//       std::cout << a[0] << " " << a[1] << " " << " " << d << "\n";
//       std::cout << Triangle(-1,  m, a, d) << " " << TriangleFull(-1,  m, af, df) << std::endl;
//     }
// 
//     if(fabs(Triangle(0,  m, a, d) - TriangleFull(0,  m, af, df)) > eps) {
//       std::cout << "volume test failed" << std::endl;
//       std::cout << a[0] << " " << a[1] << " " << d << "\n";
//       std::cout << Triangle(0,  m, a, d) << " " << TriangleFull(0,  m, af, df) << std::endl;
//     }
//   }
// 
//   for(unsigned i = 0; i < 1000; i++) {
// 
//     a[0] = 1.;
//     a[1] = 0.;
//     d = 5. * (0.5 * RAND_MAX - static_cast <myType>(rand())) / RAND_MAX;
// 
//     if(a[0] + a[1] + d <= 0) {
//       a[0] = -a[0];
//       a[1] = -a[1];
//       d = -d;
//     }
// 
//     af[0] = a[0];
//     af[1] = a[1];
//     df = d;
// 
//     if(fabs(Triangle(-1,  m, a, d) - TriangleFull(-1,  m, af, df)) > eps) {
//       std::cout << "surface test failed" << std::endl;
//       std::cout << a[0] << " " << a[1] << " " << " " << d << "\n";
//       std::cout << Triangle(-1,  m, a, d) << " " << TriangleFull(-1,  m, af, df) << std::endl;
//     }
// 
//     if(fabs(Triangle(0,  m, a, d) - TriangleFull(0,  m, af, df)) > eps) {
//       std::cout << "volume test failed" << std::endl;
//       std::cout << a[0] << " " << a[1] << " " << d << "\n";
//       std::cout << Triangle(0,  m, a, d) << " " << TriangleFull(0,  m, af, df) << std::endl;
//     }
//   }
}
