

void TestTetrahedron() {

  std::cout.precision(16);

  typedef double myTypeA;
  typedef long double myTypeB;
  //typedef boost::multiprecision::cpp_bin_float_oct myTypeB;

  std::cout << "testing the Tetrahedron \n";
  std::vector<unsigned>m = {3, 3, 3};


  std::vector<myTypeB> af = {0.0, 0.0, -1.};
  std::vector<myTypeA> a(3);

  a[0] = static_cast<myTypeA>(af[0]);
  a[1] = static_cast<myTypeA>(af[1]);
  a[2] = static_cast<myTypeA>(af[2]);

  myTypeB df = 0.00005;
  myTypeA d = static_cast<myTypeA>(df);

  int s = 0;
  std::cout << "Epsilon cases\n";
  std::cout << "s = " << s << " a = " << a[0] << " b = " << a[1] << " c = " << a[2] << " d = " << d << std::endl;
  std::cout << "s = " << s << " af= " << af[0] << " bf= " << af[1] << " cf= " << af[2] << " df= " << df << std::endl;
  std::cout << TetrahedronA(s,  m, a, d) << std::endl;
  std::cout << TetrahedronA(s,  m, af, df) << std::endl;
  std::cout << (TetrahedronA(s,  m, af, df) - TetrahedronA(s, m, a, d)) / TetrahedronA(s, m, a, d) << std::endl;
}



template <class Float1>
void TestTetrahedron(const Float1 &eps) {

  typedef typename boost::math::tools::promote_args<Float1>::type myType;
  typedef boost::multiprecision::cpp_bin_float_oct myTypeB;

 



  std::cout << "testing the Tetrahedron \n";
  std::vector<unsigned>m = {7, 5, 6};

  myType scale = 1;

  std::vector<myType> a = {0.00000 * scale, 0.000000 * scale, -1 * scale};

  std::vector<myTypeB> af(3);

  myType d = 0.1 * scale;
  myTypeB df;

  af[0] = a[0];
  af[1] = a[1];
  af[2] = a[2];
  df = d;

  std::cout.precision(16);

  int s = 4;
  std::cout << "Epsilon cases\n";
  std::cout << "s = " << s << " a = " << a[0] << " b = " << a[1] << " c = " << a[2] << " d = " << d << std::endl;
  std::cout << "s = " << s << " af= " << af[0] << " bf= " << af[1] << " cf= " << af[2] << " df= " << df << std::endl;
  myType I1 = TetrahedronA(s,  m, a, d);
  std::cout << I1 << std::endl;
  myTypeB I2 = TetrahedronA(s,  m, af, df);
  std::cout << I2 << std::endl;
  std::cout << (I1 - I2) / I2 << std::endl;

  return;

  for(unsigned i = 0; i < 1000; i++) {
    a[0] = (0.5 * RAND_MAX - static_cast <myType>(rand())) / RAND_MAX;
    a[1] = (0.5 * RAND_MAX - static_cast <myType>(rand())) / RAND_MAX;
    a[2] = sqrt(1 - a[0] * a[0] + a[1] * a[1]);
    d = 10. * (0.5 * RAND_MAX - static_cast <myType>(rand())) / RAND_MAX;

    if(a[0] + a[1] + d <= 0) {
      a[0] = -a[0];
      a[1] = -a[1];
      a[2] = -a[2];
      d = -d;
    }

    af[0] = a[0];
    af[1] = a[1];
    af[2] = a[2];
    df = d;

    if(fabs(TetrahedronA(-1,  m, a, d) - TetrahedronB(-1,  m, af, df)) > eps) {
      std::cout << "surface test failed" << std::endl;
      std::cout << a[0] << " " << a[1] << " " << a[2] << " " << d << "\n";
      std::cout << TetrahedronA(-1,  m, a, d) << " " << TetrahedronB(-1,  m, af, df) << std::endl;
    }

    if(fabs(TetrahedronA(0,  m, a, d) - TetrahedronB(0,  m, af, df)) > eps) {
      std::cout << "volume test failed" << std::endl;
      std::cout << a[0] << " " << a[1] << " " << a[2] << " " << d << "\n";
      std::cout << TetrahedronA(0,  m, a, d) << " " << TetrahedronB(0,  m, af, df) << std::endl;
    }
  }

  for(unsigned i = 0; i < 1000; i++) {

    a[0] = (0.5 * RAND_MAX - static_cast <myType>(rand())) / RAND_MAX;
    a[1] = -a[0];
    a[2] = sqrt(1 - a[0] * a[0] + a[1] * a[1]);
    d = 10. * (0.5 * RAND_MAX - static_cast <myType>(rand())) / RAND_MAX;

    if(a[0] + a[1] + d <= 0) {
      a[0] = -a[0];
      a[1] = -a[1];
      a[2] = -a[2];
      d = -d;
    }

    af[0] = a[0];
    af[1] = a[1];
    af[2] = a[2];
    df = d;

    if(fabs(TetrahedronA(-1,  m, a, d) - TetrahedronB(-1,  m, af, df)) > eps) {
      std::cout << "surface test failed" << std::endl;
      std::cout << a[0] << " " << a[1] << " " << a[2] << " " << d << "\n";
      std::cout << TetrahedronA(-1,  m, a, d) << " " << TetrahedronB(-1,  m, af, df) << std::endl;
    }

    if(fabs(TetrahedronA(0,  m, a, d) - TetrahedronB(0,  m, af, df)) > eps) {
      std::cout << "volume test failed" << std::endl;
      std::cout << a[0] << " " << a[1] << " " << a[2] << " " << d << "\n";
      std::cout << TetrahedronA(0,  m, a, d) << " " << TetrahedronB(0,  m, af, df) << std::endl;
    }
  }

  for(unsigned i = 0; i < 1000; i++) {
    a[0] = 0;
    a[1] = 0;
    a[2] = 1;
    d = 10. * (0.5 * RAND_MAX - static_cast <myType>(rand())) / RAND_MAX;

    if(a[0] + a[1] + d <= 0) {
      a[0] = -a[0];
      a[1] = -a[1];
      a[2] = -a[2];
      d = -d;
    }

    af[0] = a[0];
    af[1] = a[1];
    af[2] = a[2];
    df = d;

    if(fabs(TetrahedronA(-1,  m, a, d) - TetrahedronB(-1,  m, af, df)) > eps) {
      std::cout << "surface test failed" << std::endl;
      std::cout << a[0] << " " << a[1] << " " << a[2] << " " << d << "\n";
      std::cout << TetrahedronA(-1,  m, a, d) << " " << TetrahedronB(-1,  m, af, df) << std::endl;
    }

    if(fabs(TetrahedronA(0,  m, a, d) - TetrahedronB(0,  m, af, df)) > eps) {
      std::cout << "volume test failed" << std::endl;
      std::cout << a[0] << " " << a[1] << " " << a[2] << " " << d << "\n";
      std::cout << TetrahedronA(0,  m, a, d) << " " << TetrahedronB(0,  m, af, df) << std::endl;
    }
  }

  for(unsigned i = 0; i < 1000; i++) {
    a[0] = 0.;
    a[1] = 1.;
    a[2] = 0.;
    d = 10. * (0.5 * RAND_MAX - static_cast <myType>(rand())) / RAND_MAX;

    if(a[0] + a[1] + d <= 0) {
      a[0] = -a[0];
      a[1] = -a[1];
      a[2] = -a[2];
      d = -d;
    }

    af[0] = a[0];
    af[1] = a[1];
    af[2] = a[2];
    df = d;

    if(fabs(TetrahedronA(-1,  m, a, d) - TetrahedronB(-1,  m, af, df)) > eps) {
      std::cout << "surface test failed" << std::endl;
      std::cout << a[0] << " " << a[1] << " " << a[2] << " " << d << "\n";
      std::cout << TetrahedronA(-1,  m, a, d) << " " << TetrahedronB(-1,  m, af, df) << std::endl;
    }

    if(fabs(TetrahedronA(0,  m, a, d) - TetrahedronB(0,  m, af, df)) > eps) {
      std::cout << "volume test failed" << std::endl;
      std::cout << a[0] << " " << a[1] << " " << a[2] << " " << d << "\n";
      std::cout << TetrahedronA(0,  m, a, d) << " " << TetrahedronB(0,  m, af, df) << std::endl;
    }
  }

  for(unsigned i = 0; i < 1000; i++) {
    a[0] = 1.;
    a[1] = 0.;
    a[2] = 0.;

    d = 10. * (0.5 * RAND_MAX - static_cast <myType>(rand())) / RAND_MAX;

    if(a[0] + a[1] + d <= 0) {
      a[0] = -a[0];
      a[1] = -a[1];
      a[2] = -a[2];
      d = -d;
    }

    af[0] = a[0];
    af[1] = a[1];
    af[2] = a[2];
    df = d;

    if(fabs(TetrahedronA(-1,  m, a, d) - TetrahedronB(-1,  m, af, df)) > eps) {
      std::cout << "surface test failed" << std::endl;
      std::cout << a[0] << " " << a[1] << " " << a[2] << " " << d << "\n";
      std::cout << TetrahedronA(-1,  m, a, d) << " " << TetrahedronB(-1,  m, af, df) << std::endl;
    }

    if(fabs(TetrahedronA(0,  m, a, d) - TetrahedronB(0,  m, af, df)) > eps) {
      std::cout << "volume test failed" << std::endl;
      std::cout << a[0] << " " << a[1] << " " << a[2] << " " << d << "\n";
      std::cout << TetrahedronA(0,  m, a, d) << " " << TetrahedronB(0,  m, af, df) << std::endl;
    }
  }


  for(unsigned i = 0; i < 1000; i++) {

    a[0] = (0.5 * RAND_MAX - static_cast <myType>(rand())) / RAND_MAX;
    a[1] = sqrt(1. - a[0] * a[0]);
    a[2] = 0.;
    d = 10. * (0.5 * RAND_MAX - static_cast <myType>(rand())) / RAND_MAX;

    if(a[0] + a[1] + d <= 0) {
      a[0] = -a[0];
      a[1] = -a[1];
      a[2] = -a[2];
      d = -d;
    }

    af[0] = a[0];
    af[1] = a[1];
    af[2] = a[2];
    df = d;

    if(fabs(TetrahedronA(-1,  m, a, d) - TetrahedronB(-1,  m, af, df)) > eps) {
      std::cout << "surface test failed" << std::endl;
      std::cout << a[0] << " " << a[1] << " " << a[2] << " " << d << "\n";
      std::cout << TetrahedronA(-1,  m, a, d) << " " << TetrahedronB(-1,  m, af, df) << std::endl;
      std::cout << TetrahedronA(-1,  m, af, df) << " " << TetrahedronB(-1,  m, a, d) << std::endl;
    }

    if(fabs(TetrahedronA(0,  m, a, d) - TetrahedronB(0,  m, af, df)) > eps) {
      std::cout << "volume test failed" << std::endl;
      std::cout << a[0] << " " << a[1] << " " << a[2] << " " << d << "\n";
      std::cout << TetrahedronA(0,  m, a, d) << " " << TetrahedronB(0,  m, af, df) << std::endl;
      std::cout << TetrahedronA(0,  m, af, df) << " " << TetrahedronB(0,  m, a, d) << std::endl;
    }
  }

  for(unsigned i = 0; i < 1000; i++) {

    a[0] = (0.5 * RAND_MAX - static_cast <myType>(rand())) / RAND_MAX;
    a[1] = 0.;
    a[2] = sqrt(1. - a[0] * a[0]);;
    d = 10. * (0.5 * RAND_MAX - static_cast <myType>(rand())) / RAND_MAX;

    if(a[0] + a[1] + d <= 0) {
      a[0] = -a[0];
      a[1] = -a[1];
      a[2] = -a[2];
      d = -d;
    }

    af[0] = a[0];
    af[1] = a[1];
    af[2] = a[2];
    df = d;

    if(fabs(TetrahedronA(-1,  m, a, d) - TetrahedronB(-1,  m, af, df)) > eps) {
      std::cout << "surface test failed" << std::endl;
      std::cout << a[0] << " " << a[1] << " " << a[2] << " " << d << "\n";
      std::cout << TetrahedronA(-1,  m, a, d) << " " << TetrahedronB(-1,  m, af, df) << std::endl;
    }

    if(fabs(TetrahedronA(0,  m, a, d) - TetrahedronB(0,  m, af, df)) > eps) {
      std::cout << "volume test failed" << std::endl;
      std::cout << a[0] << " " << a[1] << " " << a[2] << " " << d << "\n";
      std::cout << TetrahedronA(0,  m, a, d) << " " << TetrahedronB(0,  m, af, df) << std::endl;
    }
  }


  for(unsigned i = 0; i < 1000; i++) {

    a[0] = 0.;
    a[1] = (0.5 * RAND_MAX - static_cast <myType>(rand())) / RAND_MAX;
    a[2] = sqrt(1. - a[0] * a[0]);;
    d = 10. * (0.5 * RAND_MAX - static_cast <myType>(rand())) / RAND_MAX;

    if(a[0] + a[1] + d <= 0) {
      a[0] = -a[0];
      a[1] = -a[1];
      a[2] = -a[2];
      d = -d;
    }

    af[0] = a[0];
    af[1] = a[1];
    af[2] = a[2];
    df = d;

    if(fabs(TetrahedronA(-1,  m, a, d) - TetrahedronB(-1,  m, af, df)) > eps) {
      std::cout << "surface test failed" << std::endl;
      std::cout << a[0] << " " << a[1] << " " << a[2] << " " << d << "\n";
      std::cout << TetrahedronA(-1,  m, a, d) << " " << TetrahedronB(-1,  m, af, df) << std::endl;
    }

    if(fabs(TetrahedronA(0,  m, a, d) - TetrahedronB(0,  m, af, df)) > eps) {
      std::cout << "volume test failed" << std::endl;
      std::cout << a[0] << " " << a[1] << " " << a[2] << " " << d << "\n";
      std::cout << TetrahedronA(0,  m, a, d) << " " << TetrahedronB(0,  m, af, df) << std::endl;
    }
  }
}


