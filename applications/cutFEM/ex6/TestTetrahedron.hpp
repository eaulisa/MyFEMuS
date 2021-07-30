template <class Float1>
void TestTetrahedron( const Float1 &eps) {
    
  typedef typename boost::math::tools::promote_args<Float1>::type myType;    
  typedef boost::multiprecision::cpp_bin_float_oct myTypeB;

  std::cout << "testing the Tetrahedron \n";
  std::vector<unsigned>m = {6, 6, 6};
  std::vector<myType>a = {0., 0., 0.};
  std::vector<myTypeB> af = {0., 0., 0.};


  myType d;
  myTypeB df;

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

    if(fabs(Tetrahedron(-1,  m, a, d) - TetrahedronB(-1,  m, af, df)) > eps) {
      std::cout << "surface test failed" << std::endl;
      std::cout << a[0] << " " << a[1] << " " << a[2] << " " << d << "\n";
      std::cout << Tetrahedron(-1,  m, a, d) << " " << TetrahedronB(-1,  m, af, df) << std::endl;
    }

    if(fabs(Tetrahedron(0,  m, a, d) - TetrahedronB(0,  m, af, df)) > eps) {
      std::cout << "volume test failed" << std::endl;
      std::cout << a[0] << " " << a[1] << " " << a[2] << " " << d << "\n";
      std::cout << Tetrahedron(0,  m, a, d) << " " << TetrahedronB(0,  m, af, df) << std::endl;
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

    if(fabs(Tetrahedron(-1,  m, a, d) - TetrahedronB(-1,  m, af, df)) > eps) {
      std::cout << "surface test failed" << std::endl;
      std::cout << a[0] << " " << a[1] << " " << a[2] << " " << d << "\n";
      std::cout << Tetrahedron(-1,  m, a, d) << " " << TetrahedronB(-1,  m, af, df) << std::endl;
    }

    if(fabs(Tetrahedron(0,  m, a, d) - TetrahedronB(0,  m, af, df)) > eps) {
      std::cout << "volume test failed" << std::endl;
      std::cout << a[0] << " " << a[1] << " " << a[2] << " " << d << "\n";
      std::cout << Tetrahedron(0,  m, a, d) << " " << TetrahedronB(0,  m, af, df) << std::endl;
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

    if(fabs(Tetrahedron(-1,  m, a, d) - TetrahedronB(-1,  m, af, df)) > eps) {
      std::cout << "surface test failed" << std::endl;
      std::cout << a[0] << " " << a[1] << " " << a[2] << " " << d << "\n";
      std::cout << Tetrahedron(-1,  m, a, d) << " " << TetrahedronB(-1,  m, af, df) << std::endl;
    }

    if(fabs(Tetrahedron(0,  m, a, d) - TetrahedronB(0,  m, af, df)) > eps) {
      std::cout << "volume test failed" << std::endl;
      std::cout << a[0] << " " << a[1] << " " << a[2] << " " << d << "\n";
      std::cout << Tetrahedron(0,  m, a, d) << " " << TetrahedronB(0,  m, af, df) << std::endl;
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

    if(fabs(Tetrahedron(-1,  m, a, d) - TetrahedronB(-1,  m, af, df)) > eps) {
      std::cout << "surface test failed" << std::endl;
      std::cout << a[0] << " " << a[1] << " " << a[2] << " " << d << "\n";
      std::cout << Tetrahedron(-1,  m, a, d) << " " << TetrahedronB(-1,  m, af, df) << std::endl;
    }

    if(fabs(Tetrahedron(0,  m, a, d) - TetrahedronB(0,  m, af, df)) > eps) {
      std::cout << "volume test failed" << std::endl;
      std::cout << a[0] << " " << a[1] << " " << a[2] << " " << d << "\n";
      std::cout << Tetrahedron(0,  m, a, d) << " " << TetrahedronB(0,  m, af, df) << std::endl;
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

    if(fabs(Tetrahedron(-1,  m, a, d) - TetrahedronB(-1,  m, af, df)) > eps) {
      std::cout << "surface test failed" << std::endl;
      std::cout << a[0] << " " << a[1] << " " << a[2] << " " << d << "\n";
      std::cout << Tetrahedron(-1,  m, a, d) << " " << TetrahedronB(-1,  m, af, df) << std::endl;
    }

    if(fabs(Tetrahedron(0,  m, a, d) - TetrahedronB(0,  m, af, df)) > eps) {
      std::cout << "volume test failed" << std::endl;
      std::cout << a[0] << " " << a[1] << " " << a[2] << " " << d << "\n";
      std::cout << Tetrahedron(0,  m, a, d) << " " << TetrahedronB(0,  m, af, df) << std::endl;
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

    if(fabs(Tetrahedron(-1,  m, a, d) - TetrahedronB(-1,  m, af, df)) > eps) {
      std::cout << "surface test failed" << std::endl;
      std::cout << a[0] << " " << a[1] << " " << a[2] << " " << d << "\n";
      std::cout << Tetrahedron(-1,  m, a, d) << " " << TetrahedronB(-1,  m, af, df) << std::endl;
      std::cout << Tetrahedron(-1,  m, af, df) << " " << TetrahedronB(-1,  m, a, d) << std::endl;
    }

    if(fabs(Tetrahedron(0,  m, a, d) - TetrahedronB(0,  m, af, df)) > eps) {
      std::cout << "volume test failed" << std::endl;
      std::cout << a[0] << " " << a[1] << " " << a[2] << " " << d << "\n";
      std::cout << Tetrahedron(0,  m, a, d) << " " << TetrahedronB(0,  m, af, df) << std::endl;
      std::cout << Tetrahedron(0,  m, af, df) << " " << TetrahedronB(0,  m, a, d) << std::endl;
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

    if(fabs(Tetrahedron(-1,  m, a, d) - TetrahedronB(-1,  m, af, df)) > eps) {
      std::cout << "surface test failed" << std::endl;
      std::cout << a[0] << " " << a[1] << " " << a[2] << " " << d << "\n";
      std::cout << Tetrahedron(-1,  m, a, d) << " " << TetrahedronB(-1,  m, af, df) << std::endl;
    }

    if(fabs(Tetrahedron(0,  m, a, d) - TetrahedronB(0,  m, af, df)) > eps) {
      std::cout << "volume test failed" << std::endl;
      std::cout << a[0] << " " << a[1] << " " << a[2] << " " << d << "\n";
      std::cout << Tetrahedron(0,  m, a, d) << " " << TetrahedronB(0,  m, af, df) << std::endl;
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

    if(fabs(Tetrahedron(-1,  m, a, d) - TetrahedronB(-1,  m, af, df)) > eps) {
      std::cout << "surface test failed" << std::endl;
      std::cout << a[0] << " " << a[1] << " " << a[2] << " " << d << "\n";
      std::cout << Tetrahedron(-1,  m, a, d) << " " << TetrahedronB(-1,  m, af, df) << std::endl;
    }

    if(fabs(Tetrahedron(0,  m, a, d) - TetrahedronB(0,  m, af, df)) > eps) {
      std::cout << "volume test failed" << std::endl;
      std::cout << a[0] << " " << a[1] << " " << a[2] << " " << d << "\n";
      std::cout << Tetrahedron(0,  m, a, d) << " " << TetrahedronB(0,  m, af, df) << std::endl;
    }
  }
}


