
template <class Float1>
void TestHexahedron(const Float1 &eps) {

  typedef typename boost::math::tools::promote_args<Float1>::type myType;
  typedef boost::multiprecision::cpp_bin_float_oct myTypeB;


  std::cout << "testing the Hexahedron \n";
  std::vector<unsigned>m = {6, 6, 6};
  std::vector<myType>a = {0., 0., 0.};
  std::vector<myTypeB> af = {0., 0., 0.};

  myType d;
  myTypeB df;

  for(unsigned i = 0; i < 1000; i++) {
    a[0] = (0.5 * RAND_MAX - static_cast <myType>(rand())) / RAND_MAX;
    a[1] = (0.5 * RAND_MAX - static_cast <myType>(rand())) / RAND_MAX;
    a[2] = (0.5 * RAND_MAX - static_cast <myType>(rand())) / RAND_MAX;

    myType det = sqrt(a[0] * a[0] + a[1] * a[1] + a[2] * a[2]);
    a[0] /= det;
    a[1] /= det;
    a[2] /= det;

    d = 5. * (0.5 * RAND_MAX - static_cast <myType>(rand())) / RAND_MAX;

    af[0] = a[0];
    af[1] = a[1];
    af[2] = a[2];
    df = d;

    if(fabs(HyperCube(-1, m, a, d) - HyperCubeAold(-1, m, af, df)) > eps) {
      std::cout << "surface test failed" << std::endl;
      std::cout << a[0] << " " << a[1] << " " << a[2] << " " << d << "\n";
      std::cout << HyperCube(-1, m, a, d) << " " << HyperCubeAold(-1, m, af, df) << std::endl;
    }

    if(fabs(HyperCube(0, m, a, d) - HyperCubeAold(0, m, af, df)) > eps) {
      std::cout << "volume test failed" << std::endl;
      std::cout << a[0] << " " << a[1] << " " << a[2] << " " << d << "\n";
      std::cout << HyperCube(0, m, a, d) << " " << HyperCubeAold(0, m, af, df) << std::endl;
    }
  }

  for(unsigned i = 0; i < 1000; i++) {
    a[0] = 0.;
    a[1] = 0.;
    a[2] = (0.5 * RAND_MAX - static_cast <myType>(rand())) / RAND_MAX;

    myType det = sqrt(a[0] * a[0] + a[1] * a[1] + a[2] * a[2]);
    a[0] /= det;
    a[1] /= det;
    a[2] /= det;

    d = 5. * (0.5 * RAND_MAX - static_cast <myType>(rand())) / RAND_MAX;

    af[0] = a[0];
    af[1] = a[1];
    af[2] = a[2];
    df = d;

    if(fabs(HyperCube(-1, m, a, d) - HyperCubeAold(-1, m, af, df)) > eps) {
      std::cout << "surface test failed" << std::endl;
      std::cout << a[0] << " " << a[1] << " " << a[2] << " " << d << "\n";
      std::cout << HyperCube(-1, m, a, d) << " " << HyperCubeAold(-1, m, af, df) << std::endl;
    }

    if(fabs(HyperCube(0, m, a, d) - HyperCubeAold(0, m, af, df)) > eps) {
      std::cout << "volume test failed" << std::endl;
      std::cout << a[0] << " " << a[1] << " " << a[2] << " " << d << "\n";
      std::cout << HyperCube(0, m, a, d) << " " << HyperCubeAold(0, m, af, df) << std::endl;
    }
  }

  for(unsigned i = 0; i < 1000; i++) {
    a[0] = 0.;
    a[1] = (0.5 * RAND_MAX - static_cast <myType>(rand())) / RAND_MAX;
    a[2] = 0.;

    myType det = sqrt(a[0] * a[0] + a[1] * a[1] + a[2] * a[2]);
    a[0] /= det;
    a[1] /= det;
    a[2] /= det;

    d = 5. * (0.5 * RAND_MAX - static_cast <myType>(rand())) / RAND_MAX;

    af[0] = a[0];
    af[1] = a[1];
    af[2] = a[2];
    df = d;

    if(fabs(HyperCube(-1, m, a, d) - HyperCubeAold(-1, m, af, df)) > eps) {
      std::cout << "surface test failed" << std::endl;
      std::cout << a[0] << " " << a[1] << " " << a[2] << " " << d << "\n";
      std::cout << HyperCube(-1, m, a, d) << " " << HyperCubeAold(-1, m, af, df) << std::endl;
    }

    if(fabs(HyperCube(0, m, a, d) - HyperCubeAold(0, m, af, df)) > eps) {
      std::cout << "volume test failed" << std::endl;
      std::cout << a[0] << " " << a[1] << " " << a[2] << " " << d << "\n";
      std::cout << HyperCube(0, m, a, d) << " " << HyperCubeAold(0, m, af, df) << std::endl;
    }
  }

  for(unsigned i = 0; i < 1000; i++) {
    a[0] = (0.5 * RAND_MAX - static_cast <myType>(rand())) / RAND_MAX;
    a[1] = 0.;
    a[2] = 0.;

    myType det = sqrt(a[0] * a[0] + a[1] * a[1] + a[2] * a[2]);
    a[0] /= det;
    a[1] /= det;
    a[2] /= det;

    d = 5. * (0.5 * RAND_MAX - static_cast <myType>(rand())) / RAND_MAX;

    af[0] = a[0];
    af[1] = a[1];
    af[2] = a[2];
    df = d;

    if(fabs(HyperCube(-1, m, a, d) - HyperCubeAold(-1, m, af, df)) > eps) {
      std::cout << "surface test failed" << std::endl;
      std::cout << a[0] << " " << a[1] << " " << a[2] << " " << d << "\n";
      std::cout << HyperCube(-1, m, a, d) << " " << HyperCubeAold(-1, m, af, df) << std::endl;
    }

    if(fabs(HyperCube(0, m, a, d) - HyperCubeAold(0, m, af, df)) > eps) {
      std::cout << "volume test failed" << std::endl;
      std::cout << a[0] << " " << a[1] << " " << a[2] << " " << d << "\n";
      std::cout << HyperCube(0, m, a, d) << " " << HyperCubeAold(0, m, af, df) << std::endl;
    }
  }

  for(unsigned i = 0; i < 1000; i++) {
    a[0] = 0.;
    a[1] = (0.5 * RAND_MAX - static_cast <myType>(rand())) / RAND_MAX;
    a[2] = (0.5 * RAND_MAX - static_cast <myType>(rand())) / RAND_MAX;

    myType det = sqrt(a[0] * a[0] + a[1] * a[1] + a[2] * a[2]);
    a[0] /= det;
    a[1] /= det;
    a[2] /= det;

    d = 5. * (0.5 * RAND_MAX - static_cast <myType>(rand())) / RAND_MAX;

    af[0] = a[0];
    af[1] = a[1];
    af[2] = a[2];
    df = d;

    if(fabs(HyperCube(-1, m, a, d) - HyperCubeAold(-1, m, af, df)) > eps) {
      std::cout << "surface test failed" << std::endl;
      std::cout << a[0] << " " << a[1] << " " << a[2] << " " << d << "\n";
      std::cout << HyperCube(-1, m, a, d) << " " << HyperCubeAold(-1, m, af, df) << std::endl;
    }

    if(fabs(HyperCube(0, m, a, d) - HyperCubeAold(0, m, af, df)) > eps) {
      std::cout << "volume test failed" << std::endl;
      std::cout << a[0] << " " << a[1] << " " << a[2] << " " << d << "\n";
      std::cout << HyperCube(0, m, a, d) << " " << HyperCubeAold(0, m, af, df) << std::endl;
    }
  }

  for(unsigned i = 0; i < 1000; i++) {
    a[0] = (0.5 * RAND_MAX - static_cast <myType>(rand())) / RAND_MAX;
    a[1] = 0.;
    a[2] = (0.5 * RAND_MAX - static_cast <myType>(rand())) / RAND_MAX;

    myType det = sqrt(a[0] * a[0] + a[1] * a[1] + a[2] * a[2]);
    a[0] /= det;
    a[1] /= det;
    a[2] /= det;

    d = 5. * (0.5 * RAND_MAX - static_cast <myType>(rand())) / RAND_MAX;

    af[0] = a[0];
    af[1] = a[1];
    af[2] = a[2];
    df = d;

    if(fabs(HyperCube(-1, m, a, d) - HyperCubeAold(-1, m, af, df)) > eps) {
      std::cout << "surface test failed" << std::endl;
      std::cout << a[0] << " " << a[1] << " " << a[2] << " " << d << "\n";
      std::cout << HyperCube(-1, m, a, d) << " " << HyperCubeAold(-1, m, af, df) << std::endl;
    }

    if(fabs(HyperCube(0, m, a, d) - HyperCubeAold(0, m, af, df)) > eps) {
      std::cout << "volume test failed" << std::endl;
      std::cout << a[0] << " " << a[1] << " " << a[2] << " " << d << "\n";
      std::cout << HyperCube(0, m, a, d) << " " << HyperCubeAold(0, m, af, df) << std::endl;
    }
  }

  for(unsigned i = 0; i < 1000; i++) {
    a[0] = (0.5 * RAND_MAX - static_cast <myType>(rand())) / RAND_MAX;
    a[1] = (0.5 * RAND_MAX - static_cast <myType>(rand())) / RAND_MAX;
    a[2] = 0.;

    myType det = sqrt(a[0] * a[0] + a[1] * a[1] + a[2] * a[2]);
    a[0] /= det;
    a[1] /= det;
    a[2] /= det;

    d = 5. * (0.5 * RAND_MAX - static_cast <myType>(rand())) / RAND_MAX;

    af[0] = a[0];
    af[1] = a[1];
    af[2] = a[2];
    df = d;

    if(fabs(HyperCube(-1, m, a, d) - HyperCubeAold(-1, m, af, df)) > eps) {
      std::cout << "surface test failed" << std::endl;
      std::cout << a[0] << " " << a[1] << " " << a[2] << " " << d << "\n";
      std::cout << HyperCube(-1, m, a, d) << " " << HyperCubeAold(-1, m, af, df) << std::endl;
    }

    if(fabs(HyperCube(0, m, a, d) - HyperCubeAold(0, m, af, df)) > eps) {
      std::cout << "volume test failed" << std::endl;
      std::cout << a[0] << " " << a[1] << " " << a[2] << " " << d << "\n";
      std::cout << HyperCube(0, m, a, d) << " " << HyperCubeAold(0, m, af, df) << std::endl;
    }
  }
}

template <class Float1>
void TestHexahedronTime(const Float1 &eps) {

  typedef typename boost::math::tools::promote_args<Float1>::type myType;

  std::cout << "testing the Hexahedron Time \n";
  std::vector<unsigned>m = {6, 6, 6};
  std::vector<myType>a(3);

  myType d;

  for(unsigned i = 0; i < 10000; i++) {
    a[0] = (0.5 * RAND_MAX - static_cast <myType>(rand())) / RAND_MAX;
    a[1] = (0.5 * RAND_MAX - static_cast <myType>(rand())) / RAND_MAX;
    a[2] = (0.5 * RAND_MAX - static_cast <myType>(rand())) / RAND_MAX;

    myType det = sqrt(a[0] * a[0] + a[1] * a[1] + a[2] * a[2]);
    a[0] /= det;
    a[1] /= det;
    a[2] /= det;

    d = 5. * (0.5 * RAND_MAX - static_cast <myType>(rand())) / RAND_MAX;

    HyperCube(-1, m, a, d);
    HyperCube(0, m, a, d);
    HyperCube(4, m, a, d);
  }
}
