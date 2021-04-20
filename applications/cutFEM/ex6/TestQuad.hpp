
template <class Float1>
void TestQuad(const Float1 &eps) {

  typedef typename boost::math::tools::promote_args<Float1>::type myType;
  typedef boost::multiprecision::cpp_bin_float_oct myTypeB;

  std::cout << "testing the Quadrilateral \n";
  std::vector<unsigned> m = {6, 7};
  std::vector<myType> a = {0.0000001, 0.0000001};
  std::vector<myTypeB> af = {0., 0.};

  myType d = 0.5;
  myTypeB df;

  for(unsigned i = 0; i < a.size();i++){
    af[0] = a[0];
    af[1] = a[1];
    //af[2] = a[2];
  }
  df = d;

  std::cout << HyperCubeB(6, a.size() - 1, m, a, d) << std::endl;
  std::cout << HyperCubeB(6, af.size() - 1, m, af, df) << std::endl;
  std::cout << HyperCube(6, m, a, d) << std::endl;
  std::cout << HyperCube(6, m, af, df) << std::endl;
  

  for(unsigned i = 0; i < 1000; i++) {
    a[0] = (0.5 * RAND_MAX - static_cast <myType>(rand())) / RAND_MAX;
    a[1] = (0.5 * RAND_MAX - static_cast <myType>(rand())) / RAND_MAX;

    myType det = sqrt(a[0] * a[0] + a[1] * a[1]);
    a[0] /= 1.e12 * det;
    a[1] /= 1. * det;

    d = 5 * (0.5 * RAND_MAX - static_cast <myType>(rand())) / RAND_MAX;

    af[0] = a[0];
    af[1] = a[1];
    df = d;

    if(fabs(HyperCube(-1, m, a, d) - HyperCubeB(-1, 1, m, af, df)) > eps) {
      std::cout << "surface test failed" << std::endl;
      std::cout << a[0] << " " << a[1] << " " << d << "\n";
      std::cout << HyperCube(-1, m, a, d) << " " << HyperCubeB(-1, 1, m, af, df) << std::endl;
    }

    if(fabs(HyperCube(0, m, a, d) - HyperCubeB(0, 1, m, af, df)) > eps) {
      std::cout << "volume test failed" << std::endl;
      std::cout << a[0] << " " << a[1] << " " << d << "\n";
      std::cout << HyperCube(0, m, a, d) << " " << HyperCubeB(0, 1, m, af, df) << std::endl;
    }
    
    if(fabs(HyperCube(5, m, a, d) - HyperCubeB(5, 1, m, af, df)) > eps) {
      std::cout << "s test failed" << std::endl;
      std::cout << a[0] << " " << a[1] << " " << d << "\n";
      std::cout << HyperCube(5, m, a, d) << " " << HyperCubeB(5, 1, m, af, df) <<  " " << HyperCube(5, m, af, df) << std::endl;
    }
  }

  for(unsigned i = 0; i < 1000; i++) {
    a[0] = 0.;
    a[1] = (0.5 * RAND_MAX - static_cast <myType>(rand())) / RAND_MAX;

    myType det = sqrt(a[0] * a[0] + a[1] * a[1]);
    a[0] /= det;
    a[1] /= det;

    d = 5. * (0.5 * RAND_MAX - static_cast <myType>(rand())) / RAND_MAX;

    af[0] = a[0];
    af[1] = a[1];
    df = d;

    if(fabs(HyperCube(-1, m, a, d) - HyperCubeB(-1, 1, m, af, df)) > eps) {
      std::cout << "surface test failed" << std::endl;
      std::cout << a[0] << " " << a[1] << " " << d << "\n";
      std::cout << HyperCube(-1, m, a, d) << " " << HyperCubeB(-1, 1, m, af, df) << std::endl;
    }

    if(fabs(HyperCube(0, m, a, d) - HyperCubeB(0, 1, m, af, df)) > eps) {
      std::cout << "volume test failed" << std::endl;
      std::cout << a[0] << " " << a[1] << " " << d << "\n";
      std::cout << HyperCube(0, m, a, d) << " " << HyperCubeB(0, 1, m, af, df) << std::endl;
    }
  }

  for(unsigned i = 0; i < 1000; i++) {
    a[0] = (0.5 * RAND_MAX - static_cast <myType>(rand())) / RAND_MAX;
    a[1] = 0.;

    myType det = sqrt(a[0] * a[0] + a[1] * a[1]);
    a[0] /= det;
    a[1] /= det;

    d = 5. * (0.5 * RAND_MAX - static_cast <myType>(rand())) / RAND_MAX;

    af[0] = a[0];
    af[1] = a[1];
    df = d;

    if(fabs(HyperCube(-1, m, a, d) - HyperCubeB(-1, 1, m, af, df)) > eps) {
      std::cout << "surface test failed" << std::endl;
      std::cout << a[0] << " " << a[1] << " " << d << "\n";
      std::cout << HyperCube(-1, m, a, d) << " " << HyperCubeB(-1, 1, m, af, df) << std::endl;
    }

    if(fabs(HyperCube(0, m, a, d) - HyperCubeB(0, 1, m, af, df)) > eps) {
      std::cout << "volume test failed" << std::endl;
      std::cout << a[0] << " " << a[1] << " " << d << "\n";
      std::cout << HyperCube(0, m, a, d) << " " << HyperCubeB(0, 1, m, af, df) << std::endl;
    }
  }
}

