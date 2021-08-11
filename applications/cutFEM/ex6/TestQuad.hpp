
template <class Float1>
void TestQuad(const Float1 &eps) {

  typedef typename boost::math::tools::promote_args<Float1>::type myType;
  typedef boost::multiprecision::cpp_bin_float_oct myTypeB;

  std::cout << "testing the Quadrilateral \n";
  std::vector<unsigned> m = {6, 6};
  std::vector<myType> a = {0.0000000001, 0.0000000001};
  std::vector<myTypeB> af = {0., 0.};

  myType d = 0.5;
  myTypeB df;

  for(unsigned i = 0; i < a.size(); i++) {
    af[0] = a[0];
    af[1] = a[1];
    //af[2] = a[2];
  }
  df = d;

  unsigned s = 5;

  for(unsigned j = 0; j < 2; j++) {
    myType eps1 = 1.0e-5;
    myType c1 = (j == 0) ? eps1 : 1.;
    myType c2 = (j == 0) ? 1. : eps1;

    std::cout << "Epsilon cases\n";

    std::cout << "\na = 0, b = +-eps\n";

    a = {0.0, -c1};
    af[0] = a[0];
    af[1] = a[1];
    df = d = -c2;


    myType I1 = HyperCube(s,  m, a, d);
    myTypeB I2 = HyperCube(s,  m, af, df);
    if((I2 != 0. &&  fabs((I1 - I2) / I2) < eps) || I1 == 0.) {
      std::cout << "passed ";
      std::cout << "s = " << s << " a = " << a[0] << " b = " << a[1] << " d = " << d << " I = " << I1 << std::endl;
    }
    else {
      std::cout << "Warning failed ";
      std::cout << "s = " << s << " a = " << a[0] << " b = " << a[1] << " d = " << d << " I1 =" << I1 << " I2 = " << I2 << std::endl;
    }

    a = { 0.0, -c1};
    af[0] = a[0];
    af[1] = a[1];
    df = d = c2;

    I1 = HyperCube(s,  m, a, d);
    I2 = HyperCube(s,  m, af, df);
    if((I2 != 0. &&  fabs((I1 - I2) / I2) < eps) || I1 == 0.) {
      std::cout << "passed ";
      std::cout << "s = " << s << " a = " << a[0] << " b = " << a[1] << " d = " << d << " I = " << I1 << std::endl;
    }
    else {
      std::cout << "Warning failed ";
      std::cout << "s = " << s << " a = " << a[0] << " b = " << a[1] << " d = " << d << " I1 =" << I1 << " I2 = " << I2 << std::endl;
    }

    a = { 0.0, c1};
    af[0] = a[0];
    af[1] = a[1];
    df = d = -c2;

    I1 = HyperCube(s,  m, a, d);
    I2 = HyperCube(s,  m, af, df);
    if((I2 != 0. &&  fabs((I1 - I2) / I2) < eps) || I1 == 0.) {
      std::cout << "passed ";
      std::cout << "s = " << s << " a = " << a[0] << " b = " << a[1] << " d = " << d << " I = " << I1 << std::endl;
    }
    else {
      std::cout << "Warning failed ";
      std::cout << "s = " << s << " a = " << a[0] << " b = " << a[1] << " d = " << d << " I1 =" << I1 << " I2 = " << I2 << std::endl;
    }

    a = { 0.0, c1};
    af[0] = a[0];
    af[1] = a[1];
    df = d = c2;

    I1 = HyperCube(s,  m, a, d);
    I2 = HyperCube(s,  m, af, df);
    if((I2 != 0. &&  fabs((I1 - I2) / I2) < eps) || I1 == 0.) {
      std::cout << "passed ";
      std::cout << "s = " << s << " a = " << a[0] << " b = " << a[1] << " d = " << d << " I = " << I1 << std::endl;
    }
    else {
      std::cout << "Warning failed ";
      std::cout << "s = " << s << " a = " << a[0] << " b = " << a[1] << " d = " << d << " I1 =" << I1 << " I2 = " << I2 << std::endl;
    }


    //////////////////////////////////////////////
    std::cout << "\na = +-eps, b = 0\n";

    a = {-c1, 0.0};
    af[0] = a[0];
    af[1] = a[1];
    df = d = -c2;
    I1 = HyperCube(s,  m, a, d);
    I2 = HyperCube(s,  m, af, df);
    if((I2 != 0. &&  fabs((I1 - I2) / I2) < eps) || I1 == 0.) {
      std::cout << "passed ";
      std::cout << "s = " << s << " a = " << a[0] << " b = " << a[1] << " d = " << d << " I = " << I1 << std::endl;
    }
    else {
      std::cout << "Warning failed ";
      std::cout << "s = " << s << " a = " << a[0] << " b = " << a[1] << " d = " << d << " I1 =" << I1 << " I2 = " << I2 << std::endl;
    }

    a = {-c1, 0.0};
    af[0] = a[0];
    af[1] = a[1];
    df = d = c2;
    I1 = HyperCube(s,  m, a, d);
    I2 = HyperCube(s,  m, af, df);
    if((I2 != 0. &&  fabs((I1 - I2) / I2) < eps) || I1 == 0.) {
      std::cout << "passed ";
      std::cout << "s = " << s << " a = " << a[0] << " b = " << a[1] << " d = " << d << " I = " << I1 << std::endl;
    }
    else {
      std::cout << "Warning failed ";
      std::cout << "s = " << s << " a = " << a[0] << " b = " << a[1] << " d = " << d << " I1 =" << I1 << " I2 = " << I2 << std::endl;
    }

    a = {c1, 0.0};
    af[0] = a[0];
    af[1] = a[1];
    df = d = -c2;
    I1 = HyperCube(s,  m, a, d);
    I2 = HyperCube(s,  m, af, df);
    if((I2 != 0. &&  fabs((I1 - I2) / I2) < eps) || I1 == 0.) {
      std::cout << "passed ";
      std::cout << "s = " << s << " a = " << a[0] << " b = " << a[1] << " d = " << d << " I = " << I1 << std::endl;
    }
    else {
      std::cout << "Warning failed ";
      std::cout << "s = " << s << " a = " << a[0] << " b = " << a[1] << " d = " << d << " I1 =" << I1 << " I2 = " << I2 << std::endl;
    }

    a = {c1, 0.0};
    af[0] = a[0];
    af[1] = a[1];
    df = d = c2;
    I1 = HyperCube(s,  m, a, d);
    I2 = HyperCube(s,  m, af, df);
    if((I2 != 0. &&  fabs((I1 - I2) / I2) < eps) || I1 == 0.) {
      std::cout << "passed ";
      std::cout << "s = " << s << " a = " << a[0] << " b = " << a[1] << " d = " << d << " I = " << I1 << std::endl;
    }
    else {
      std::cout << "Warning failed ";
      std::cout << "s = " << s << " a = " << a[0] << " b = " << a[1] << " d = " << d << " I1 =" << I1 << " I2 = " << I2 << std::endl;
    }


    //////////////////////////////////////////////
    std::cout << "\na = -+eps, b = +-eps\n";

    a = { c1, -c1};
    af[0] = a[0];
    af[1] = a[1];
    df = d = c2;
    I1 = HyperCube(s,  m, a, d);
    I2 = HyperCube(s,  m, af, df);
    if((I2 != 0. &&  fabs((I1 - I2) / I2) < eps) || I1 == 0.) {
      std::cout << "passed ";
      std::cout << "s = " << s << " a = " << a[0] << " b = " << a[1] << " d = " << d << " I = " << I1 << std::endl;
    }
    else {
      std::cout << "Warning failed ";
      std::cout << "s = " << s << " a = " << a[0] << " b = " << a[1] << " d = " << d << " I1 =" << I1 << " I2 = " << I2 << std::endl;
    }

    a = { c1, -c1};
    af[0] = a[0];
    af[1] = a[1];
    df = d = -c2;
    I1 = HyperCube(s,  m, a, d);
    I2 = HyperCube(s,  m, af, df);
    if((I2 != 0. &&  fabs((I1 - I2) / I2) < eps) || I1 == 0.) {
      std::cout << "passed ";
      std::cout << "s = " << s << " a = " << a[0] << " b = " << a[1] << " d = " << d << " I = " << I1 << std::endl;
    }
    else {
      std::cout << "Warning failed ";
      std::cout << "s = " << s << " a = " << a[0] << " b = " << a[1] << " d = " << d << " I1 =" << I1 << " I2 = " << I2 << std::endl;
    }

    a = {-c1, c1};
    af[0] = a[0];
    af[1] = a[1];
    df = d = c2;
    I1 = HyperCube(s,  m, a, d);
    I2 = HyperCube(s,  m, af, df);
    if((I2 != 0. &&  fabs((I1 - I2) / I2) < eps) || I1 == 0.) {
      std::cout << "passed ";
      std::cout << "s = " << s << " a = " << a[0] << " b = " << a[1] << " d = " << d << " I = " << I1 << std::endl;
    }
    else {
      std::cout << "Warning failed ";
      std::cout << "s = " << s << " a = " << a[0] << " b = " << a[1] << " d = " << d << " I1 =" << I1 << " I2 = " << I2 << std::endl;
    }


    a = {-c1, c1};
    af[0] = a[0];
    af[1] = a[1];
    df = d = -c2;
    I1 = HyperCube(s,  m, a, d);
    I2 = HyperCube(s,  m, af, df);
    if((I2 != 0. &&  fabs((I1 - I2) / I2) < eps) || I1 == 0.) {
      std::cout << "passed ";
      std::cout << "s = " << s << " a = " << a[0] << " b = " << a[1] << " d = " << d << " I = " << I1 << std::endl;
    }
    else {
      std::cout << "Warning failed ";
      std::cout << "s = " << s << " a = " << a[0] << " b = " << a[1] << " d = " << d << " I1 =" << I1 << " I2 = " << I2 << std::endl;
    }

    //////////////////////////////////////////////
    std::cout << "\na = +-eps, b = +-eps\n";

    a = { c1, c1};
    af[0] = a[0];
    af[1] = a[1];
    df = d = c2;
    I1 = HyperCube(s,  m, a, d);
    I2 = HyperCube(s,  m, af, df);
    if((I2 != 0. &&  fabs((I1 - I2) / I2) < eps) || I1 == 0.) {
      std::cout << "passed ";
      std::cout << "s = " << s << " a = " << a[0] << " b = " << a[1] << " d = " << d << " I = " << I1 << std::endl;
    }
    else {
      std::cout << "Warning failed ";
      std::cout << "s = " << s << " a = " << a[0] << " b = " << a[1] << " d = " << d << " I1 =" << I1 << " I2 = " << I2 << std::endl;
    }

    a = { c1, c1};
    af[0] = a[0];
    af[1] = a[1];
    df = d = -c2;
    I1 = HyperCube(s,  m, a, d);
    I2 = HyperCube(s,  m, af, df);
    if((I2 != 0. &&  fabs((I1 - I2) / I2) < eps) || I1 == 0.) {
      std::cout << "passed ";
      std::cout << "s = " << s << " a = " << a[0] << " b = " << a[1] << " d = " << d << " I = " << I1 << std::endl;
    }
    else {
      std::cout << "Warning failed ";
      std::cout << "s = " << s << " a = " << a[0] << " b = " << a[1] << " d = " << d << " I1 =" << I1 << " I2 = " << I2 << std::endl;
    }

    a = {-c1, -c1};
    af[0] = a[0];
    af[1] = a[1];
    df = d = c2;
    I1 = HyperCube(s,  m, a, d);
    I2 = HyperCube(s,  m, af, df);
    if((I2 != 0. &&  fabs((I1 - I2) / I2) < eps) || I1 == 0.) {
      std::cout << "passed ";
      std::cout << "s = " << s << " a = " << a[0] << " b = " << a[1] << " d = " << d << " I = " << I1 << std::endl;
    }
    else {
      std::cout << "Warning failed ";
      std::cout << "s = " << s << " a = " << a[0] << " b = " << a[1] << " d = " << d << " I1 =" << I1 << " I2 = " << I2 << std::endl;
    }

    a = {-c1, -c1};
    af[0] = a[0];
    af[1] = a[1];
    df = d = -c2;
    I1 = HyperCube(s,  m, a, d);
    I2 = HyperCube(s,  m, af, df);
    if((I2 != 0. &&  fabs((I1 - I2) / I2) < eps) || I1 == 0.) {
      std::cout << "passed ";
      std::cout << "s = " << s << " a = " << a[0] << " b = " << a[1] << " d = " << d << " I = " << I1 << std::endl;
    }
    else {
      std::cout << "Warning failed ";
      std::cout << "s = " << s << " a = " << a[0] << " b = " << a[1] << " d = " << d << " I1 =" << I1 << " I2 = " << I2 << std::endl;
    }

    //////////////////////////////////////////////
    std::cout << "\na = +-1, b = +-eps\n";

    a = { 1.0, -c1};
    af[0] = a[0];
    af[1] = a[1];
    df = d = -c2;
    I1 = HyperCube(s,  m, a, d);
    I2 = HyperCube(s,  m, af, df);
    if((I2 != 0. &&  fabs((I1 - I2) / I2) < eps) || I1 == 0.) {
      std::cout << "passed ";
      std::cout << "s = " << s << " a = " << a[0] << " b = " << a[1] << " d = " << d << " I = " << I1 << std::endl;
    }
    else {
      std::cout << "Warning failed ";
      std::cout << "s = " << s << " a = " << a[0] << " b = " << a[1] << " d = " << d << " I1 =" << I1 << " I2 = " << I2 << std::endl;
    }

    a = { 1.0, -c1};
    af[0] = a[0];
    af[1] = a[1];
    df = d = c2;
    I1 = HyperCube(s,  m, a, d);
    I2 = HyperCube(s,  m, af, df);
    if((I2 != 0. &&  fabs((I1 - I2) / I2) < eps) || I1 == 0.) {
      std::cout << "passed ";
      std::cout << "s = " << s << " a = " << a[0] << " b = " << a[1] << " d = " << d << " I = " << I1 << std::endl;
    }
    else {
      std::cout << "Warning failed ";
      std::cout << "s = " << s << " a = " << a[0] << " b = " << a[1] << " d = " << d << " I1 =" << I1 << " I2 = " << I2 << std::endl;
    }

    a = { 1.0, c1};
    af[0] = a[0];
    af[1] = a[1];
    df = d = -c2;
    I1 = HyperCube(s,  m, a, d);
    I2 = HyperCube(s,  m, af, df);
    if((I2 != 0. &&  fabs((I1 - I2) / I2) < eps) || I1 == 0.) {
      std::cout << "passed ";
      std::cout << "s = " << s << " a = " << a[0] << " b = " << a[1] << " d = " << d << " I = " << I1 << std::endl;
    }
    else {
      std::cout << "Warning failed ";
      std::cout << "s = " << s << " a = " << a[0] << " b = " << a[1] << " d = " << d << " I1 =" << I1 << " I2 = " << I2 << std::endl;
    }


    a = { 1.0, c1};
    af[0] = a[0];
    af[1] = a[1];
    df = d = c2;
    I1 = HyperCube(s,  m, a, d);
    I2 = HyperCube(s,  m, af, df);
    if((I2 != 0. &&  fabs((I1 - I2) / I2) < eps) || I1 == 0.) {
      std::cout << "passed ";
      std::cout << "s = " << s << " a = " << a[0] << " b = " << a[1] << " d = " << d << " I = " << I1 << std::endl;
    }
    else {
      std::cout << "Warning failed ";
      std::cout << "s = " << s << " a = " << a[0] << " b = " << a[1] << " d = " << d << " I1 =" << I1 << " I2 = " << I2 << std::endl;
    }


    a = { -1.0, -c1};
    af[0] = a[0];
    af[1] = a[1];
    df = d = -c2;
    I1 = HyperCube(s,  m, a, d);
    I2 = HyperCube(s,  m, af, df);
    if((I2 != 0. &&  fabs((I1 - I2) / I2) < eps) || I1 == 0.) {
      std::cout << "passed ";
      std::cout << "s = " << s << " a = " << a[0] << " b = " << a[1] << " d = " << d << " I = " << I1 << std::endl;
    }
    else {
      std::cout << "Warning failed ";
      std::cout << "s = " << s << " a = " << a[0] << " b = " << a[1] << " d = " << d << " I1 =" << I1 << " I2 = " << I2 << std::endl;
    }

    a = { -1.0, -c1};
    af[0] = a[0];
    af[1] = a[1];
    df = d = c2;
    I1 = HyperCube(s,  m, a, d);
    I2 = HyperCube(s,  m, af, df);
    if((I2 != 0. &&  fabs((I1 - I2) / I2) < eps) || I1 == 0.) {
      std::cout << "passed ";
      std::cout << "s = " << s << " a = " << a[0] << " b = " << a[1] << " d = " << d << " I = " << I1 << std::endl;
    }
    else {
      std::cout << "Warning failed ";
      std::cout << "s = " << s << " a = " << a[0] << " b = " << a[1] << " d = " << d << " I1 =" << I1 << " I2 = " << I2 << std::endl;
    }

    a = { -1.0, c1};
    af[0] = a[0];
    af[1] = a[1];
    df = d = -c2;
    I1 = HyperCube(s,  m, a, d);
    I2 = HyperCube(s,  m, af, df);
    if((I2 != 0. &&  fabs((I1 - I2) / I2) < eps) || I1 == 0.) {
      std::cout << "passed ";
      std::cout << "s = " << s << " a = " << a[0] << " b = " << a[1] << " d = " << d << " I = " << I1 << std::endl;
    }
    else {
      std::cout << "Warning failed ";
      std::cout << "s = " << s << " a = " << a[0] << " b = " << a[1] << " d = " << d << " I1 =" << I1 << " I2 = " << I2 << std::endl;
    }

    a = { -1.0, c1};
    af[0] = a[0];
    af[1] = a[1];
    df = d = c2;
    I1 = HyperCube(s,  m, a, d);
    I2 = HyperCube(s, m, af, df);
    if((I2 != 0. &&  fabs((I1 - I2) / I2) < eps) || I1 == 0.) {
      std::cout << "passed ";
      std::cout << "s = " << s << " a = " << a[0] << " b = " << a[1] << " d = " << d << " I = " << I1 << std::endl;
    }
    else {
      std::cout << "Warning failed ";
      std::cout << "s = " << s << " a = " << a[0] << " b = " << a[1] << " d = " << d << " I1 =" << I1 << " I2 = " << I2 << std::endl;
    }
  }




  std::cout << HyperCubeBOld(6, a.size() - 1, m, a, d) << std::endl;
  std::cout << HyperCubeBOld(6, af.size() - 1, m, af, df) << std::endl;
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

    if(fabs(HyperCube(-1, m, a, d) - HyperCubeBOld(-1, 1, m, af, df)) > eps) {
      std::cout << "surface test failed" << std::endl;
      std::cout << a[0] << " " << a[1] << " " << d << "\n";
      std::cout << HyperCube(-1, m, a, d) << " " << HyperCubeBOld(-1, 1, m, af, df) << std::endl;
    }

    if(fabs(HyperCube(0, m, a, d) - HyperCubeBOld(0, 1, m, af, df)) > eps) {
      std::cout << "volume test failed" << std::endl;
      std::cout << a[0] << " " << a[1] << " " << d << "\n";
      std::cout << HyperCube(0, m, a, d) << " " << HyperCubeBOld(0, 1, m, af, df) << std::endl;
    }

    if(fabs(HyperCube(5, m, a, d) - HyperCubeBOld(5, 1, m, af, df)) > eps) {
      std::cout << "s test failed" << std::endl;
      std::cout << a[0] << " " << a[1] << " " << d << "\n";
      std::cout << HyperCube(5, m, a, d) << " " << HyperCubeBOld(5, 1, m, af, df) <<  " " << HyperCube(5, m, af, df) << std::endl;
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

    if(fabs(HyperCube(-1, m, a, d) - HyperCubeBOld(-1, 1, m, af, df)) > eps) {
      std::cout << "surface test failed" << std::endl;
      std::cout << a[0] << " " << a[1] << " " << d << "\n";
      std::cout << HyperCube(-1, m, a, d) << " " << HyperCubeBOld(-1, 1, m, af, df) << std::endl;
    }

    if(fabs(HyperCube(0, m, a, d) - HyperCubeBOld(0, 1, m, af, df)) > eps) {
      std::cout << "volume test failed" << std::endl;
      std::cout << a[0] << " " << a[1] << " " << d << "\n";
      std::cout << HyperCube(0, m, a, d) << " " << HyperCubeBOld(0, 1, m, af, df) << std::endl;
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

    if(fabs(HyperCube(-1, m, a, d) - HyperCubeBOld(-1, 1, m, af, df)) > eps) {
      std::cout << "surface test failed" << std::endl;
      std::cout << a[0] << " " << a[1] << " " << d << "\n";
      std::cout << HyperCube(-1, m, a, d) << " " << HyperCubeBOld(-1, 1, m, af, df) << std::endl;
    }

    if(fabs(HyperCube(0, m, a, d) - HyperCubeBOld(0, 1, m, af, df)) > eps) {
      std::cout << "volume test failed" << std::endl;
      std::cout << a[0] << " " << a[1] << " " << d << "\n";
      std::cout << HyperCube(0, m, a, d) << " " << HyperCubeBOld(0, 1, m, af, df) << std::endl;
    }
  }
}

