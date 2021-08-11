template <class Float1>
void TestTriangle(const Float1 &eps) {

  typedef typename boost::math::tools::promote_args<Float1>::type myType;
  typedef boost::multiprecision::cpp_bin_float_oct myTypeB;

  myType eps1 = 1.0e-12;

  std::cout << "testing the Triangle \n";
  std::vector<unsigned>m = {13, 5};
  //std::vector<myType>a = { -1e-05, -1e-5};
  std::vector<myType>a = { 0., 0. };
  std::vector<myTypeB> af(2);
  af[0] = a[0];
  af[1] = a[1];

  myType d = 1;
  myTypeB df = d;

  int s = 6;

  std::cout.precision(14);

  std::cout << "a+b+d = "<<a[0] + a[1] +d <<std::endl;
  std::cout << "d = " << d << std::endl;
  std::cout << "d/b = "<<d/a[1] << std::endl;
  std::cout << "|a+b|/|a-b| = "  <<fabs((a[0]+a[1])/(a[0]-a[1]))<<std::endl;
  
  std::cout << "s + m + n= " << s+m[0]+m[1] << " a = " << a[0] << " b = " << a[1] << " d = " << d << std::endl;
  std::cout << TriangleA(s, m, af, df) << std::endl;
  std::cout << TriangleA(s, m, a, d) << std::endl;
  std::cout << fabs((TriangleA(s, m, af, df)- TriangleA(s, m, a, d))/TriangleA(s, m, af, df)) << std::endl; 
  
  std::cout << TriangleFull(s, m, af, df) << std::endl;
  std::cout << TriangleFull(s, m, a, d) << std::endl;
  std::cout << fabs((TriangleFull(s, m, af, df)- TriangleFull(s, m, a, d))/TriangleFull(s, m, af, df)) << std::endl; 
  
  
  //return;
  for(unsigned j = 0; j < 2; j++) {

    myType c1 = (j == 0) ? eps1 : 1.;
    myType c2 = (j == 0) ? 1. : eps1;

    std::cout << "Epsilon cases\n";

    std::cout << "\na = 0, b = +-eps\n";

    a = {0.0, -c1};
    af[0] = a[0];
    af[1] = a[1];
    df = d = -c2;


    myType I1 = TriangleA(s,  m, a, d);
    myTypeB I2 = TriangleA(s,  m, af, df);
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

    I1 = TriangleA(s,  m, a, d);
    I2 = TriangleA(s,  m, af, df);
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

    I1 = TriangleA(s,  m, a, d);
    I2 = TriangleA(s,  m, af, df);
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

    I1 = TriangleA(s,  m, a, d);
    I2 = TriangleA(s,  m, af, df);
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
    I1 = TriangleA(s,  m, a, d);
    I2 = TriangleA(s,  m, af, df);
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
    I1 = TriangleA(s,  m, a, d);
    I2 = TriangleA(s,  m, af, df);
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
    I1 = TriangleA(s,  m, a, d);
    I2 = TriangleA(s,  m, af, df);
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
    I1 = TriangleA(s,  m, a, d);
    I2 = TriangleA(s,  m, af, df);
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
    I1 = TriangleA(s,  m, a, d);
    I2 = TriangleA(s,  m, af, df);
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
    I1 = TriangleA(s,  m, a, d);
    I2 = TriangleA(s,  m, af, df);
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
    I1 = TriangleA(s,  m, a, d);
    I2 = TriangleA(s,  m, af, df);
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
    I1 = TriangleA(s,  m, a, d);
    I2 = TriangleA(s,  m, af, df);
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
    I1 = TriangleA(s,  m, a, d);
    I2 = TriangleA(s,  m, af, df);
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
    I1 = TriangleA(s,  m, a, d);
    I2 = TriangleA(s,  m, af, df);
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
    I1 = TriangleA(s,  m, a, d);
    I2 = TriangleA(s,  m, af, df);
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
    I1 = TriangleA(s,  m, a, d);
    I2 = TriangleA(s,  m, af, df);
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
    I1 = TriangleA(s,  m, a, d);
    I2 = TriangleA(s,  m, af, df);
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
    I1 = TriangleA(s,  m, a, d);
    I2 = TriangleA(s,  m, af, df);
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
    I1 = TriangleA(s,  m, a, d);
    I2 = TriangleA(s,  m, af, df);
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
    I1 = TriangleA(s,  m, a, d);
    I2 = TriangleA(s,  m, af, df);
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
    I1 = TriangleA(s,  m, a, d);
    I2 = TriangleA(s,  m, af, df);
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
    I1 = TriangleA(s,  m, a, d);
    I2 = TriangleA(s,  m, af, df);
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
    I1 = TriangleA(s,  m, a, d);
    I2 = TriangleA(s,  m, af, df);
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
    I1 = TriangleA(s,  m, a, d);
    I2 = TriangleA(s,  m, af, df);
    if((I2 != 0. &&  fabs((I1 - I2) / I2) < eps) || I1 == 0.) {
      std::cout << "passed ";
      std::cout << "s = " << s << " a = " << a[0] << " b = " << a[1] << " d = " << d << " I = " << I1 << std::endl;
    }
    else {
      std::cout << "Warning failed ";
      std::cout << "s = " << s << " a = " << a[0] << " b = " << a[1] << " d = " << d << " I1 =" << I1 << " I2 = " << I2 << std::endl;
    }
  }









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

    if(fabs(TriangleA(-1,  m, a, d) - TriangleFull(-1,  m, af, df)) > eps) {
      std::cout << "surface test failed" << std::endl;
      std::cout << a[0] << " " << a[1] << " " << " " << d << "\n";
      std::cout << TriangleFull(-1,  m, a, d) << " " << TriangleFull(-1,  m, af, df) << std::endl;
    }

    if(fabs(TriangleA(0,  m, a, d) - TriangleFull(0,  m, af, df)) > eps) {
      std::cout << "volume test failed" << std::endl;
      std::cout << a[0] << " " << a[1] << " " << d << "\n";
      std::cout << TriangleFull(0,  m, a, d) << " " << TriangleFull(0,  m, af, df) << std::endl;
    }
  }

  for(unsigned i = 0; i < 1000; i++) {

    a[0] = -0.05;
    a[1] = -a[0];
    d = 500. * (0.5 * RAND_MAX - static_cast <myType>(rand())) / RAND_MAX;

    if(a[0] + a[1] + d <= 0) {
      a[0] = -a[0];
      a[1] = -a[1];
      d = -d;
    }

    af[0] = a[0];
    af[1] = a[1];
    df = d;

    if(fabs(TriangleA(-1,  m, a, d) - TriangleFull(-1,  m, af, df)) > eps) {
      std::cout << "surface test failed" << std::endl;
      std::cout << a[0] << " " << a[1] << " " << " " << d << "\n";
      std::cout << TriangleA(-1,  m, a, d) << " " << TriangleFull(-1,  m, af, df) << std::endl;
    }

    if(fabs(TriangleA(0,  m, a, d) - TriangleFull(0,  m, af, df)) > eps) {
      std::cout << "volume test failed" << std::endl;
      std::cout << a[0] << " " << a[1] << " " << d << "\n";
      std::cout << TriangleA(0,  m, a, d) << " " << TriangleFull(0,  m, af, df) << std::endl;
    }
  }

  for(unsigned i = 0; i < 1000; i++) {

    a[0] = 0.;
    a[1] = 1.;
    d = 5. * (0.5 * RAND_MAX - static_cast <myType>(rand())) / RAND_MAX;

    if(a[0] + a[1] + d <= 0) {
      a[0] = -a[0];
      a[1] = -a[1];
      d = -d;
    }

    af[0] = a[0];
    af[1] = a[1];
    df = d;

    if(fabs(TriangleA(-1,  m, a, d) - TriangleFull(-1,  m, af, df)) > eps) {
      std::cout << "surface test failed" << std::endl;
      std::cout << a[0] << " " << a[1] << " " << " " << d << "\n";
      std::cout << TriangleA(-1,  m, a, d) << " " << TriangleFull(-1,  m, af, df) << std::endl;
    }

    if(fabs(TriangleA(0,  m, a, d) - TriangleFull(0,  m, af, df)) > eps) {
      std::cout << "volume test failed" << std::endl;
      std::cout << a[0] << " " << a[1] << " " << d << "\n";
      std::cout << TriangleA(0,  m, a, d) << " " << TriangleFull(0,  m, af, df) << std::endl;
    }
  }

  for(unsigned i = 0; i < 1000; i++) {

    a[0] = 1.;
    a[1] = 0.;
    d = 5. * (0.5 * RAND_MAX - static_cast <myType>(rand())) / RAND_MAX;

    if(a[0] + a[1] + d <= 0) {
      a[0] = -a[0];
      a[1] = -a[1];
      d = -d;
    }

    af[0] = a[0];
    af[1] = a[1];
    df = d;

    if(fabs(TriangleA(-1,  m, a, d) - TriangleFull(-1,  m, af, df)) > eps) {
      std::cout << "surface test failed" << std::endl;
      std::cout << a[0] << " " << a[1] << " " << " " << d << "\n";
      std::cout << TriangleA(-1,  m, a, d) << " " << TriangleFull(-1,  m, af, df) << std::endl;
    }

    if(fabs(TriangleA(0,  m, a, d) - TriangleFull(0,  m, af, df)) > eps) {
      std::cout << "volume test failed" << std::endl;
      std::cout << a[0] << " " << a[1] << " " << d << "\n";
      std::cout << TriangleA(0,  m, a, d) << " " << TriangleFull(0,  m, af, df) << std::endl;
    }
  }
}
