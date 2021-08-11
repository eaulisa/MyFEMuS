template <class Float1>
void TestTriangle(const Float1 &eps) {

  typedef typename boost::math::tools::promote_args<Float1>::type myType;
  typedef boost::multiprecision::cpp_bin_float_oct myTypeB;

  std::cout << "testing the Triangle \n";
  std::vector<unsigned>m = {13, 5};
  std::vector<myType>a = { -1, 1};
  std::vector<myTypeB> af(2);
  af[0] = a[0];
  af[1] = a[1];


  myType d = 0.01;
  myTypeB df = d;

  int s = 6;

  std::cout.precision(14);

  std::cout << "s = " << s << " a = " << a[0] << " b = " << a[1] << " d = " << d <<std::endl;
  std::cout << TriangleFull(s, m, af, df) << std::endl;
  std::cout << TriangleA(s, m, a, d) << std::endl;

  
  
  
  myType c1 = 1.;
  myType c2 = 0.01;

  std::cout << "Epsilon cases\n";

  std::cout << "\na = 0, b = +-eps\n";

  a = {0.0, -c1};
  af[0] = a[0];
  af[1] = a[1];
  df = d = -c2;

  std::cout << "s = " << s << " a = " << a[0] << " b = " << a[1] << " d = " << d << "I = " << TriangleA(s,  m, a, d) << std::endl;
  if(fabs(TriangleA(s,  m, a, d) - TriangleFull(s,  m, af, df)) < eps)
    std::cout << "passed\n";
  else
    std::cout << "failed\n";

  a = { 0.0, -c1};
  af[0] = a[0];
  af[1] = a[1];
  df = d = c2;

  if(fabs(TriangleA(s,  m, a, d) - TriangleFull(s,  m, af, df)) < eps)
    std::cout << "passed\t";
  else
    std::cout << "failed\t";
  std::cout << "s = " << s << " a = " << a[0] << " b = " << a[1] << " d = " << d << " I = " << TriangleA(s,  m, a, d) << std::endl;


  a = { 0.0, c1};
  af[0] = a[0];
  af[1] = a[1];
  df = d = -c2;

  if(fabs(TriangleA(s,  m, a, d) - TriangleFull(s,  m, af, df)) < eps)
    std::cout << "passed\t";
  else
    std::cout << "failed\t";
  std::cout << "s = " << s << " a = " << a[0] << " b = " << a[1] << " d = " << d << " I = " << TriangleA(s,  m, a, d) << std::endl;


  a = { 0.0, c1};
  af[0] = a[0];
  af[1] = a[1];
  df = d = c2;

  if(fabs(TriangleA(s,  m, a, d) - TriangleFull(s,  m, af, df)) < eps)
    std::cout << "passed\t";
  else
    std::cout << "failed\t";
  std::cout << "s = " << s << " a = " << a[0] << " b = " << a[1] << " d = " << d << " I = " << TriangleA(s,  m, a, d) << std::endl;


  //////////////////////////////////////////////
  std::cout << "\na = +-eps, b = 0\n";

  a = {-c1, 0.0};
  af[0] = a[0];
  af[1] = a[1];
  df = d = -c2;
  if(fabs(TriangleA(s,  m, a, d) - TriangleFull(s,  m, af, df)) < eps)
    std::cout << "passed\t";
  else
    std::cout << "failed\t";
  std::cout << "s = " << s << " a = " << a[0] << " b = " << a[1] << " d = " << d << " I = " << TriangleA(s,  m, a, d) << std::endl;

  a = {-c1, 0.0};
  af[0] = a[0];
  af[1] = a[1];
  df = d = c2;
  if(fabs(TriangleA(s,  m, a, d) - TriangleFull(s,  m, af, df)) < eps)
    std::cout << "passed\t";
  else
    std::cout << "failed\t";
  std::cout << "s = " << s << " a = " << a[0] << " b = " << a[1] << " d = " << d << " I = " << TriangleA(s,  m, a, d) << std::endl;

  a = {c1, 0.0};
  af[0] = a[0];
  af[1] = a[1];
  df = d = -c2;
  if(fabs(TriangleA(s,  m, a, d) - TriangleFull(s,  m, af, df)) < eps)
    std::cout << "passed\t";
  else
    std::cout << "failed\t";
  std::cout << "s = " << s << " a = " << a[0] << " b = " << a[1] << " d = " << d << " I = " << TriangleA(s,  m, a, d) << std::endl;

  a = {c1, 0.0};
  af[0] = a[0];
  af[1] = a[1];
  df = d = c2;
  if(fabs(TriangleA(s,  m, a, d) - TriangleFull(s,  m, af, df)) < eps)
    std::cout << "passed\t";
  else
    std::cout << "failed\t";
  std::cout << "s = " << s << " a = " << a[0] << " b = " << a[1] << " d = " << d << " I = " << TriangleA(s,  m, a, d) << std::endl;

  //////////////////////////////////////////////
  std::cout << "\na = -+eps, b = +-eps\n";

  a = { c1, -c1};
  af[0] = a[0];
  af[1] = a[1];
  df = d = c2;
  if(fabs(TriangleA(s,  m, a, d) - TriangleFull(s,  m, af, df)) < eps)
    std::cout << "passed\t";
  else
    std::cout << "failed\t";
  std::cout << "s = " << s << " a = " << a[0] << " b = " << a[1] << " d = " << d << " I = " << TriangleA(s,  m, a, d) << std::endl;

  a = { c1, -c1};
  af[0] = a[0];
  af[1] = a[1];
  df = d = -c2;
  if(fabs(TriangleA(s,  m, a, d) - TriangleFull(s,  m, af, df)) < eps)
    std::cout << "passed\t";
  else
    std::cout << "failed\t";
  std::cout << "s = " << s << " a = " << a[0] << " b = " << a[1] << " d = " << d << " I = " << TriangleA(s,  m, a, d) << std::endl;

  a = {-c1, c1};
  af[0] = a[0];
  af[1] = a[1];
  df = d = c2;
  if(fabs(TriangleA(s,  m, a, d) - TriangleFull(s,  m, af, df)) < eps)
    std::cout << "passed\t";
  else
    std::cout << "failed\t";
  std::cout << "s = " << s << " a = " << a[0] << " b = " << a[1] << " d = " << d << " I = " << TriangleA(s,  m, a, d) << std::endl;


  a = {-c1, c1};
  af[0] = a[0];
  af[1] = a[1];
  df = d = -c2;
  if(fabs(TriangleA(s,  m, a, d) - TriangleFull(s,  m, af, df)) < eps)
    std::cout << "passed\t";
  else
    std::cout << "failed\t";
  std::cout << "s = " << s << " a = " << a[0] << " b = " << a[1] << " d = " << d << " I = " << TriangleA(s,  m, a, d) << std::endl;

  //////////////////////////////////////////////
  std::cout << "\na = +-eps, b = +-eps\n";

  a = { c1, c1};
  af[0] = a[0];
  af[1] = a[1];
  df = d = c2;
  if(fabs(TriangleA(s,  m, a, d) - TriangleFull(s,  m, af, df)) < eps)
    std::cout << "passed\t";
  else
    std::cout << "failed\t";
  std::cout << "s = " << s << " a = " << a[0] << " b = " << a[1] << " d = " << d << " I = " << TriangleA(s,  m, a, d) << std::endl;

  a = { c1, c1};
  af[0] = a[0];
  af[1] = a[1];
  df = d = -c2;
  if(fabs(TriangleA(s,  m, a, d) - TriangleFull(s,  m, af, df)) < eps)
    std::cout << "passed\t";
  else
    std::cout << "failed\t";
  std::cout << "s = " << s << " a = " << a[0] << " b = " << a[1] << " d = " << d << " I = " << TriangleA(s,  m, a, d) << std::endl;

  a = {-c1, -c1};
  af[0] = a[0];
  af[1] = a[1];
  df = d = c2;

  if(fabs(TriangleA(s,  m, a, d) - TriangleFull(s,  m, af, df)) < eps)
    std::cout << "passed\t";
  else
    std::cout << "failed\t";
  std::cout << "s = " << s << " a = " << a[0] << " b = " << a[1] << " d = " << d << " I = " << TriangleA(s,  m, a, d) << std::endl;

  a = {-c1, -c1};
  af[0] = a[0];
  af[1] = a[1];
  df = d = -c2;
  if(fabs(TriangleA(s,  m, a, d) - TriangleFull(s,  m, af, df)) < eps)
    std::cout << "passed\t";
  else
    std::cout << "failed\t";
  std::cout << "s = " << s << " a = " << a[0] << " b = " << a[1] << " d = " << d << " I = " << TriangleA(s,  m, a, d) << std::endl;

  //////////////////////////////////////////////
  std::cout << "\na = +-1, b = +-eps\n";

  a = { 1.0, -c1};
  af[0] = a[0];
  af[1] = a[1];
  df = d = -c2;
  if(fabs(TriangleA(s,  m, a, d) - TriangleFull(s,  m, af, df)) < eps)
    std::cout << "passed\t";
  else
    std::cout << "failed\t";
  std::cout << "s = " << s << " a = " << a[0] << " b = " << a[1] << " d = " << d << " I = " << TriangleA(s,  m, a, d) << std::endl;

  a = { 1.0, -c1};
  af[0] = a[0];
  af[1] = a[1];
  df = d = c2;
  if(fabs(TriangleA(s,  m, a, d) - TriangleFull(s,  m, af, df)) < eps)
    std::cout << "passed\t";
  else
    std::cout << "failed\t";
  std::cout << "s = " << s << " a = " << a[0] << " b = " << a[1] << " d = " << d << " I = " << TriangleA(s,  m, a, d) << std::endl;

  a = { 1.0, c1};
  af[0] = a[0];
  af[1] = a[1];
  df = d = -c2;
  if(fabs(TriangleA(s,  m, a, d) - TriangleFull(s,  m, af, df)) < eps)
    std::cout << "passed\t";
  else
    std::cout << "failed\t";
  std::cout << "s = " << s << " a = " << a[0] << " b = " << a[1] << " d = " << d << " I = " << TriangleA(s,  m, a, d) << std::endl;


  a = { 1.0, c1};
  af[0] = a[0];
  af[1] = a[1];
  df = d = c2;
  if(fabs(TriangleA(s,  m, a, d) - TriangleFull(s,  m, af, df)) < eps)
    std::cout << "passed\t";
  else
    std::cout << "failed\t";
  std::cout << "s = " << s << " a = " << a[0] << " b = " << a[1] << " d = " << d << " I = " << TriangleA(s,  m, a, d) << std::endl;


  a = { -1.0, -c1};
  af[0] = a[0];
  af[1] = a[1];
  df = d = -c2;
  if(fabs(TriangleA(s,  m, a, d) - TriangleFull(s,  m, af, df)) < eps)
    std::cout << "passed\t";
  else
    std::cout << "failed\t";
  std::cout << "s = " << s << " a = " << a[0] << " b = " << a[1] << " d = " << d << " I = " << TriangleA(s,  m, a, d) << std::endl;

  a = { -1.0, -c1};
  af[0] = a[0];
  af[1] = a[1];
  df = d = c2;
  if(fabs(TriangleA(s,  m, a, d) - TriangleFull(s,  m, af, df)) < eps)
    std::cout << "passed\t";
  else
    std::cout << "failed\t";
  std::cout << "s = " << s << " a = " << a[0] << " b = " << a[1] << " d = " << d << " I = " << TriangleA(s,  m, a, d) << std::endl;

  a = { -1.0, c1};
  af[0] = a[0];
  af[1] = a[1];
  df = d = -c2;
  if(fabs(TriangleA(s,  m, a, d) - TriangleFull(s,  m, af, df)) < eps)
    std::cout << "passed\t";
  else
    std::cout << "failed\t";
  std::cout << "s = " << s << " a = " << a[0] << " b = " << a[1] << " d = " << d << " I = " << TriangleA(s,  m, a, d) << std::endl;

  a = { -1.0, c1};
  af[0] = a[0];
  af[1] = a[1];
  df = d = c2;
  if(fabs(TriangleA(s,  m, a, d) - TriangleFull(s,  m, af, df)) < eps)
    std::cout << "passed\t";
  else
    std::cout << "failed\t";
  std::cout << "s = " << s << " a = " << a[0] << " b = " << a[1] << " d = " << d << " I = " << TriangleA(s,  m, a, d) << std::endl;


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
