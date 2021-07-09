template <class Float1>
void TestLine( const Float1 &eps) {
    
  typedef typename boost::math::tools::promote_args<Float1>::type myType;    
  typedef boost::multiprecision::cpp_bin_float_oct myTypeB;
  
  myType a;
  myType d;

  myTypeB af;
  myTypeB df;

  unsigned m = 5;
  int ss = 4;

  std::cout << "testing the Line Integrals \n";

  for(unsigned i = 0; i < 1000; i++) {
    a = 0.1 * (0.5 * RAND_MAX - static_cast <myType>(rand())) / RAND_MAX;;
    d = 10. * (0.5 * RAND_MAX - static_cast <myType>(rand())) / RAND_MAX;

    af = a;
    df = d;

    if(fabs(Intm1to1LimLiC(-1, m, af, df) - Intm1to1LimLi(-1, m, a, d)) > eps) {
      std::cout << "surface test failed" << std::endl;
      std::cout << "normal = " << a << " " << d << "\n";
      std::cout << Intm1to1LimLiC(-1, m, af, df) << " " << Intm1to1LimLi(-1, m, a, d) << std::endl;
    }

    if(fabs(Intm1to1LimLiC(0, m, af, df) - Intm1to1LimLi(0, m, a, d)) > eps) {
      std::cout << "volume test failed" << std::endl;
      std::cout << "normal = " << a << " " << d << "\n";
      std::cout << Intm1to1LimLiC(0, m, af, df) << " " << Intm1to1LimLi(0, m, a, d) << std::endl;
    }

    if(fabs(Intm1to1LimLiC(ss, m, a, d) - Intm1to1LimLi(ss, m, a, d)) > eps) {
      std::cout << "ss test failed" << std::endl;
      std::cout << "normal = " << a << " " << d << "\n";
      std::cout << Intm1to1LimLiC(ss, m, a, d) << " " << Intm1to1LimLi(ss, m, a, d) << std::endl;
    }

    //////////////////////////////////////////////////////////////////

    if(fabs(Int0to1LimLiC(-1, m, af, df) - Int0to1LimLi(-1, m, a, d)) > eps) {
      std::cout << "surface test failed" << std::endl;
      std::cout << "normal = " << a << " " << d << "\n";
      std::cout << Int0to1LimLiC(-1, m, af, df) << " " << Int0to1LimLi(-1, m, a, d) << std::endl;
    }

    if(fabs(Int0to1LimLiC(0, m, af, df) - Int0to1LimLi(0, m, a, d)) > eps) {
      std::cout << "volume test failed" << std::endl;
      std::cout << "normal = " << a << " " << d << "\n";
      std::cout << Int0to1LimLiC(0, m, af, df) << " " << Int0to1LimLi(0, m, a, d) << std::endl;
    }

    if(fabs(Int0to1LimLiC(ss, m, a, d) - Int0to1LimLi(ss, m, a, d)) > eps) {
      std::cout << "ss test failed" << std::endl;
      std::cout << "normal = " << a << " " << d << "\n";
      std::cout << Int0to1LimLiC(ss, m, a, d) << " " << Int0to1LimLi(ss, m, a, d) << std::endl;
    }
  }

}
