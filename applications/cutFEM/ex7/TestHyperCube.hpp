
template <class Float1>
void TestQuad(Float1 &eps) {

  typedef typename boost::math::tools::promote_args<Float1>::type myType;
  
  //typedef double myTypeB;
  //typedef boost::multiprecision::cpp_bin_float_quad myType;
  typedef boost::multiprecision::cpp_bin_float_oct myTypeB;
  
  

  std::cout << "testing the Quadrilateral \n";
  std::vector<unsigned> m = {13, 5};
  std::vector<myType> a = {0.0000000001, 0.0000000001};
  std::vector<myTypeB> af = {0., 0.};

  myType d = 0.5;
  myTypeB df;

  std::cout.precision(14);

  for(unsigned i = 0; i < a.size(); i++) {
    af[0] = a[0];
    af[1] = a[1];
    //af[2] = a[2];
  }
  df = d;

  unsigned s = 5;

  myType eps1 = 1.0e-10;
  std::vector<std::vector<myType>> epsCut{{0, 0, 0},
    {eps1, 0, 0}, {-eps1, 0, 0}, {0, eps1, 0}, {0, -eps1, 0}, {0, 0, eps1}, {0, 0, -eps1},
    {eps1, eps1, 0}, {eps1, -eps1, 0}, {-eps1, eps1, 0}, {-eps1, -eps1, 0},
    {eps1, -eps1, eps1}, {-eps1, -eps1, eps1}, {eps1, -eps1, -eps1}, {-eps1, -eps1, -eps1}
  };

//   std::vector<std::vector<myType>> smallCut{{0, 0, 0}, {0, 0, 1}, {0, 0, -1},
//     {-1, -1, 0}, {-1, 0, 0}, {-1, 1, -1}, {0, 1, -1},
//     {1, 1, -2}, {1, 0, -1}, {1, -1, -1}, {0, -1, 0},
//     {1, 1, 0}, {1, 0, 0}, {1, -1, 1}, {0, -1, 1},
//     {-1, -1, 2}, {-1, 0, 1}, {-1, 1, 1}, {0, 1, 0}
// };

  std::vector<std::vector<myType>> smallCut{{0, 0, 0}, {0, 0, 1}, {0, 0, -1},
    {-1, -1, -2}, {-1, 0, -1}, {-1, 1, -2}, {0, 1, -1},
    {1, 1, -2}, {1, 0, -1}, {1, -1, -2}, {0, -1, -1},
    {1, 1, 2}, {1, 0, 1}, {1, -1, 2}, {0, -1, 1},
    {-1, -1, 2}, {-1, 0, 1}, {-1, 1, 2}, {0, 1, 1}
  };


  for(unsigned i = 0; i < smallCut.size(); i++) {
    for(unsigned j = 0; j < epsCut.size(); j++) {

      a[0] = smallCut[i][0] + epsCut[j][0];
      a[1] = smallCut[i][1] + epsCut[j][1];
      d = smallCut[i][2] + epsCut[j][2];

      for(unsigned i = 0; i < a.size(); i++) {
        d -= a[i];
        a[i] *= 2;
      }

      af[0] = static_cast<myType>(a[0]);
      af[1] = static_cast<myType>(a[1]);
      df = static_cast<myType>(d);

      myType I1 = HyperCube(s,  m, a, d);
      myTypeB I2 = HyperCube(s,  m, af, df);
      if((I2 != 0. &&  fabs((I1 - I2) / I2) < eps) || I1 == 0.) {
        //std::cout << "passed ";
        //std::cout << "s = " << s << " a = " << a[0] << " b = " << a[1] << " d = " << d << " I = " << I1 << std::endl;
      }
      else {
        std::cout << "Warning failed ";
        std::cout << "s = " << s << " a = " << a[0] << " b = " << a[1] << " d = " << d << " I1 =" << I1 << " I2 = " << I2 << std::endl;
      }
    }
    //std::cout << std::endl;
  }

  eps1 = 1.0e-10;
  epsCut = { {0, 0, 0, 0},
    {eps1, 0, 0, 0}, {-eps1, 0, 0, 0},
    {0, eps1, 0, 0}, {0, -eps1, 0, 0},
    {0, 0, eps1, 0}, {0, 0, -eps1, 0},
    {0, 0, 0, eps1}, {0, 0, 0, -eps1},
    {eps1, eps1, eps1, 0}, {eps1, -eps1, eps1, 0},
    {eps1, eps1, -eps1, 0},{eps1, -eps1, -eps1, 0},
    {-eps1, eps1, eps1, 0}, {-eps1, -eps1, eps1, 0},
    {-eps1, eps1, -eps1, 0}, {-eps1, -eps1, -eps1, eps1},
    {eps1, eps1, eps1, eps1}, {eps1, -eps1, eps1, eps1},
    {eps1, eps1, -eps1, eps1}, {eps1, -eps1, -eps1, eps1},
    {-eps1, eps1, eps1, eps1}, {-eps1, -eps1, eps1, eps1},
    {-eps1, eps1, -eps1, eps1}, {-eps1, -eps1, -eps1, eps1},
    {eps1, eps1, eps1, -eps1}, {eps1, -eps1, eps1, -eps1},
    {eps1, eps1, -eps1, -eps1}, {eps1, -eps1, -eps1, -eps1},
    {-eps1, eps1, eps1, -eps1}, {-eps1, -eps1, eps1, -eps1},
    {-eps1, eps1, -eps1, -eps1}, {-eps1, -eps1, -eps1, -eps1},
  };

  smallCut = {
    {0, 0, 0, 0}, {0, 0, 0, 1}, {0, 0, 0, -1},
    {-1, -1, -1, -3}, {-1, 1, -1, -3}, {1, 1, -1, -3}, {1, -1, -1, -3},
    {-1, -1, 1, -3}, {-1, 1, 1, -3}, {1, 1, 1, -3}, {1, -1, 1, -3},
    {1, 1, 1, 3}, {1, -1, 1, 3}, {-1, -1, 1, 3}, {-1, 1, 1, 3},
    {1, 1, -1, 3}, {1, -1, -1, 3}, {-1, -1, -1, 3}, {-1, 1, -1, 3},
    {-1, 0, 0, -1}, {1, 0, 0, - 1}, {0, -1, 0, -1}, {0, 1, 0, -1}, {0, 0, -1, -1}, {0, 0, 1, -1},
    {1, 0, 0, 1}, {-1, 0, 0, 1}, {0, 1, 0, 1}, {0, -1, 0, 1}, {0, 0, 1, 1}, {0, 0, -1, 1},
  };

  a.resize(3);  
  af.resize(3);  
  m = {5,6,7};
  s = 0;
  eps =1.0e-10;
//   for(unsigned i = 0; i < smallCut.size(); i++) {
//     for(unsigned j = 0; j < epsCut.size(); j++) {
//      
//       a[0] = smallCut[i][0] + epsCut[j][0];
//       a[1] = smallCut[i][1] + epsCut[j][1];
//       a[2] = smallCut[i][2] + epsCut[j][2];
//       d = smallCut[i][3] + epsCut[j][3];
// 
//       for(unsigned i = 0; i < a.size(); i++) {
//         d -= a[i];
//         a[i] *= 2;
//       }
// 
//       af[0] = static_cast<myType>(a[0]);
//       af[1] = static_cast<myType>(a[1]);
//       af[2] = static_cast<myType>(a[2]);
//       df = static_cast<myType>(d);
// 
//       myType I1 = HyperCube(s,  m, a, d);
//       myTypeB I2 = HyperCube(s,  m, af, df);
//       if((I2 != 0. &&  fabs((I1 - I2) / I2) < eps) || I1 == 0.) {
//         //std::cout << "passed ";
//         //std::cout << "s = " << s << " a = " << a[0] << " b = " << a[1] << " c = " << a[2] << " d = " << d << " I = " << I1 << std::endl;
//       }
//       else {
//         std::cout << "Warning failed ";
//         std::cout << "s = " << s << " a = " << a[0] << " b = " << a[1] << " c = " << a[2] << " d = " << d << " I1 =" << I1 << " I2 = " << I2 << std::endl;
//       }
//     }
//     //std::cout << std::endl;
//   }




}

template <class myType>
void TestHex(myType &eps) {

  //typedef typename boost::math::tools::promote_args<Float1>::type myType;
  
  //typedef double myTypeB;
  //typedef boost::multiprecision::cpp_bin_float_quad myType;
  typedef boost::multiprecision::cpp_bin_float_oct myTypeB;
  
  

  std::cout << "testing the Hexahedron \n";
  std::vector<unsigned> m = {5, 6, 7};
  std::vector<myType> a;
  std::vector<myTypeB> af;
  
  a.reserve(100);
  af.reserve(100);

  myType d = 0.5;
  myTypeB df = 0.5;

  std::cout.precision(14);


  double eps1 = 1.0e-10;
  std::vector<std::vector<double>> epsCut = { {0, 0, 0, 0},
    {eps1, 0, 0, 0}, {-eps1, 0, 0, 0},
    {0, eps1, 0, 0}, {0, -eps1, 0, 0},
    {0, 0, eps1, 0}, {0, 0, -eps1, 0},
    {0, 0, 0, eps1}, {0, 0, 0, -eps1},
    {eps1, eps1, eps1, 0}, {eps1, -eps1, eps1, 0},
    {eps1, eps1, -eps1, 0},{eps1, -eps1, -eps1, 0},
    {-eps1, eps1, eps1, 0}, {-eps1, -eps1, eps1, 0},
    {-eps1, eps1, -eps1, 0}, {-eps1, -eps1, -eps1, eps1},
    {eps1, eps1, eps1, eps1}, {eps1, -eps1, eps1, eps1},
    {eps1, eps1, -eps1, eps1}, {eps1, -eps1, -eps1, eps1},
    {-eps1, eps1, eps1, eps1}, {-eps1, -eps1, eps1, eps1},
    {-eps1, eps1, -eps1, eps1}, {-eps1, -eps1, -eps1, eps1},
    {eps1, eps1, eps1, -eps1}, {eps1, -eps1, eps1, -eps1},
    {eps1, eps1, -eps1, -eps1}, {eps1, -eps1, -eps1, -eps1},
    {-eps1, eps1, eps1, -eps1}, {-eps1, -eps1, eps1, -eps1},
    {-eps1, eps1, -eps1, -eps1}, {-eps1, -eps1, -eps1, -eps1},
  };

  std::vector<std::vector<double>> smallCut = {
    {0, 0, 0, 0}, {0, 0, 0, 1}, {0, 0, 0, -1},
    {-1, -1, -1, -3}, {-1, 1, -1, -3}, {1, 1, -1, -3}, {1, -1, -1, -3},
    {-1, -1, 1, -3}, {-1, 1, 1, -3}, {1, 1, 1, -3}, {1, -1, 1, -3},
    {1, 1, 1, 3}, {1, -1, 1, 3}, {-1, -1, 1, 3}, {-1, 1, 1, 3},
    {1, 1, -1, 3}, {1, -1, -1, 3}, {-1, -1, -1, 3}, {-1, 1, -1, 3},
    {-1, 0, 0, -1}, {1, 0, 0, - 1}, {0, -1, 0, -1}, {0, 1, 0, -1}, {0, 0, -1, -1}, {0, 0, 1, -1},
    {1, 0, 0, 1}, {-1, 0, 0, 1}, {0, 1, 0, 1}, {0, -1, 0, 1}, {0, 0, 1, 1}, {0, 0, -1, 1},
  };

 
  int s = 0;
  eps =1.0e-10;
  
  //std::cout << HyperCube(s,  m, af, df) << "\n";
  
  for(unsigned i = 0; i < smallCut.size(); i++) {
    for(unsigned j = 0; j < epsCut.size(); j++) {
      a.resize(3);
      af.resize(3); 
      a[0] = myType(smallCut[i][0] + epsCut[j][0]);
      a[1] = myType(smallCut[i][1] + epsCut[j][1]);
      a[2] = myType(smallCut[i][2] + epsCut[j][2]);
      d = myType(smallCut[i][3] + epsCut[j][3]);

      for(unsigned i = 0; i < a.size(); i++) {
        d -= a[i];
        a[i] *= 2;
      }

      af[0] = myTypeB(a[0]);
      af[1] = myTypeB(a[1]);
      af[2] = myTypeB(a[2]);
      df = myTypeB(d);

      myType I1 = HyperCube(s,  m, a, d);
      myTypeB I2 = HyperCube(s,  m, af, df);
      if((I2 != 0. &&  fabs((I1 - I2) / I2) < eps) || I1 == 0.) {
        //std::cout << "passed ";
        //std::cout << "s = " << s << " a = " << a[0] << " b = " << a[1] << " c = " << a[2] << " d = " << d << " I = " << I1 << std::endl;
      }
      else {
        std::cout << "Warning failed " << i << " " << j<<" ";
        std::cout << "s = " << s << " a = " << a[0] << " b = " << a[1] << " c = " << a[2] << " d = " << d << " I1 =" << I1 << " I2 = " << I2 << std::endl;
      }
    }
  }

}


