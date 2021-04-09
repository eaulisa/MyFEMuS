#include <boost/multiprecision/cpp_dec_float.hpp>
#include <boost/multiprecision/cpp_bin_float.hpp>

using boost::multiprecision::cpp_dec_float_50;
using boost::multiprecision::cpp_dec_float_100;

namespace boost {
  namespace multiprecision {
    typedef number<cpp_dec_float<200> > cpp_dec_float_200;

    typedef number < backends::cpp_bin_float < 24, backends::digit_base_2, void, boost::int16_t, -126, 127 >, et_off >         cpp_bin_float_single;
    typedef number < backends::cpp_bin_float < 53, backends::digit_base_2, void, boost::int16_t, -1022, 1023 >, et_off >       cpp_bin_float_double;
    typedef number < backends::cpp_bin_float < 64, backends::digit_base_2, void, boost::int16_t, -16382, 16383 >, et_off >     cpp_bin_float_double_extended;
    typedef number < backends::cpp_bin_float < 113, backends::digit_base_2, void, boost::int16_t, -16382, 16383 >, et_off >    cpp_bin_float_quad;
    typedef number < backends::cpp_bin_float < 237, backends::digit_base_2, void, boost::int32_t, -262142, 262143 >, et_off >  cpp_bin_float_oct;
  }
} // namespaces

//typedef double myType;
//   typedef double myType;
//typedef boost::multiprecision::cpp_bin_float_oct myType;


//typedef cpp_dec_float_50 myType;
//typedef cpp_dec_float_100 myType;
//typedef boost::multiprecision::cpp_dec_float_200 myType;
//typedef boost::multiprecision::cpp_bin_float_double myType;
//typedef boost::multiprecision::cpp_bin_float_double_extended myType; //long double
//typedef boost::multiprecision::cpp_bin_float_quad myType;
//typedef long double myTypeB;


template <class Float1, class Float2>
typename boost::math::tools::promote_args<Float1, Float2>::type
Int0to1LimLiA(const int &s, const unsigned &m, const Float1 &a, const Float2 &d) {
  //this function can be called only if s>0, (a+d) and d are both non-negative
  typedef typename boost::math::tools::promote_args<Float1, Float2>::type Type;
  Type INT(0.);
  Type x(a + d);

  for(unsigned i = 1; i <= s; i++) {
    INT -= pow(-a, s - i) / factorial<Type>(m + 1 + s - i) * LimLi(i, x) ;
  }
  INT += pow(-a, s) / factorial<Type>(m + 1 + s);
  INT *= factorial<Type>(m);
  return INT;

}

template <class Float1, class Float2>
typename boost::math::tools::promote_args<Float1, Float2>::type
Int0to1LimLiB(const int &s, const unsigned &m, const Float1 &a, const Float2 &d) {
  //this function can be called for any s, (a+d) and d, however if they are both non-negative and a->0 it suffers of significant digit cancellation
  //for s > 0, and (a+d) and d non-negative, we then call the alternative version Int0to1LimLiA
  //for s=-1 and s=0, we need to be sure that d is negative, so that cancellation won't occurs, this is done apriori.
  typedef typename boost::math::tools::promote_args<Float1, Float2>::type Type;

  Type INT = 0.;
  Type x(a + d);
  if(x >= 0.) {
    Type g =  1. / (-a);
    for(unsigned i = 1; i <= m + 1; i++) {
      INT += g * LimLi(s + i, x);
      g *= (m + 1 - i) / (-a);
    }
  }
  INT += pow(-1. / a, m) * LimLi(s + m + 1, d) / a * factorial<Type>(m);
  return INT;
}

template <class Float1, class Float2>
typename boost::math::tools::promote_args<Float1, Float2>::type
Int0to1LimLiC(const int &s, const unsigned &m, const Float1 &a, const Float2 &d) {

  if(s > 0 && d >= 0 && a + d >= 0) {
    return Int0to1LimLiA(s, m, a, d);
  }
  else {
    return Int0to1LimLiB(s, m, a, d);
  }
}


template <class Float1, class Float2>
typename boost::math::tools::promote_args<Float1, Float2>::type
Intm1to1LimLiA(const int &s, const unsigned &m, const Float1 &a, const Float2 &d) {
  //this function can be called only if s > 0, (a+d) and (-a+d) are both non-negative

  typedef typename boost::math::tools::promote_args<Float1, Float2>::type Type;
  Type INT(0.);
  Type x(a + d);
  Type y(-a + d);

  for(unsigned i = 1; i <= s; i++) {
    INT += pow(-a, s - i) / factorial<Type>(m + 1 + s - i) *
           (pow(x, i) + pow(-1, m + s - i) * pow(y, i)) / factorial<Type>(i);
  }
  INT += pow(-a, s) * (1 + pow(-1, m + s)) / factorial<Type>(m + 1 + s);
  INT *= factorial<Type>(m);
  return INT;

}

template <class Float1, class Float2>
typename boost::math::tools::promote_args<Float1, Float2>::type
Intm1to1LimLiB(const int &s, const unsigned &m, const Float1 &a, const Float2 &d) {
  //this function can be called for any s, (a+d) and (-a+d), however if they are both non-negative and a->0 it suffers of significant digit cancellation
  //for s > 0, and (a+d) and (-a+d) non-negative, we then call the alternative version Intm1to1LimLiA
  //for s=-1 and s=0, we need to be sure that d is negative, so that cancellation won't occurs, this is done apriori.
  typedef typename boost::math::tools::promote_args<Float1, Float2>::type Type;
  Type INT = 0.;
  Type x(a + d);
  if(x >= 0.) {
    Type g =  1. / (-a);
    for(unsigned i = 1; i <= m + 1; i++) {
      INT += g * LimLi(s + i, x);
      g *= (m + 1 - i) / (-a);
    }
  }
  x = (-a + d);
  if(x >= 0.) {
    Type g =  pow(-1., m) / a;
    for(unsigned i = 1; i <= m + 1; i++) {
      INT += g * LimLi(s + i, x);
      g *= (m + 1 - i) / a;
    }
  }
  return INT;
}

template <class Float1, class Float2>
typename boost::math::tools::promote_args<Float1, Float2>::type
Intm1to1LimLiC(const int &s, const unsigned &m, const Float1 &a, const Float2 &d) {

  if(s > 0 && -a + d >= 0 && a + d >= 0) {
    return Intm1to1LimLiA(s, m, a, d);
  }
  else {
    return Intm1to1LimLiB(s, m, a, d);
  }
}

