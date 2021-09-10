#ifndef __femus_cut_fem_TET_hpp__
#define __femus_cut_fem_TET_hpp__

#include "Triangle.hpp"

unsigned cast = 0;
unsigned uncast = 0;


template <class Type>
Type TetrahedronA(const int &s, const std::vector<unsigned> &m, const std::vector <Type> &a, const Type & d);

template <class Type>
Type TetrahedronCast(const int &s, const std::vector<unsigned> &m, const std::vector <Type> &a, const Type & d) {
  if(false && //(a[2] < 0 && fabs((a[0] + a[1]) * (fabs(a[0] / a[2]) + fabs(a[1] / a[2]) + fabs(d / a[2])))  < 1.0e-10) ||
      (a[2] > 0 && (fabs(a[0] / a[2]) + fabs(a[1] / a[2]) + fabs((a[2] + d) / a[2]) < 1.0e-5))) {
    //std::cout << "c";
    cast++;
    cpp_bin_float_oct df = d;
    return static_cast<Type>(TetrahedronA< cpp_bin_float_oct >(s, m, {a[0], a[1], a[2]}, df));
  }
  else {
    //std::cout << "u";
    uncast++;
    return TetrahedronA< Type >(s, m, a, d);
  }
};





template <class Type>
Type TetrahedronB(const int &s, const std::vector<unsigned> &m_input, const std::vector <Type> &a_input, const Type & d);

template <class Type>
Type TetrahedronB(const int &s, const std::vector<unsigned> &m_input, const std::vector <Type> &a_input, const Type & d) {

  //std::cout << "B ";

  const Type &a = a_input[0];
  const Type &b = a_input[1];
  const Type &c = a_input[2];
  const unsigned &m = m_input[0];
  const unsigned &n = m_input[1];
  const unsigned &o = m_input[2];

  Type TET = 0.;
  if(fabs(c) >= fabs(b)) {
    for(unsigned i = 1; i <= o + 1; i++)  {
      TET -= TriangleA<Type>(s + i, {m, n + o + 1u - i}, {a, b + c}, d) / (factorial<Type> (o + 1u - i) * pow(-c, i));
    }
    TET += TriangleA<Type>(s + o + 1, {m, n}, {a, b}, d) / pow(-c, o + 1);
    TET * factorial<Type>(o);
  }
  else {
    std::cout<<"^";  
    for(unsigned i = 1; i <= n + 1; i++)  {
      TET = (-TriangleA<Type>(s + i, {m + n + 1u - i, o}, {a + b, c}, d) +
             TriangleA<Type>(s + i, {m, n + o + 1u - i}, {a, b + c}, d)) / (factorial<Type> (n + 1u - i) * pow(-b, i));
    }
  }
  return TET;
}

template <class Type>
Type TetrahedronC(const int &s, const std::vector<unsigned> &m_input, const std::vector <Type> &a_input, const Type & d) {

  //std::cout << "C ";

  const Type &a = a_input[0];
  const Type &b = a_input[1];
  const Type &c = a_input[2];
  const unsigned &m = m_input[0];
  const unsigned &n = m_input[1];
  const unsigned &o = m_input[2];

  Type TET = 0;
  for(unsigned i = 1; i <= s + 1; i++)  {
    TET += TriangleA<Type>(s + 1 - i, {m, n + o + i}, {a, b + c}, d) * pow(-c, i - 1) / factorial<Type>(o + i);
  }
  TET += TetrahedronA<Type>(-1, {m, n, o + s + 1u}, a_input, d) * pow(-c, s + 1u) / factorial<Type>(o + s + 1);

  TET *= factorial <Type> (o);

  return TET;
}


template <class Type>
Type TetrahedronA(const int &s, const std::vector<unsigned> &m, const std::vector <Type> &a, const Type & d) {

//   if(fabs(a[2]) < fabs(a[1])) {
//     std::swap(a[2], a[1]);
//     std::swap(m[2], m[1]);
//   }

  switch(s) {
    case -1:
      if(a[0] + a[1] + d <= 0) return TetrahedronB(-1, m, a, d);
      else return TetrahedronB<Type>(-1, m, {-a[0], -a[1], -a[2]}, -d);
      break;
    default:
      if(a[0] + a[1] + d <= fabs(a[2])) {
        return TetrahedronB<Type>(s, m, a, d);
      }
      else {
        return TetrahedronC<Type>(s, m, a, d);
      }
  }
}


template <class Type>
Type Tetrahedron(const int &s, const std::vector<unsigned> &m, const std::vector <Type> &a, const Type & d) {

  if(fabs(a[0]) <= fabs(a[1]) && fabs(a[0]) <= fabs(a[2])) {
    if(fabs(a[1]) <= fabs(a[2]))
      return TetrahedronCast<Type>(s, {m[0], m[1], m[2]}, {-a[0], a[1], a[2]}, d + a[0]) ;
    else
      return TetrahedronCast<Type>(s, {m[0], m[2], m[1]}, {-a[0], a[2], a[1]}, d + a[0]) ;
  }
  else if(fabs(a[1]) <= fabs(a[2])) {
    if(fabs(a[0]) <= fabs(a[2]))
      return TetrahedronCast<Type>(s, {m[1], m[0], m[2]}, {-a[1], a[0], a[2]}, d + a[1]) ;
    else
      return TetrahedronCast<Type>(s, {m[1], m[2], m[0]}, {-a[1], a[2], a[0]}, d + a[1]) ;
  }
  else {
    if(fabs(a[0]) <= fabs(a[1]))
      return TetrahedronCast<Type>(s, {m[2], m[0], m[1]}, {-a[2], a[0], a[1]}, d + a[2]) ;
    else
      return TetrahedronCast<Type>(s, {m[2], m[1], m[0]}, {-a[2], a[1], a[0]}, d + a[2]) ;
  }

}

















template <class Type>
Type TetrahedronB1(const int &s, const std::vector<unsigned> &m_input, const std::vector <Type> &a_input, const Type & d) {

  std::cout << "B1 ";

  const Type &a = a_input[0];
  const Type &b = a_input[1];
  const Type &c = a_input[2];
  const unsigned &m = m_input[0];
  const unsigned &n = m_input[1];
  const unsigned &o = m_input[2];

  Type S0 = TriangleA(s + o + 1, std::vector<unsigned> {m, n}, std::vector<Type> {a, b}, d) / pow(-c, o + 1);

  Type fn = factorial<Type>(n);
  Type cMb = c - b;
  Type cMbOc = (c - b) / c;
  Type bMcOc = -cMbOc;

  Type S1 = Type(0);
  Type fi = 1 / fn;

  Type fj = fn;
  Type Sj = Type(0);
  Type cMbOci = Type(1);
  for(unsigned i = 1; i <= o + 1; i++) {
    fi = -fi / (n + i);
    cMbOci *= cMbOc;
    Sj += fj;
    S1 -= fi * TriangleA(s, std::vector<unsigned> {m + o + 1 - i, n + i}, std::vector<Type> {a + c, b - c}, d) /
          factorial<Type>(o + 1 - i) * cMbOci * Sj;
    fj = fj / cMbOc * (n + i) / i;
  }

  Type S2 = 0;
  for(unsigned i = 0; i <= o ; i++) {
    Type S2i = 0.;
    for(unsigned k = 1; k <= o + 1 - i; k++) {
      Type S2ik = 0.;
      Type f =  fn / factorial<Type>(n + k);
      for(unsigned j = 0; j <= o + 1 - i - k; j++) {
        //std::cout << i << " " << k << " " << j << " " <<  o + 1 - i - k << std::endl;
        S2ik += f / factorial<Type>(o + 1 - i - k - j);
        f = - f * Type(n + j + 1) / Type((j + 1) * (n + k + j + 1)) ;
      }
      S2i += pow(bMcOc, k - 1) * S2ik;
    }
    S2 -= S2i * LSI(s + 1 + i, m + n + o + 1 - i, a + b, d) / pow(-c, i);
  }

  return (S0 + S1 + S2) * factorial<Type>(o);

}





#endif
