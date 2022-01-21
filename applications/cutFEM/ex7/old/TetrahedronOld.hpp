#ifndef __femus_cut_fem_TET_Old_hpp__
#define __femus_cut_fem_TET__Old_hpp__

#include "TriangleOld.hpp"

template <class Type>
Type TetrahedronA(const int &s, const std::vector<unsigned> &m, const std::vector <Type> &a, const Type & d);

template <class Type>
Type TetrahedronB(const int &s, const std::vector<unsigned> &m_input, const std::vector <Type> &a_input, const Type & d);

template <class Type>
Type TetrahedronB(const int &s, const std::vector<unsigned> &m_input, const std::vector <Type> &a_input, const Type & d) {

  const Type &a = a_input[0];
  const Type &b = a_input[1];
  const Type &c = a_input[2];
  const unsigned &m = m_input[0];
  const unsigned &n = m_input[1];
  const unsigned &o = m_input[2];

  Type TET = 0;
  if(fabs(c) > fabs(b)) {
    for(unsigned i = 1; i <= o + 1; i++)  {
      TET -= TriangleA<Type>(s + i, {m, n + o + 1u - i}, {a, b + c}, d) / (factorial<Type> (o + 1u - i) * pow(-c, i));
    }
    TET += TriangleA<Type>(s + o + 1, {m, n}, {a, b}, d) / pow(-c, o + 1);
    TET *= factorial<Type>(o);
  }
  else {
    for(unsigned i = 1; i <= n + 1; i++)  {
      TET += (-TriangleA<Type>(s + i, {m + n + 1u - i, o}, {a + b, c}, d) +
              TriangleA<Type>(s + i, {m, n + o + 1u - i}, {a, b + c}, d)) / (factorial<Type> (n + 1u - i) * pow(-b, i));
    }
    TET *= factorial<Type>(n);
  }
  return TET;
}

template <class Type>
Type TetrahedronC(const int &s, const std::vector<unsigned> &m_input, const std::vector <Type> &a_input, const Type & d) {

  const Type &a = a_input[0];
  const Type &b = a_input[1];
  const Type &c = a_input[2];
  const unsigned &m = m_input[0];
  const unsigned &n = m_input[1];
  const unsigned &o = m_input[2];

  Type TET = 0;
  if(fabs(c) > fabs(b)) {
    for(unsigned i = 0; i <= s; i++)  {
      TET += TriangleA<Type>(s - i, {m, n + o + i + 1u}, {a, b + c}, d) * pow(-c, i) / factorial<Type>(o + i + 1);
    }
    TET += TetrahedronA<Type>(-1, {m, n, o + s + 1u}, a_input, d) * pow(-c, s + 1u) / factorial<Type>(o + s + 1);

    TET *= factorial <Type> (o);
  }
  else {
    for(unsigned i = 0; i <= s; i++)  {
      TET += (TriangleA<Type>(s - i, {m + n + i + 1u, o}, {a + b, c}, d)
              - TriangleA<Type>(s - i, {m, n + o + i + 1u}, {a, b + c}, d)) * pow(-b, i) / factorial<Type>(n + i + 1);
    }
    TET += TetrahedronA<Type>(-1, {m, n + s + 1u, o}, a_input, d) * pow(-b, s + 1u) / factorial<Type>(n + s + 1);

    TET *= factorial <Type> (n);
  }


  return TET;
}


template <class Type>
Type TetrahedronA(const int &s, const std::vector<unsigned> &m, const std::vector <Type> &a, const Type & d) {
   
  switch(s) {
    case -1:
      if(a[0] + a[1] + a[2] + d <= 0) return TetrahedronB(-1, m, a, d);
      else return TetrahedronB<Type>(-1, m, {-a[0], -a[1], -a[2]}, -d);
      break;
    default:
      if(a[0] + a[1] + a[2] + d <= std::max(fabs(a[1]), fabs(a[2]))) {
        return TetrahedronB<Type>(s, m, a, d);
      }
      else {
        return TetrahedronC<Type>(s, m, a, d);
      }
  }
}


template <class Type, class Type2>
Type Tetrahedron(const int &s, const std::vector<unsigned> &m, const std::vector <Type> &a, const Type & d) {

  Type m1 = std::max(fabs(a[1] + a[0]) , fabs(a[2] - a[1]) );  
  Type m2 = std::max(fabs(a[2] + a[1]) , fabs(a[0] - a[2]) ); 
  Type m3 = std::max(fabs(a[0] + a[2]) , fabs(a[1] - a[0]) );  
    
  if( m1 > m2 && m1 > m3 ) {
    return static_cast<Type>(TetrahedronA<Type2>(s, {m[0], m[1], m[2]}, {a[0], a[1] + a[0], a[2] - a[1]}, d));
  }
  else if(m2 > m3) {
    return static_cast<Type>(TetrahedronA<Type2>(s, {m[1], m[2], m[0]}, {a[1], a[2] + a[1], a[0] - a[2]}, d));
  }
  else {
    return static_cast<Type>(TetrahedronA<Type2>(s, {m[2], m[0], m[1]}, {a[2], a[0] + a[2], a[1] - a[0]}, d));
  }

}





#endif
