#ifndef __femus_cut_fem_PRI_hpp__
#define __femus_cut_fem_PRI_hpp__

#include "./old/TriangleOld.hpp"
#include "./old/HyperCubeOld.hpp"

template <class Type>
Type PrismA(const int &s, const std::vector<unsigned> &m, const std::vector <Type> &a, const Type & d);

template <class Type>
Type PrismB(const int &s, const std::vector<unsigned> &m_input, const std::vector <Type> &a_input, const Type & d) {

  const Type &a = a_input[0];
  const Type &b = a_input[1];
  const Type &c = a_input[2];
  const unsigned &m = m_input[0];
  const unsigned &n = m_input[1];
  const unsigned &o = m_input[2];

  Type PRI = 0;
  if(fabs(c) >= std::max(fabs(a), fabs(b))) {
    //std::cout << "^";
    for(unsigned i = 1; i <= o + 1; i++)  {
      PRI -= TriangleA<Type>(s + i, {m, n}, {a, b}, d + c) / (factorial<Type> (o + 1u - i) * pow(-c, i));
    }
    PRI += TriangleA<Type>(s + o + 1, {m, n}, {a, b}, d) / pow(-c, o + 1);
    PRI *= factorial<Type>(o);
  }
  else if(fabs(b) >= fabs(a)) {
    //std::cout << "!";
    for(unsigned i = 1; i <= n + 1; i++)  {
      PRI -= HyperCubeA<Type>(s + i, {m + n + 1u - i, o}, {a + b, c}, d) / (factorial<Type> (n + 1u - i) * pow(-b, i));
    }
    PRI += HyperCubeA<Type>(s + n + 1, {m, o}, {a, c}, d) / pow(-b, n + 1);
    PRI *= factorial<Type>(n);
  }
  else {
    //std::cout << "*";
    for(unsigned i = 1; i <= m + 1; i++)  {
      PRI += (-HyperCubeA<Type>(s + i, {n, o}, {b, c}, d + a)
              + HyperCubeA<Type>(s + i, {m + n + 1u - i, o}, {a + b, c}, d)) / (factorial<Type> (m + 1u - i) * pow(-a, i));
    }
    PRI *= factorial<Type>(m);
  }
  return PRI;
}

template <class Type>
Type PrismC(const int &s, const std::vector<unsigned> &m_input, const std::vector <Type> &a_input, const Type & d) {

  const Type &a = a_input[0];
  const Type &b = a_input[1];
  const Type &c = a_input[2];
  const unsigned &m = m_input[0];
  const unsigned &n = m_input[1];
  const unsigned &o = m_input[2];

  Type PRI = 0;
  if(fabs(c) >= std::max(fabs(a), fabs(b))) {
    //std::cout << "@";
    for(unsigned i = 0; i <= s; i++)  {
      PRI += TriangleA<Type>(s - i, {m, n}, {a, b}, d + c) * pow(-c, i) / factorial<Type>(o + i + 1);
    }
    PRI += PrismA<Type>(-1, {m, n, o + s + 1u}, a_input, d) * pow(-c, s + 1u) / factorial<Type>(o + s + 1);

    PRI *= factorial <Type> (o);
  }
  else if( fabs(b) >= fabs(a)) {
    //std::cout << "#";
    for(unsigned i = 0; i <= s; i++)  {
      PRI += HyperCubeA<Type>(s - i, {m + n + i + 1u, o}, {a + b, c}, d) * pow(-b, i) / factorial<Type>(n + i + 1);
    }
    PRI += PrismA<Type>(-1, {m, n + s + 1u, o}, a_input, d) * pow(-b, s + 1u) / factorial<Type>(n + s + 1);

    PRI *= factorial <Type> (n);
  }
  else {
    //std::cout << "%";
    for(unsigned i = 0; i <= s; i++)  {
      PRI += (HyperCubeA<Type>(s - i, {n, o}, {b, c}, d + a)
              - HyperCubeA<Type>(s - i, {m + n + i + 1u, o}, {a + b, c}, d)) * pow(-a, i) / factorial<Type>(m + i + 1);
    }
    PRI += PrismA<Type>(-1, {m + s + 1u, n, o}, a_input, d) * pow(-a, s + 1u) / factorial<Type>(m + s + 1);
    PRI *= factorial <Type> (m);
  }
  return PRI;
}

template <class Type>
Type PrismA(const int &s, const std::vector<unsigned> &m, const std::vector <Type> &a, const Type & d) {

  switch(s) {
    case -1:
      if(a[0] + a[1] + a[2] + d <= 0) return PrismB(-1, m, a, d);
      else return PrismB<Type>(-1, m, {-a[0], -a[1], -a[2]}, -d);
      break;
    default:
      if( a[0] + a[1] + a[2] + d <= std::max(std::max(fabs(a[0]), fabs(a[1])), fabs(a[2]))) {
        return PrismB<Type>(s, m, a, d);
      }
      else {
        return PrismC<Type>(s, m, a, d);
      }
  }
}

template <class Type, class Type2>
Type Prism(const int &s, const std::vector<unsigned> &m, const std::vector <Type> &a, const Type & d) {
  return static_cast<Type>(/*2 **/ PrismA<Type2>(s, {m[0], m[1], m[2]}, {-a[0], a[1], 2 * a[2]}, d + a[0] - a[2]));
}

#endif
