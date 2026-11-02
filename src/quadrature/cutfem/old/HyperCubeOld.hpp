
#ifndef __femus_cut_fem_HCI_old_hpp__
#define __femus_cut_fem_HCI_old_hpp__

#include "LineOld.hpp"

template <class Type>
Type HyperCubeA(const int &s, std::vector<unsigned> m, std::vector <Type> a, const Type &d);

template <class Type>
Type HyperCubeA1(const unsigned & n, const int &s, std::vector<unsigned> &m,
                const std::vector <Type> &a, const std::vector <Type> &ma,
                const Type & d, const Type & md);

template <class Type>
Type HyperCubeB(const unsigned & n, const int &s, std::vector<unsigned> &m,
                const std::vector <Type> &a, const std::vector <Type> &ma,
                const Type & d, const Type & md);

template <class Type>
Type HyperCubeC(const unsigned & n, const int &s, std::vector<unsigned>& m,
                const std::vector <Type> &a, const std::vector <Type> &ma,
                const Type & d, const Type & md);

template <class Type, class Type2>
Type HyperCube(const int &s, const std::vector<unsigned> &m, const std::vector <Type> &a, const Type &d) {
  std::vector<Type2> a2(a.size());
  Type d1 = d;
  for(unsigned i = 0; i < a.size(); i++) {
    d1 -= a[i];
    a2[i] = static_cast<Type2>(2 * a[i]);
  }
  Type2 d2 = static_cast<Type2>(d1);
  return static_cast<Type>(pow(Type(2), a2.size()) * HyperCubeA<Type2>(s, m, a2, d2));
}

template <class Type>
Type HyperCubeA(const int &s, std::vector<unsigned> m, std::vector <Type> a, const Type &d) {

  Type HCI = 1;
  
  for(int i = a.size() - 1; i >= 0; i--) {
    if(a[i] == 0) {
      HCI *= 1 / Type(m[i] + 1);
      a.erase(a.begin() + i);
      m.erase(m.begin() + i);
    }
  }

  if(a.size() > 0) {
    // all the left a[i] coefficients \ne 0
    for(unsigned i = 0; i < a.size() - 1; i++) { // reorder iterated integral from the smallest to the largest coefficient a[i]
      for(unsigned j = i + 1; j < a.size(); j++) {
        if(fabs(a[i]) > fabs(a[j])) {
          std::swap(a[i], a[j]);
          std::swap(m[i], m[j]);
        }
      }
    }

    std::vector <Type> ma(a.size());
    for(unsigned i = 0; i < a.size(); i++)  ma[i] = -a[i];

    return HCI * HyperCubeA1(-1 + a.size(), s, m, a, ma, d, -d);
  }
  else {
    return -HCI * LimLi(s, d);
  }
}

template <class Type>
Type HyperCubeA1(const unsigned & n, const int &s, std::vector<unsigned> &m,
                const std::vector <Type> &a, const std::vector <Type> &ma,
                const Type & d, const Type & md) {

  if(n == 0)  {
    return LSI(s, m[0], a[0], d);
  }

  Type sum = d;
  for(unsigned i = 0; i <= n; i++) sum += a[i];

  switch(s) {
    case -1: // interface integral
      if(sum <= 0) {
        return HyperCubeB(n, -1, m, a, ma, d, md);
      }
      else {
        return HyperCubeB(n, -1, m, ma, a, md, d);
      }
      break;
    default:
      if(sum <= fabs(a[n])) {
        return HyperCubeB(n, s, m, a, ma, d, md);
      }
      else {
        return HyperCubeC(n, s, m, a, ma, d, md);
      }
  }
}

template <class Type>
Type HyperCubeB(const unsigned &n, const int &s, std::vector<unsigned> &m,
                const std::vector <Type> &a, const std::vector <Type> &ma,
                const Type &d, const Type &md) {

  Type aI = 1 / a[n];
  Type c = aI; // this is m!/(m+1-j)! 1/a^j for j = 1,...,m + 1

  Type I1 = 0;
  for(int j = 1; j <= m[n] + 1;  c *= -aI * (m[n] + 1 - j), j++) {
    I1 += HyperCubeA1(n - 1, s + j, m, a, ma, d + a[n], -d - a[n]) * c;
  }
  Type I2 = HyperCubeA1(n - 1, s + m[n] + 1, m, a, ma, d, -d) * factorial<Type>(m[n]) * pow(-aI, m[n] + 1);
  return I1 + I2;

}

template <class Type>
Type HyperCubeC(const unsigned & n, const int &s, std::vector<unsigned>& m,
                const std::vector <Type> &a, const std::vector <Type> &ma,
                const Type & d, const Type & md) { //alternative formula at s = -1

  Type an = a[n];
  unsigned mn = m[n];
  
  m[n] += (s + 1);
  Type I1 = HyperCubeA1(n, -1, m, a, ma, d, md) * pow(-an, s + 1) * factorial<Type>(mn) / factorial<Type>(mn + s + 1) ;
  m[n] -= (s + 1);

  Type I2 = 0;
  Type c1 = 1 / Type(mn + 1);
  for(unsigned i = 0; i <= s; c1 *= -an / (mn + 2 + i), i++) {
    I2 += HyperCubeA1(n - 1, s - i, m, a, ma, d + an, -(d + an)) * c1;
  }
  return I1 + I2;

}

#endif
