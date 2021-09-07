
#ifndef __femus_cut_fem_HCI_hpp__
#define __femus_cut_fem_HCI_hpp__

#include "Line.hpp"

template <class Type>
Type HyperCube(const int &s, std::vector<unsigned> m, std::vector <Type> a, Type d);

template <class Type>
Type HyperCubeA(const unsigned & n, const int &s, std::vector<unsigned> &m,
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

template <class Type>
Type HyperCube(const int &s, std::vector<unsigned> m, std::vector <Type> a, Type d) {

//   for(unsigned i = 0; i < a.size(); i++) {
//     d -= a[i];
//     a[i] *= 2;
//   }

  Type HCI = pow(Type(2), a.size());


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

    return HCI * HyperCubeA(-1 + a.size(), s, m, a, ma, d, -d);
  }
  else {
    return -HCI * LimLi(s, d);
  }
}

template <class Type>
Type HyperCubeA(const unsigned & n, const int &s, std::vector<unsigned> &m,
                const std::vector <Type> &a, const std::vector <Type> &ma,
                const Type & d, const Type & md) {

  //std::cout << "A";
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
        //std::cout << "!";
        return HyperCubeB(n, s, m, a, ma, d, md);
      }
      else {
        //std::cout << "*";
        return HyperCubeC(n, s, m, a, ma, d, md);
      }
  }
}

template <class Type>
Type HyperCubeB(const unsigned & n, const int &s, std::vector<unsigned> &m,
                const std::vector <Type> &a, const std::vector <Type> &ma,
                const Type & d, const Type & md) {

  Type aI = 1. / a[n];
  Type c = aI; // this is m!/(m+1-j)! 1/a^j for j = 1,...,m + 1

  Type I1 = 0;

  for(int j = 1; j <= m[n] + 1;  c *= -aI * (m[n] + 1 - j), j++) {
    I1 += HyperCubeA(n - 1, s + j, m, a, ma, d + a[n], -d - a[n]) * c;
  }

  Type I2 = HyperCubeA(n - 1, s + m[n] + 1, m, a, ma, d, -d) * factorial<Type>(m[n]) * pow(-aI, m[n] + 1);

  return I1 + I2;


}

template <class Type>
Type HyperCubeC(const unsigned & n, const int &s, std::vector<unsigned>& m,
                const std::vector <Type> &a, const std::vector <Type> &ma,
                const Type & d, const Type & md) { //alternative formula at s = -1

  Type HCI = 0.;
  Type an = a[n];

  unsigned mn = m[n];
  Type mnf = factorial<Type>(mn);
  m[n] += (s + 1);
  Type I1 = HyperCubeA(n, -1, m, a, ma, d, md) * mnf * pow(-an, s + 1) / factorial<Type>(mn + s + 1) ;
  m[n] -= (s + 1);

  Type I2 = 0;
  Type c1 = 1 / Type(mn + 1);

  for(unsigned i = 0; i <= s; c1 *= -an / (mn + 2 + i), i++) {
    I2 += HyperCubeA(n - 1, s - i, m, a, ma, d + an, -(d + an)) * c1;
  }
  return I1 + I2;
//   if(std::max(fabs(I1), fabs(I2)) == 0 || fabs(I1 / I2 + 1)  > 1.0e-3) {
//     return I1 + I2;
//   }
//   else  {
//     std::cout << "?";
//     return HyperCubeB(n, s, m, a, ma, d, md);
//   }

}



//These below are alternative formulations. They work too!


template <class Type>
Type HyperCubeC0(const unsigned & n, const int &s, std::vector<unsigned>& m,
                 const std::vector <Type> &a, const std::vector <Type> &ma,
                 const Type & d, const Type & md) { //alternative formula at s = 0

  //std::cout << "C";
  Type HCI = 0.;
  Type an = a[n];

  unsigned mn = m[n];
  Type mnf = factorial<Type>(mn);
//   m[n] += (s + 1);
//   Type I1 = mnf / factorial<Type>(mn + s + 1) * pow(-an, s + 1) * HyperCubeA(n, -1, m, a, ma, d, md) ;
//   m[n] -= (s + 1);

  m[n] += s;
  Type I1 = mnf / factorial<Type>(mn + s) * pow(-an, s) * HyperCubeA(n, 0, m, a, ma, d, md) ;
  m[n] -= s;

  Type I2 = 0;
  Type c1 = 1 / Type(mn + 1);

  for(unsigned i = 0; i < s; c1 *= -an / (mn + 2 + i), i++) {
    I2 += c1 * HyperCubeA(n - 1, s - i, m, a, ma, d + an, -(d + an));
  }
  if(std::max(fabs(I1), fabs(I2)) == 0 || fabs(I1 / I2 + 1)  > 1.0e-3) {
    return I1 + I2;
  }
  else  {
    std::cout << "?";
    return HyperCubeB(n, s, m, a, ma, d, md);
  }

}




template <class Type>
Type HyperCubeA0(const unsigned & n, const int &s, std::vector<unsigned> &m,
                 const std::vector <Type> &a, const std::vector <Type> &ma,
                 const Type & d, const Type & md) {

  //std::cout << "A";
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
    case 0: // step function integral
      if(sum <= fabs(a[n])) {
        //std::cout << "%";
        return HyperCubeB(n, 0, m, a, ma, d, md);
      }
      else {
        //std::cout << "^";
        Type I1 = 1;
        for(unsigned i = 0; i <= n; i++) {
          I1 *= 1 / Type(m[i] + 1);
        }
        Type I2 = -HyperCubeB(n, 0, m, ma, a, md, d);
        if(std::max(fabs(I1), fabs(I2)) == 0 || fabs(I1 / I2 + 1)  > 1.0e-3) {
          return I1 + I2;
        }
        else  {
          std::cout << "@";
          return HyperCubeB(n, 0, m, a, ma, d, md);
        }
      }
      break;
    default: // all other cases, whatever they mean
      if(sum <= fabs(a[n])) {
        //std::cout << "!";
        return HyperCubeB(n, s, m, a, ma, d, md);
      }
      else {
        //std::cout << "*";
        return HyperCubeC(n, s, m, a, ma, d, md);
      }
  }
}

#endif
