
#ifndef __femus_cut_fem_TRI_hpp__
#define __femus_cut_fem_TRI_hpp__

#include "Line.hpp"

template <class Type>
Type TriangleReduced(const int &s, const std::vector<unsigned> &m_input, const std::vector <Type> &a_input, const Type &d) {

  const Type &a = a_input[0];
  const Type &b = a_input[1];
  const unsigned &m = m_input[0];
  const unsigned &n = m_input[1];

  Type TRI = 0.;
  if(b == 0) { //parallel to right edge
    TRI = (-LimLi(s + 1, a + d) +
           LimLi(s + m + n + 2, d) * factorial<Type>(m + n + 1) * pow(-1 / a, m + n + 1)
          ) / (a * (n + 1));
  }
  else if(a == 0) { //parallel to bottom edge
    TRI = (LimLi(s + n + 1, d) * factorial<Type>(n) * pow(-1. / b, n) -
           LimLi(s + m + n + 2, d) * factorial<Type>(m + n + 1) * pow(-1. / b, m + n + 1)
          ) / (b * (m + 1.));
  }
  else if(a + b == 0) { // parallel to y=x edge//TODO
    if(a + d  >= 0.) {
      for(unsigned i = 1; i <= m + 1; i++) {
        TRI += LimLi(s + n + 1 + i, a + d) * pow(-1. / a, i) / factorial<Type>(m + 1 - i);
      }
      TRI *= factorial<Type>(n) * factorial<Type>(m) * pow(1. / a, n + 1);
    }
    TRI += LimLi(s + 1, d) / (a * (m + n + 1));
  }
  else { // general case
    if(fabs(b) > fabs(a)) {
      if(d > 0) {
        for(unsigned i = 1; i <= n + 1; i++) {
          TRI += factorial<Type>(m + n + 1 - i) / factorial<Type>(n + 1 - i) * pow((a + b) / b, i);
        }
        TRI *= LimLi(s + m + n + 2, d) / pow(-a - b, m + n + 2) ;
      }
      TRI += LSI(s + n + 1, m, a, d) * pow(-1 / b, n + 1);
      TRI *= factorial<Type>(n);
    }
    else {
      if(d > 0) {
        for(unsigned i = 1; i <= m + 1; i++) {
          TRI += factorial<Type>(m + n + 1 - i) / factorial<Type>(m + 1 - i) * pow((a + b) / a, i);
        }
        TRI *= - LimLi(s + m + n + 2, d) / pow(-a - b, m + n + 2) ;
      }
      if(a + d > 0.) {
        Type TRI2 = 0.;
        for(unsigned i = 1; i <= m + 1; i++) {
          TRI2 +=  LimLi(s + n + i + 1, a + d) / (factorial<Type>(m + 1 - i) * pow(-a, i));
        }
        TRI += TRI2 * factorial<Type>(n) / pow(-b, n + 1) ;
      }
      TRI *= factorial<Type>(m);
    }

  }
  return TRI;
}



template <class Type>
Type TriangleFull(const int &s, const std::vector<unsigned> &m_input, const std::vector <Type> &a_input, const Type &d) {

  const Type &a = a_input[0];
  const Type &b = a_input[1];
  const unsigned &m = m_input[0];
  const unsigned &n = m_input[1];

  Type TRI = 0;
  if(b == 0) { //parallel to right edge
    if(a == 0) TRI = -LimLi(s, d) / ((m + n + 2) * (n + 1));
    else TRI = LSI(s, m + n + 1, a, d) / (n + 1);
  }
  else if(a == 0) { //parallel to bottom edge
    TRI = (LSI(s, n, b, d) - LSI(s, m + n + 1, b, d)) / (m + 1);
  }
  else if(a + b == 0) { //parallel to y = x edge //TODO
    for(unsigned i = 1; i <= n + 1; i++) {
      TRI += LimLi(s + i, d) * pow(1. / a, i) / (factorial<Type>(n + 1 - i) * (m + n + 2 - i));
    }
    TRI += LSI(s + n + 1, m, a, d) / pow(a, n + 1);
    TRI *= factorial<Type>(n);
  }
  else { //generic case
    if(fabs(a) < fabs(b)) { // |a| < |b|
      for(unsigned j = 1; j <= n + 1; j++) {
        TRI -= LSI(s + j, m + n + 1 - j, a + b, d) * pow(-1 / b, j)
               / factorial<Type>(n + 1 - j);
      }
      TRI += LSI(s + n + 1, m, a, d) * pow(-1 / b, n + 1);
      TRI *= factorial<Type>(n);
    }
    else { // |b| < |a|
      for(unsigned j = 1; j <= m + 1; j++) {
        TRI += (LSI(s + j, m + n + 1 - j, a + b, d) -
                LSI(s + j, n, b, a + d)) * pow(-1 / a, j)
               / factorial<Type>(m + 1 - j);
      }
      TRI *= factorial<Type>(m);
    }
  }

  return TRI;
}



template <class Type>
Type TriangleA(const int &s, const std::vector<unsigned> &m, const std::vector <Type> &a, const Type &d) {

  if(a[0] == 0 && a[1] == 0) return -LimLi(s, d) / Type((m[0] + m[1] + 2) * (m[1] + 1)); // only from higher dimensions calls n = <0,0,c>

  Type sum = a[0] + a[1] + d;
  switch(s) {
    case -1:
      if(sum <= 0) {
        return TriangleReduced(-1, m, a, d);
      }
//       else if(sum <= std::max(fabs(a[1]), fabs(a[0]))) {
//         return TriangleFull(-1, m, a, d);
//       }
      else return TriangleReduced(-1, m, std::vector<Type> {-a[0], -a[1]}, -d);
      break;
    default:
      if(sum <= 0) {
        return TriangleReduced(s, m, a, d);
      }
      else if(sum <= std::max(fabs(a[1]), fabs(a[0]))) {
        return TriangleFull(s, m, a, d);
      }
      else {
        Type I1, I2;
        if(fabs(a[0]) <= fabs(a[1])) {
          Type m1f = factorial<Type>(m[1]);

          I1 = m1f / factorial<Type>(m[1] + s + 1) * pow(-a[1], s + 1)
               * TriangleReduced(-1, std::vector<unsigned> {m[0], m[1] + s + 1}, std::vector<Type> {-a[0], -a[1]}, -d);

          I2 = Type(0);
          for(int i = s; i >= 0; i--) {
            I2 += m1f / factorial<Type>(m[1] + 1 + i) * pow(-a[1], i)
                  * LSI(s - i, m[0] + m[1] + 1 + i, a[0] + a[1], d);
          }
        }
        else {
          Type m0f = factorial<Type>(m[0]);
          I1 = m0f / factorial<Type>(m[0] + s + 1) * pow(-a[0], s + 1)
               * TriangleReduced(-1, std::vector<unsigned> {m[0] + s + 1, m[1]}, std::vector<Type> {-a[0], -a[1]}, -d);
          I2 = Type(0);
          for(int i = 0; i <= s; i++) {
            I2 += m0f / factorial<Type>(m[0] + 1 + i) * pow(-a[0], i)
                  * (LSI(s - i, m[1], a[1], a[0] + d) -
                     LSI(s - i, m[0] + m[1] + 1 + i, a[0] + a[1], d));
          }

          //std::cout<<"@"<<m0f / factorial<Type>(m[0] + s + 1) * pow(-a[0], s + 1)<<" "<<TriangleReduced(-1, std::vector<unsigned> {m[0] + s + 1, m[1]}, std::vector<Type> {-a[0], -a[1]}, -d)<<"\n";


        }
        if(std::max(fabs(I1), fabs(I2)) == 0 || fabs(I1 / I2 + 1)  > 1.0e-1) { // || fabs(I1 + I2) > 1.0e-14) {
          return I1 + I2;
        }
        else  {

          Type I3 = TriangleFull(s, m, a, d);
          //std::cout << "@" << I1 << " " << I2 << " " << I3 << " " << sum << " " << a[0] << " " << a[1] << " " << d << std::endl;
          //std::cout << "!"<<fabs(I1 / I2 + 1)<<" ";// << " " << I1 + I2 << " " << I3 / (I1 + I2) << " " << I3 << std::endl;
          return I3;
        }
      }
  }
}

template <class Type>
Type Triangle(const int &s, const std::vector<unsigned> &m, const std::vector <Type> &a, const Type &d) {
  return TriangleA(s, m,  std::vector<Type> {-a[0], a[1]}, d + a[0]);
}

#endif

