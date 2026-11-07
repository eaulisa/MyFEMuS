#ifndef FEMUS_PLANAR_MOMENTS_LS_HPP
#define FEMUS_PLANAR_MOMENTS_LS_HPP

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <utility>
#include <vector>
#include <boost/math/special_functions/factorials.hpp>

namespace CutFemWeightLSDetail {

template <class TypeA>
class CutFEMmap {
  public:
    virtual void clear() = 0;
    virtual TypeA operator()(const int &s, const std::vector<unsigned> &m, const std::vector<TypeA> &a, const TypeA &d) = 0;
    virtual ~CutFEMmap(){};
    virtual void SetBaseType (const unsigned &value){
    };
};

template <class Type>
class LimLimap {

  protected:
    LimLimap() {

    };
    ~LimLimap() {

    };

    void clear() {

    };

    Type limLi(const int &s, const Type &d) {

      return this->LimLi(s, d);

    }

  private:
    Type LimLi(const int &s, const Type &x);

};

template <class Type>
Type LimLimap<Type>::LimLi(const int &s, const Type &x) {

  if(x > 0) return -pow(x, s) / boost::math::factorial<Type>(s);
  else if(x < 0) return Type(0);
  else return (s == 0) * Type(-0.5);

}

template <class Type>
class LSImap : public CutFEMmap<Type>, public LimLimap <Type> {
  protected:

    LSImap() : LimLimap<Type>() {

    };
    ~LSImap() {};

    void clear() {
      LimLimap<Type>::clear();

    };

    Type LSIm1(const int &m, const Type &a, Type d);
    Type LSI(const int &s, const unsigned &m, const Type &a, const Type &d);

    Type lsi(const int &s, const unsigned &m, const Type &a, const Type &d) {

      return this->LSI(s, m, a, d);

    }

    Type lsi(const int &s, const unsigned &m, const std::pair<Type, Type> &key) {

      return this->LSI(s, m, key.first, key.second);

    }

};

template <class Type>
Type LSImap<Type>::LSIm1(const int &m, const Type &a, Type d) {

  if(a == 0) {
    throw std::invalid_argument("A planar line delta requires a nonzero slope.");
  }

  d = -d / a;
  if(d > 0 && d < 1) {
    return pow(d, m) / fabs(a);
  }
  else if(d == 0 || d == 1) {
    return (m == 0) ? 0.5 / fabs(a) : 0.5 * pow(d, m) / fabs(a);
  }
  else {
    return Type(0);
  }
}

template <class Type>
Type LSImap<Type>::LSI(const int &s, const unsigned &m, const Type &a, const Type &d) {

  switch(s) {
    case -1:
      return this->LSIm1(m, a, d);
      break;

    default:

      if(a == 0) {
        return -this->limLi(s, d) / Type(m + 1);
      }

      Type INT(0);
      Type x(a + d);

      if(x <= 0 || d <= 0) {

        if(x > 0) {
          Type g =  1 / (-a);
          for(unsigned i = 1; i <= m + 1; i++) {
            INT += g * this->limLi(s + i, x);
            g *= (m + 1 - i) / (-a);
          }
        }

        else if(d > 0) {
          INT -= this->limLi(s + m + 1, d) * boost::math::factorial<Type>(m) / pow(-a, m + 1);
        }
      }
      else {

        for(unsigned i = 0; i < s; i++) {
          INT -= pow(-a, i) / boost::math::factorial<Type>(m + 1 + i) * this->limLi(s - i, x) ;
        }

        INT += pow(-a, s) / boost::math::factorial<Type>(m + 1 + s);
        INT *= boost::math::factorial<Type>(m);
      }
      return INT;
  }
}

template <class TypeIO, class TypeA>
class HCImap :  public LSImap <TypeA> {
  public:
    HCImap() : LSImap<TypeA>() {

      _fast = false;
    };

    ~HCImap() {
      clear();
    };

    void clear() {
      LSImap<TypeA>::clear();

    };

    TypeIO operator()(const int &s, const std::vector<unsigned> &m, const std::vector<TypeIO> &a, const TypeIO &d);

  protected:

    TypeA hcia1(const unsigned & n, const int &s, std::vector<unsigned> &m, const std::vector <TypeA> &a, const std::vector <TypeA> &ma, const TypeA & d, const TypeA & md) {

      return this->HyperCubeA1(n, s, m, a, ma, d, md);

    }

  private:
    TypeA hcia(const int &s, const std::vector<unsigned> &m, const std::vector<TypeA> &a, const TypeA &d) {

      return this->HyperCubeA(s, m, a, d);

    }

    TypeA HyperCubeA(const int &s, const std::vector<unsigned> &m, const std::vector <TypeA> &a, const TypeA &d);
    TypeA HyperCubeA1(const unsigned & n, const int &s, std::vector<unsigned> &m, const std::vector <TypeA> &a, const std::vector <TypeA> &ma, const TypeA & d, const TypeA & md);
    TypeA HyperCubeB(const unsigned &n, const int &s, std::vector<unsigned> &m, const std::vector <TypeA> &a, const std::vector <TypeA> &ma, const TypeA &d, const TypeA &md);
    TypeA HyperCubeC(const unsigned & n, const int &s, std::vector<unsigned>& m, const std::vector <TypeA> &a, const std::vector <TypeA> &ma, const TypeA & d, const TypeA & md);

    TypeA _sum;

    std::vector <unsigned> _m;
    std::vector <TypeA> _a;
    std::vector <TypeA> _ma;
    std::vector<TypeA> _a2;
    bool _fast;

};

template <class TypeIO, class TypeA>
TypeIO HCImap<TypeIO, TypeA>::operator()(const int &s, const std::vector<unsigned> &m, const std::vector <TypeIO> &a, const TypeIO &d) {
  _a2.resize(a.size());
  TypeIO d1 = d;
  for(unsigned i = 0; i < a.size(); i++) {
    d1 -= a[i];
    _a2[i] = static_cast<TypeA>(2 * a[i]);
  }
  TypeA d2 = static_cast<TypeA>(d1);
  return static_cast<TypeIO>(this->hcia(s, m, _a2, d2));
}

template <class TypeIO, class TypeA>
TypeA HCImap<TypeIO, TypeA>::HyperCubeA(const int &s, const std::vector<unsigned> &m, const std::vector <TypeA> &a, const TypeA &d) {

  _a = a;
  _m = m;

  TypeA HCI = 1;

  for(int i = _a.size() - 1; i >= 0; i--) {
    if(_a[i] == 0) {
      HCI /= TypeA(_m[i] + 1);
      _a.erase(_a.begin() + i);
      _m.erase(_m.begin() + i);
    }
  }

  if(_a.size() > 0) {

    for(unsigned i = 0; i < _a.size() - 1; i++) {
      for(unsigned j = i + 1; j < _a.size(); j++) {
        if(fabs(_a[i]) > fabs(_a[j])) {
          std::swap(_a[i], _a[j]);
          std::swap(_m[i], _m[j]);
        }
      }
    }

    _ma.resize(_a.size());
    for(unsigned i = 0; i < _a.size(); i++)  _ma[i] = -_a[i];

    return HCI * this->hcia1(_a.size() - 1, s, _m, _a, _ma, d, -d);
  }
  else {
    return -HCI * this->limLi(s, d);
  }
}

template <class TypeIO, class TypeA>
TypeA HCImap<TypeIO, TypeA>::HyperCubeA1(const unsigned & n, const int &s, std::vector<unsigned> &m, const std::vector <TypeA> &a, const std::vector <TypeA> &ma, const TypeA & d, const TypeA & md) {

  if(_fast) goto fast;

  if(n == 0)  {
    return this->lsi(s, m[0], a[0], d);
  }

  _sum = d;
  for(unsigned i = 0; i <= n; i++) _sum += a[i];

  switch(s) {
    case -1:
      if(_sum <= 0) {
        return this->HyperCubeB(n, -1, m, a, ma, d, md);
      }
      else {
      fast:
        _fast = false;
        return this->HyperCubeB(n, -1, m, ma, a, md, d);

      }
      break;
    default:
      if(_sum <= fabs(a[n])) {
        return this->HyperCubeB(n, s, m, a, ma, d, md);
      }
      else {
        return this->HyperCubeC(n, s, m, a, ma, d, md);
      }
  }
}

template <class TypeIO, class TypeA>
TypeA HCImap<TypeIO, TypeA>::HyperCubeB(const unsigned &n, const int &s, std::vector<unsigned> &m, const std::vector <TypeA> &a, const std::vector <TypeA> &ma, const TypeA &d, const TypeA &md) {

  TypeA dpa = d + a[n];
  TypeA mdpa = -dpa;
  unsigned mp1 = m[n] + 1;
  unsigned nm1 = n - 1;
  int spmp1 = s + mp1;

  TypeA aI = 1 / a[n];
  TypeA c = aI;

  TypeA HCIb = 0;
  for(unsigned j = 1; j <= m[n]; c *= -aI * (mp1 - j), j++) {

    HCIb += this->hcia1(nm1, s + j, m, a, ma, dpa, mdpa) * c;
  }
  HCIb += this->hcia1(nm1, spmp1, m, a, ma, dpa, mdpa) * c;
  HCIb -= this->hcia1(nm1, spmp1, m, a, ma, d, md) * c;
  return HCIb;

}

template <class TypeIO, class TypeA>
TypeA HCImap<TypeIO, TypeA>::HyperCubeC(const unsigned & n, const int &s, std::vector<unsigned>& m, const std::vector <TypeA> &a, const std::vector <TypeA> &ma, const TypeA & d, const TypeA & md) {

  TypeA dpa = d + a[n];
  TypeA mdpa = -dpa;
  unsigned mp1 = m[n] + 1;
  unsigned nm1 = n - 1;

  TypeA c = 1 / TypeA(mp1);
  TypeA HCIc = 0;
  for(unsigned i = 0; i < s; i++, c *= -a[n] / (mp1 + i)) {
    HCIc += this->hcia1(nm1, s - i, m, a, ma, dpa, mdpa) * c;
  }
  HCIc += this->hcia1(nm1, 0, m, a, ma, dpa, mdpa) * c;
  c *= -a[n];

  _fast = true;
  m[n] += (s + 1);
  HCIc += this->hcia1(n, -1, m, a, ma, d, md) * c ;
  m[n] -= (s + 1);

  return HCIc;

}

template <class TypeIO, class TypeA>
class TRImap :  public LSImap <TypeA> {
  public:

    TRImap() : LSImap<TypeA>() {

    };
    ~TRImap() {};

    void clear() {
      LSImap<TypeA>::clear();

    };

    TypeIO operator()(const int &s, const std::vector<unsigned> &m, const std::vector<TypeIO> &a, const TypeIO &d) {
      return static_cast<TypeIO>(this->tria(s, m, {static_cast<TypeA>(-a[0]), static_cast<TypeA>(a[1])}, static_cast<TypeA>(d + a[0])));
    }

  protected:
    TypeA tria(const int &s, const std::vector<unsigned> &m, const std::vector<TypeA> &a, const TypeA &d) {

      return TriangleA(s, m, a, d);

    }

  private:
    TypeA TriangleA(const int &s, const std::vector<unsigned> &m, const std::vector <TypeA> &a, const TypeA &d);
    TypeA TriangleBReduced(const int &s, const std::vector<unsigned> &m_input, const std::vector <TypeA> &a_input, const TypeA &d);
    TypeA TriangleBFull(const int &s, const std::vector<unsigned> &m, const std::vector <TypeA> &a, const TypeA &d);
    TypeA TriangleC(const int &s, const std::vector<unsigned> &m, const std::vector <TypeA> &a, const TypeA &d);

};

template <class TypeIO, class TypeA>
TypeA TRImap<TypeIO, TypeA>::TriangleA(const int &s, const std::vector<unsigned> &m, const std::vector <TypeA> &a, const TypeA &d) {

  if(a[0] == 0 && a[1] == 0) return -this->limLi(s, d) / TypeA((m[0] + m[1] + 2) * (m[1] + 1));

  TypeA sum = a[0] + a[1] + d;
  switch(s) {
    case -1:
      if(sum <= 0) return TriangleBReduced(-1, m, a, d);
      else return TriangleBReduced(-1, m, std::vector<TypeA> {-a[0], -a[1]}, -d);
      break;
    default:
      if(sum <= 0) {
        return TriangleBReduced(s, m, a, d);
      }
      else if(sum <= std::max(fabs(a[1]), fabs(a[0]))) {
        return TriangleBFull(s, m, a, d);
      }
      else {
        return TriangleC(s, m, a, d);
      }
  }
}

template <class TypeIO, class TypeA>
TypeA TRImap<TypeIO, TypeA>::TriangleBReduced(const int &s, const std::vector<unsigned> &m_input, const std::vector <TypeA> &a_input, const TypeA &d) {

  const TypeA &a = a_input[0];
  const TypeA &b = a_input[1];
  const unsigned &m = m_input[0];
  const unsigned &n = m_input[1];

  TypeA TRI = 0.;
  if(b == 0) {
    TRI = (-this->limLi(s + 1, a + d) + this->limLi(s + m + n + 2, d) * boost::math::factorial<TypeA>(m + n + 1) * pow(-1 / a, m + n + 1) ) / (a * (n + 1));
  }
  else if(a == 0) {
    TRI = (this->limLi(s + n + 1, d) * boost::math::factorial<TypeA>(n) * pow(-1. / b, n) - this->limLi(s + m + n + 2, d) * boost::math::factorial<TypeA>(m + n + 1) * pow(-1 / b, m + n + 1) ) / (b * (m + 1));
  }
  else if(a + b == 0) {
    if(a + d  >= 0.) {
      for(unsigned i = 1; i <= m + 1; i++) {
        TRI += this->limLi(s + n + 1 + i, a + d) * pow(-1 / a, i) / boost::math::factorial<TypeA>(m + 1 - i);
      }
      TRI *= boost::math::factorial<TypeA>(n) * boost::math::factorial<TypeA>(m) * pow(1 / a, n + 1);
    }
    TRI += this->limLi(s + 1, d) / (a * (m + n + 1));
  }
  else {
    if(fabs(b) > fabs(a)) {
      if(d > 0) {
        for(unsigned i = 1; i <= n + 1; i++) {
          TRI += boost::math::factorial<TypeA>(m + n + 1 - i) / boost::math::factorial<TypeA>(n + 1 - i) * pow((a + b) / b, i);
        }
        TRI *= this->limLi(s + m + n + 2, d) / pow(-a - b, m + n + 2) ;
      }
      TRI += this->lsi(s + n + 1, m, a, d) * pow(-1 / b, n + 1);
      TRI *= boost::math::factorial<TypeA>(n);
    }
    else {
      if(d > 0) {
        for(unsigned i = 1; i <= m + 1; i++) {
          TRI += boost::math::factorial<TypeA>(m + n + 1 - i) / boost::math::factorial<TypeA>(m + 1 - i) * pow((a + b) / a, i);
        }
        TRI *= - this->limLi(s + m + n + 2, d) / pow(-a - b, m + n + 2) ;
      }
      if(a + d > 0) {
        TypeA TRI2 = 0.;
        for(unsigned i = 1; i <= m + 1; i++) {
          TRI2 +=  this->limLi(s + n + i + 1, a + d) / (boost::math::factorial<TypeA>(m + 1 - i) * pow(-a, i));
        }
        TRI += TRI2 * boost::math::factorial<TypeA>(n) / pow(-b, n + 1) ;
      }
      TRI *= boost::math::factorial<TypeA>(m);
    }

  }
  return TRI;
}

template <class TypeIO, class TypeA>
TypeA TRImap<TypeIO, TypeA>::TriangleBFull(const int &s, const std::vector<unsigned> &m_input, const std::vector <TypeA> &a_input, const TypeA &d) {

  const TypeA &a = a_input[0];
  const TypeA &b = a_input[1];
  const unsigned &m = m_input[0];
  const unsigned &n = m_input[1];

  TypeA TRI = 0;
  if(b == 0) {
    if(a == 0) TRI = -this->limLi(s, d) / ((m + n + 2) * (n + 1));
    else TRI = this->lsi(s, m + n + 1, a, d) / (n + 1);
  }
  else if(a == 0) {
    TRI = (this->lsi(s, n, b, d) - this->lsi(s, m + n + 1, b, d)) / (m + 1);
  }
  else if(a + b == 0) {
    for(unsigned i = 1; i <= n + 1; i++) {
      TRI += this->limLi(s + i, d) * pow(1. / a, i) / (boost::math::factorial<TypeA>(n + 1 - i) * (m + n + 2 - i));
    }
    TRI += this->lsi(s + n + 1, m, a, d) / pow(a, n + 1);
    TRI *= boost::math::factorial<TypeA>(n);
  }
  else {
    if(fabs(a) < fabs(b)) {
      for(unsigned j = 1; j <= n + 1; j++) {
        TRI -= this->lsi(s + j, m + n + 1 - j, a + b, d) * pow(-1 / b, j) / boost::math::factorial<TypeA>(n + 1 - j);
      }
      TRI += this->lsi(s + n + 1, m, a, d) * pow(-1 / b, n + 1);
      TRI *= boost::math::factorial<TypeA>(n);
    }
    else {
      for(unsigned j = 1; j <= m + 1; j++) {
        TRI += (this->lsi(s + j, m + n + 1 - j, a + b, d) - this->lsi(s + j, n, b, a + d)) * pow(-1 / a, j) / boost::math::factorial<TypeA>(m + 1 - j);
      }
      TRI *= boost::math::factorial<TypeA>(m);
    }
  }

  return TRI;
}

template <class TypeIO, class TypeA>
TypeA TRImap<TypeIO, TypeA>::TriangleC(const int &s, const std::vector<unsigned> &m_input, const std::vector <TypeA> &a_input, const TypeA &d) {

  const TypeA &a = a_input[0];
  const TypeA &b = a_input[1];
  const unsigned &m = m_input[0];
  const unsigned &n = m_input[1];

  TypeA TRI = 0;
  if(true && fabs(b) > fabs(a)) {

    for(int i = s; i >= 0; i--) {
      TRI += this->lsi(s - i, m + n + 1 + i, a + b, d) * pow(-b, i) / boost::math::factorial<TypeA>(n + 1 + i);
    }
    TRI += this->TriangleBReduced(-1, {m, n + s + 1}, {-a, -b}, -d) * pow(-b, s + 1) / boost::math::factorial<TypeA>(n + s + 1);
    TRI *= boost::math::factorial<TypeA>(n);
  }
  else {

    for(int i = 0; i <= s; i++) {
      TRI += (this->lsi(s - i, n, b, a + d) - this->lsi(s - i, m + n + 1 + i, a + b, d)) * pow(-a, i) / boost::math::factorial<TypeA>(m + 1 + i);
    }
    TRI += TriangleBReduced(-1, {m + s + 1, n}, {-a, -b}, -d) * pow(-a, s + 1) / boost::math::factorial<TypeA>(m + s + 1);
    TRI *= boost::math::factorial<TypeA>(m);
  }
  return TRI;
}

template <class TypeIO, class TypeA>
class TTImap : public TRImap <TypeA, TypeA> {
  public:

    TTImap() : TRImap<TypeA, TypeA>() {

      _baseType = 0;

    };

    ~TTImap() {
      clear();
    };

    void clear() {
      TRImap<TypeA, TypeA>::clear();

    };

    TypeIO operator()(const int &s, const std::vector<unsigned> &m, const std::vector<TypeIO> &a, const TypeIO &d);

    void SetBaseType(const unsigned &value) {
      _baseType = value;
    };

  protected:
    TypeA ttia(const int &s, const std::vector<unsigned> &m, const std::vector<TypeA> &a, const TypeA &d) {

      return TetrahedronA(s, m, a, d);

    }

  private:
    TypeA TetrahedronA(const int &s, const std::vector<unsigned> &m, const std::vector <TypeA> &a, const TypeA &d);
    TypeA TetrahedronB(const int &s, const std::vector<unsigned> &m, const std::vector <TypeA> &a, const TypeA &d);
    TypeA TetrahedronC(const int &s, const std::vector<unsigned> &m, const std::vector <TypeA> &a, const TypeA &d);

    unsigned _baseType;

};

template <class TypeIO, class TypeA>
TypeIO TTImap<TypeIO, TypeA>::operator()(const int &s, const std::vector<unsigned> &m, const std::vector <TypeIO> &a, const TypeIO & d) {

  if(a[0] == a[1] && a[1] == a[2]) {
    return this->lsi(s, m[0] + m[1] + m[2] + 2, std::make_pair<TypeA, TypeA>(static_cast<TypeA>(a[0]), static_cast<TypeA>(d))) / ((m[2] + 1) * (m[1] + m[2] + 2));
  }

  TypeIO m1 = std::max(fabs(a[1] - a[0]), fabs(a[2] - a[1]));
  TypeIO m2 = fabs(a[0] - a[2]);

  if((_baseType == 0 && m1 >= m2) || _baseType == 1) {

    return static_cast<TypeIO>( this->ttia(s, {m[0], m[1], m[2]}, {static_cast<TypeA>(a[0]), static_cast<TypeA>(a[1] - a[0]), static_cast<TypeA>(a[2] - a[1])}, static_cast<TypeA>(d)));
  }

  else {

    return static_cast<TypeIO>( this->ttia(s, {m[0], m[1], m[2]}, {static_cast<TypeA>(a[1]), static_cast<TypeA>(a[2] - a[1]), static_cast<TypeA>(a[0] - a[2])}, static_cast<TypeA>(d)));
  }
}

template <class TypeIO, class TypeA>
TypeA TTImap<TypeIO, TypeA>::TetrahedronA(const int &s, const std::vector<unsigned> &m, const std::vector <TypeA> &a, const TypeA & d) {

  switch(s) {
    case -1:
      if(a[0] + a[1] + a[2] + d <= 0) return TetrahedronB(-1, m, a, d);
      else return TetrahedronB(-1, m, {-a[0], -a[1], -a[2]}, -d);
      break;
    default:
      if(a[0] + a[1] + a[2] + d <= std::max(fabs(a[1]), fabs(a[2]))) {
        return TetrahedronB(s, m, a, d);
      }
      else {
        return TetrahedronC(s, m, a, d);
      }
  }
}

template <class TypeIO, class TypeA>
TypeA TTImap<TypeIO, TypeA>::TetrahedronB(const int &s, const std::vector<unsigned> &m_input, const std::vector <TypeA> &a_input, const TypeA & d) {

  const TypeA &a = a_input[0];
  const TypeA &b = a_input[1];
  const TypeA &c = a_input[2];
  const unsigned &m = m_input[0];
  const unsigned &n = m_input[1];
  const unsigned &o = m_input[2];

  TypeA TET = 0;
  if(fabs(c) > fabs(b)) {
    for(unsigned i = 1; i <= o + 1; i++)  {
      TET -= this->tria(s + i, {m, n + o + 1u - i}, {a, b + c}, d) / (boost::math::factorial<TypeA> (o + 1u - i) * pow(-c, i));
    }
    TET += this->tria(s + o + 1, {m, n}, {a, b}, d) / pow(-c, o + 1);
    TET *= boost::math::factorial<TypeA>(o);
  }
  else {
    for(unsigned i = 1; i <= n + 1; i++)  {
      TET += (-this->tria(s + i, {m + n + 1u - i, o}, {a + b, c}, d) + this->tria(s + i, {m, n + o + 1u - i}, {a, b + c}, d)) / (boost::math::factorial<TypeA> (n + 1u - i) * pow(-b, i));
    }
    TET *= boost::math::factorial<TypeA>(n);
  }
  return TET;
}

template <class TypeIO, class TypeA>
TypeA TTImap<TypeIO, TypeA>::TetrahedronC(const int &s, const std::vector<unsigned> &m_input, const std::vector <TypeA> &a_input, const TypeA & d) {

  const TypeA &a = a_input[0];
  const TypeA &b = a_input[1];
  const TypeA &c = a_input[2];
  const unsigned &m = m_input[0];
  const unsigned &n = m_input[1];
  const unsigned &o = m_input[2];

  TypeA TET = 0;
  if(fabs(c) > fabs(b)) {
    for(unsigned i = 0; i <= s; i++)  {
      TET += this->tria(s - i, {m, n + o + i + 1u}, {a, b + c}, d) * pow(-c, i) / boost::math::factorial<TypeA>(o + i + 1);
    }
    TET += this->ttia(-1, {m, n, o + s + 1u}, a_input, d) * pow(-c, s + 1u) / boost::math::factorial<TypeA>(o + s + 1);

    TET *= boost::math::factorial<TypeA> (o);
  }
  else {
    for(unsigned i = 0; i <= s; i++)  {
      TET += (this->tria(s - i, {m + n + i + 1u, o}, {a + b, c}, d) - this->tria(s - i, {m, n + o + i + 1u}, {a, b + c}, d)) * pow(-b, i) / boost::math::factorial<TypeA>(n + i + 1);
    }
    TET += this->ttia(-1, {m, n + s + 1u, o}, a_input, d) * pow(-b, s + 1u) / boost::math::factorial<TypeA>(n + s + 1);

    TET *= boost::math::factorial<TypeA> (n);
  }

  return TET;
}

}

#endif
