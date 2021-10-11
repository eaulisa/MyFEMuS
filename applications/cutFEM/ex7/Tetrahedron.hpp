#ifndef __femus_cut_fem_TET_hpp__
#define __femus_cut_fem_TET_hpp__

#include "Triangle.hpp"

template <class TypeIO, class TypeA>
class TTImap : public TRImap <TypeA, TypeA> {
  public:

    TTImap(const unsigned &mMax, const unsigned &sMax = 0, const unsigned &ds = 0) : TRImap <TypeA, TypeA> (mMax + 1, sMax, mMax + 1) {

      _TTImap.resize(2u + sMax + ds);
      unsigned max = 2u + mMax + sMax;

      for(unsigned s = 0; s < _TTImap.size(); s++) {
        _TTImap[s].resize(max);
        for(unsigned i = 0; i < _TTImap[s].size(); i++) {
          _TTImap[s][i].resize(max - i);
          for(unsigned j = 0; j < _TTImap[s][i].size(); j++) {
            _TTImap[s][i][j].resize(max - i - j);
          }
        }
      }
    };

    ~TTImap() {};

    void clear() {
      TRImap<TypeA, TypeA>::clear();
      for(unsigned s = 0; s < _TTImap.size(); s++) {
        for(unsigned i = 0; i < _TTImap[s].size(); i++) {
          for(unsigned j = 0; j < _TTImap[s][i].size(); j++) {
            for(unsigned k = 0; k < _TTImap[s][i][j].size(); k++) {
              _TTImap[s][i][j][k].clear();
            }
          }
        }
      }
      cnt = 0;
    };

    void printCounter() {
      TRImap<TypeIO, TypeA>::printCounter();
      std::cout << "TTI counter = " << cnt << std::endl;
    }

    TypeA TetrahedronA(const int &s, const std::vector<unsigned> &m, const std::vector <TypeA> &a, const TypeA &d);
    TypeA TetrahedronB(const int &s, const std::vector<unsigned> &m, const std::vector <TypeA> &a, const TypeA &d);
    TypeA TetrahedronC(const int &s, const std::vector<unsigned> &m, const std::vector <TypeA> &a, const TypeA &d);
    TypeIO Tetrahedron(const int &s, const std::vector<unsigned> &m, const std::vector <TypeIO> &a, const TypeIO &d);

    TypeA ttia(const int &s, const std::vector<unsigned> &m, const std::vector<TypeA> &a, const TypeA &d) {
      _index = std::make_pair(a, d);
      _it = _TTImap[s + 1][m[0]][m[1]][m[2]].find(_index);
      if(_it == _TTImap[s + 1][m[0]][m[1]][m[2]].end()) {
        cnt++;  
        //std::cout << "n3 s = " << s << " m = " << m[0] << ", " << m[1] << ", " << m[2] << " a = " << a[0] << ", " << a[1] << ", " << a[1] << " d = " << d << std::endl;;
        _TTImap[s + 1][m[0]][m[1]][m[2]][_index] = TetrahedronA(s, m, a, d);
      }
      return _TTImap[s + 1][m[0]][m[1]][m[2]][_index];
    }

    TypeIO operator()(const int &s, const std::vector<unsigned> &m, const std::vector<TypeIO> &a, const TypeIO &d) {
      return this->Tetrahedron(s, m, a, d);
    }

  private:
    std::vector<std::vector<std::vector<std::vector<std::map < std::pair<std::vector<TypeA>, TypeA>, TypeA > > > > >_TTImap;
    typename std::map < std::pair<std::vector<TypeA>, TypeA>, TypeA >::iterator _it;
    std::pair<std::vector<TypeA>, TypeA> _index;

    static unsigned cnt;
};

template <class TypeIO, class TypeA> unsigned TTImap <TypeIO, TypeA>::cnt = 0;

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
      TET -= this->tria(s + i, {m, n + o + 1u - i}, {a, b + c}, d) / (factorial<TypeA> (o + 1u - i) * pow(-c, i));
    }
    TET += this->tria(s + o + 1, {m, n}, {a, b}, d) / pow(-c, o + 1);
    TET *= factorial<TypeA>(o);
  }
  else {
    for(unsigned i = 1; i <= n + 1; i++)  {
      TET += (-this->tria(s + i, {m + n + 1u - i, o}, {a + b, c}, d) +
              this->tria(s + i, {m, n + o + 1u - i}, {a, b + c}, d)) / (factorial<TypeA> (n + 1u - i) * pow(-b, i));
    }
    TET *= factorial<TypeA>(n);
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
      TET += this->tria(s - i, {m, n + o + i + 1u}, {a, b + c}, d) * pow(-c, i) / factorial<TypeA>(o + i + 1);
    }
    TET += this->ttia(-1, {m, n, o + s + 1u}, a_input, d) * pow(-c, s + 1u) / factorial<TypeA>(o + s + 1);

    TET *= factorial <TypeA> (o);
  }
  else {
    for(unsigned i = 0; i <= s; i++)  {
      TET += (this->tria(s - i, {m + n + i + 1u, o}, {a + b, c}, d)
              - this->tria(s - i, {m, n + o + i + 1u}, {a, b + c}, d)) * pow(-b, i) / factorial<TypeA>(n + i + 1);
    }
    TET += this->ttia(-1, {m, n + s + 1u, o}, a_input, d) * pow(-b, s + 1u) / factorial<TypeA>(n + s + 1);

    TET *= factorial <TypeA> (n);
  }


  return TET;
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
TypeIO TTImap<TypeIO, TypeA>::Tetrahedron(const int &s, const std::vector<unsigned> &m, const std::vector <TypeIO> &a, const TypeIO & d) {

  TypeIO m1 = std::max(fabs(a[1] + a[0]), fabs(a[2] - a[1]));
  TypeIO m2 = std::max(fabs(a[2] + a[1]), fabs(a[0] - a[2]));
  TypeIO m3 = std::max(fabs(a[0] + a[2]), fabs(a[1] - a[0]));

  if(m1 > m2 && m1 > m3) {
    return static_cast<TypeIO>(
    this->ttia(s, {m[0], m[1], m[2]},
    {static_cast<TypeA>(a[0]), static_cast<TypeA>(a[1] + a[0]), static_cast<TypeA>(a[2] - a[1])},
    static_cast<TypeA>(d)));
  }
  else if(m2 > m3) {
    return static_cast<TypeIO>(
    this->ttia(s, {m[1], m[2], m[0]},
    {static_cast<TypeA>(a[1]), static_cast<TypeA>(a[2] + a[1]), static_cast<TypeA>(a[0] - a[2])},
    static_cast<TypeA>(d)));
  }
  else {
    return static_cast<TypeIO>(
    this->ttia(s, {m[2], m[0], m[1]},
    {static_cast<TypeA>(a[2]), static_cast<TypeA>(a[0] + a[2]), static_cast<TypeA>(a[1] - a[0])},
    static_cast<TypeA>(d)));
  }

}





#endif
