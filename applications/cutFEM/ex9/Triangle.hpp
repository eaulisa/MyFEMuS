
#ifndef __femus_cut_fem_TRI_hpp__
#define __femus_cut_fem_TRI_hpp__

#include "Line.hpp"
#include "CutFem.hpp"


template <class TypeIO, class TypeA>
class TRImap : /*public CutFEMmap<TypeIO, TypeA>,*/ public LSImap <TypeA> {
  public:

    TRImap(const unsigned &mMax, const unsigned &sMax = 0, const unsigned &ds = 0) : LSImap <TypeA> (mMax + 1, sMax, mMax + 1) {

      _TRImap.resize(2u + sMax + ds);
      unsigned max = 2u + mMax + sMax;

      for(unsigned s = 0; s < _TRImap.size(); s++) {
        _TRImap[s].resize(max);
        for(unsigned i = 0; i < _TRImap[s].size(); i++) {
          _TRImap[s][i].resize(max - i);
        }
      }
      _cnt = 0;

    };
    ~TRImap() {};

    void clear() {
      LSImap<TypeA>::clear();
      for(unsigned s = 0; s < _TRImap.size(); s++) {
        for(unsigned i = 0; i < _TRImap[s].size(); i++) {
          for(unsigned j = 0; j < _TRImap[s][i].size(); j++) {
            _TRImap[s][i][j].clear();
          }
        }
      }
      _cnt = 0;
    };

    void printCounter() {
      LSImap<TypeA>::printCounter();
      std::cout << "TRI counter = " << _cnt << std::endl;
    }

    TypeIO operator()(const int &s, const std::vector<unsigned> &m, const std::vector<TypeIO> &a, const TypeIO &d) {
      return static_cast<TypeIO>(this->tria(s, m, {static_cast<TypeA>(-a[0]), static_cast<TypeA>(a[1])}, static_cast<TypeA>(d + a[0])));
    }
    
protected:
    TypeA tria(const int &s, const std::vector<unsigned> &m, const std::vector<TypeA> &a, const TypeA &d) {
      _key = std::make_pair(a, d);
      _it = _TRImap[s + 1][m[0]][m[1]].find(_key);
      if(_it == _TRImap[s + 1][m[0]][m[1]].end()) {
        _cnt++;
        _I1 = TriangleA(s, m, a, d);
        _TRImap[s + 1][m[0]][m[1]][_key] = _I1;
        return _I1;
      }
      else {
        return _it->second;
      }
    }
   
  private:
    TypeA TriangleA(const int &s, const std::vector<unsigned> &m, const std::vector <TypeA> &a, const TypeA &d);
    TypeA TriangleBReduced(const int &s, const std::vector<unsigned> &m_input, const std::vector <TypeA> &a_input, const TypeA &d);
    TypeA TriangleBFull(const int &s, const std::vector<unsigned> &m, const std::vector <TypeA> &a, const TypeA &d);
    TypeA TriangleC(const int &s, const std::vector<unsigned> &m, const std::vector <TypeA> &a, const TypeA &d);  
      
    std::vector<std::vector<std::vector<std::map < std::pair<std::vector<TypeA>, TypeA>, TypeA > > > > _TRImap;
    typename std::map < std::pair<std::vector<TypeA>, TypeA>, TypeA >::iterator _it;
    std::pair<std::vector<TypeA>, TypeA> _key;

    TypeA _I1;
    unsigned _cnt;
};

template <class TypeIO, class TypeA>
TypeA TRImap<TypeIO, TypeA>::TriangleA(const int &s, const std::vector<unsigned> &m, const std::vector <TypeA> &a, const TypeA &d) {

  if(a[0] == 0 && a[1] == 0) return -this->limLi(s, d) / TypeA((m[0] + m[1] + 2) * (m[1] + 1)); // only from higher dimensions calls n = <0,0,c>

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
  if(b == 0) { //parallel to right edge
    TRI = (-this->limLi(s + 1, a + d) +
           this->limLi(s + m + n + 2, d) * factorial<TypeA>(m + n + 1) * pow(-1 / a, m + n + 1)
          ) / (a * (n + 1));
  }
  else if(a == 0) { //parallel to bottom edge
    TRI = (this->limLi(s + n + 1, d) * factorial<TypeA>(n) * pow(-1. / b, n) -
           this->limLi(s + m + n + 2, d) * factorial<TypeA>(m + n + 1) * pow(-1 / b, m + n + 1)
          ) / (b * (m + 1));
  }
  else if(a + b == 0) { // parallel to y=x edge//TODO
    if(a + d  >= 0.) {
      for(unsigned i = 1; i <= m + 1; i++) {
        TRI += this->limLi(s + n + 1 + i, a + d) * pow(-1 / a, i) / factorial<TypeA>(m + 1 - i);
      }
      TRI *= factorial<TypeA>(n) * factorial<TypeA>(m) * pow(1 / a, n + 1);
    }
    TRI += this->limLi(s + 1, d) / (a * (m + n + 1));
  }
  else { // general case
    if(fabs(b) > fabs(a)) {
      if(d > 0) {
        for(unsigned i = 1; i <= n + 1; i++) {
          TRI += factorial<TypeA>(m + n + 1 - i) / factorial<TypeA>(n + 1 - i) * pow((a + b) / b, i);
        }
        TRI *= this->limLi(s + m + n + 2, d) / pow(-a - b, m + n + 2) ;
      }
      TRI += this->lsi(s + n + 1, m, a, d) * pow(-1 / b, n + 1);
      TRI *= factorial<TypeA>(n);
    }
    else {
      if(d > 0) {
        for(unsigned i = 1; i <= m + 1; i++) {
          TRI += factorial<TypeA>(m + n + 1 - i) / factorial<TypeA>(m + 1 - i) * pow((a + b) / a, i);
        }
        TRI *= - this->limLi(s + m + n + 2, d) / pow(-a - b, m + n + 2) ;
      }
      if(a + d > 0) {
        TypeA TRI2 = 0.;
        for(unsigned i = 1; i <= m + 1; i++) {
          TRI2 +=  this->limLi(s + n + i + 1, a + d) / (factorial<TypeA>(m + 1 - i) * pow(-a, i));
        }
        TRI += TRI2 * factorial<TypeA>(n) / pow(-b, n + 1) ;
      }
      TRI *= factorial<TypeA>(m);
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
  if(b == 0) { //parallel to right edge
    if(a == 0) TRI = -this->limLi(s, d) / ((m + n + 2) * (n + 1));
    else TRI = this->lsi(s, m + n + 1, a, d) / (n + 1);
  }
  else if(a == 0) { //parallel to bottom edge
    TRI = (this->lsi(s, n, b, d) - this->lsi(s, m + n + 1, b, d)) / (m + 1);
  }
  else if(a + b == 0) { //parallel to y = x edge //TODO
    for(unsigned i = 1; i <= n + 1; i++) {
      TRI += this->limLi(s + i, d) * pow(1. / a, i) / (factorial<TypeA>(n + 1 - i) * (m + n + 2 - i));
    }
    TRI += this->lsi(s + n + 1, m, a, d) / pow(a, n + 1);
    TRI *= factorial<TypeA>(n);
  }
  else { //generic case
    if(fabs(a) < fabs(b)) { // |a| < |b|
      for(unsigned j = 1; j <= n + 1; j++) {
        TRI -= this->lsi(s + j, m + n + 1 - j, a + b, d) * pow(-1 / b, j)
               / factorial<TypeA>(n + 1 - j);
      }
      TRI += this->lsi(s + n + 1, m, a, d) * pow(-1 / b, n + 1);
      TRI *= factorial<TypeA>(n);
    }
    else { // |b| < |a|
      for(unsigned j = 1; j <= m + 1; j++) {
        TRI += (this->lsi(s + j, m + n + 1 - j, a + b, d) -
                this->lsi(s + j, n, b, a + d)) * pow(-1 / a, j)
               / factorial<TypeA>(m + 1 - j);
      }
      TRI *= factorial<TypeA>(m);
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
    //std::cout<<"!";
    for(int i = s; i >= 0; i--) {
      TRI += this->lsi(s - i, m + n + 1 + i, a + b, d) * pow(-b, i) / factorial<TypeA>(n + 1 + i);
    }
    TRI += this->TriangleBReduced(-1, {m, n + s + 1}, {-a, -b}, -d) * pow(-b, s + 1) / factorial<TypeA>(n + s + 1);
    TRI *= factorial<TypeA>(n);
  }
  else {
    //std::cout<<"#";
    for(int i = 0; i <= s; i++) {
      TRI += (this->lsi(s - i, n, b, a + d) - this->lsi(s - i, m + n + 1 + i, a + b, d)) * pow(-a, i) / factorial<TypeA>(m + 1 + i);
    }
    TRI += TriangleBReduced(-1, {m + s + 1, n}, {-a, -b}, -d) * pow(-a, s + 1) / factorial<TypeA>(m + s + 1);
    TRI *= factorial<TypeA>(m);
  }
  return TRI;
}

#endif

