
#ifndef __femus_cut_fem_SQUARE_hpp__
#define __femus_cut_fem_SQUARE_hpp__

#include "Line.hpp"
#include "CutFem.hpp"

template <class TypeIO, class TypeA>
class SQImap :/* public CutFEMmap<TypeIO, TypeA>,*/ public LSImap <TypeA> {
  public:
    SQImap(const unsigned &mMax, const unsigned &sMax = 0, const unsigned &ds = 0) : LSImap <TypeA> (mMax + 1, sMax, mMax + 1) {

      //  std::cout << mMax << " " << sMax << "\n";

      _SQImap.resize(2u + sMax + ds);
      _SQI2.resize(2u + sMax + ds);
      _SQI2init.resize(2u + sMax + ds);
      unsigned max = 2u + mMax + sMax;

      for(unsigned s = 0; s < _SQImap.size(); s++) {
        _SQImap[s].resize(max);
        _SQI2[s].resize(max);
        _SQI2init[s].resize(max);
        for(unsigned i = 0; i < _SQImap[s].size(); i++) {
          _SQImap[s][i].resize(max - i);
          _SQI2[s][i].resize(max - i);
          _SQI2init[s][i].resize(max - i);
          for(unsigned j = 0; j < _SQImap[s][i].size(); j++) {
            _SQI2[s][i][j].resize(4);
            _SQI2init[s][i][j].assign(4, false);
          }
        }
      }
      _cnt = 0;
      _cntB = 0;
    };

    ~SQImap() {
      clear();
    };

    void printCounter() {
      LSImap<TypeA>::printCounter();
      std::cout << "SQI counters = " << _cnt << " " << _cntB << std::endl;
//       for(int s = 0; s < _SQImap.size(); s++) {
//         std::cout << "s = " << -1 + s << std::endl;
//         for(unsigned i = 0; i < _SQImap[s].size(); i++) {
//           for(unsigned j = 0; j < _SQImap[s][i].size(); j++) {
//             std::cout << _SQImap[s][i][j].size() << " ";
//           }
//           std::cout << std::endl;
//         }
//       }
    }

    void clear() {
      LSImap<TypeA>::clear();
      for(unsigned s = 0; s < _SQImap.size(); s++) {
        for(unsigned i = 0; i < _SQImap[s].size(); i++) {
          for(unsigned j = 0; j < _SQImap[s][i].size(); j++) {
            _SQImap[s][i][j].clear();
            _SQI2init[s][i][j].assign(4, false);
          }
        }
      }
      _cnt = 0;
      _cntB = 0;
    };

    //Call for the square [-1,1] x [-1,1]
    TypeIO operator()(const int &s, const std::vector<unsigned> &m, const std::vector<TypeIO> &a, const TypeIO &d);

  protected:

    TypeA sqiA(const int &s, const std::vector<unsigned> &m,
               const std::vector <TypeA> &a, const TypeA & d) {

      return this->SquareA(s, m, a, d);

      _keyA = std::make_tuple(a[0], a[1], d);
      _itA = _SQImap[s + 1][m[0]][m[1]].find(_keyA);

      if(_itA == _SQImap[s + 1][m[0]][m[1]].end()) {
        _cnt++;
        _IA = this->SquareA(s, m, a, d);
        _SQImap[s + 1][m[0]][m[1]][_keyA] = _IA;
        return _IA;
      }
      else {
        _cntB++;
        return _itA->second;
      }
    }




    TypeA sqiA1(const int &s, std::vector<unsigned> &m,
                const std::vector <TypeA> &a, const std::vector <TypeA> &ma,
                const TypeA & d, const TypeA & md) {

      //return this->SquareA1(s, m, a, ma, d, md);
      _keyA1 = std::make_tuple(a[0], a[1], d);
      _itA1 = _SQImap[s + 1][m[0]][m[1]].find(_keyA1);

      if(_itA1 == _SQImap[s + 1][m[0]][m[1]].end()) {
        _cnt++;
        _IA = this->SquareA1(s, m, a, ma, d, md);
        _SQImap[s + 1][m[0]][m[1]][_keyA1] = _IA;
        return _IA;
      }
      else {
        _cntB++;
        return _itA1->second;
      }
    }


    TypeA sqiA1(const int &s, std::vector<unsigned> &m, const std::tuple<TypeA, TypeA, TypeA> &keyA1,
                const std::vector <TypeA> &a, const std::vector <TypeA> &ma,
                const TypeA & d, const TypeA & md) {

      _itA1 = _SQImap[s + 1][m[0]][m[1]].find(keyA1);

      if(_itA1 == _SQImap[s + 1][m[0]][m[1]].end()) {
        _cnt++;
        _IA = this->SquareA1(s, m, a, ma, d, md);
        _SQImap[s + 1][m[0]][m[1]][keyA1] = _IA;
        return _IA;
      }
      else {
        _cntB++;
        return _itA1->second;
      }
    }

    TypeA sqiA2(const int &s, std::vector<unsigned> &m, const unsigned &idx,
                const std::vector <TypeA> &a, const std::vector <TypeA> &ma,
                const TypeA & d, const TypeA & md) {

      if(!_SQI2init[s + 1][m[0]][m[1]][idx]) {
        _SQI2[s + 1][m[0]][m[1]][idx] = this->SquareA1(s, m, a, ma, d, md);
        _SQI2init[s + 1][m[0]][m[1]][idx] = true;
      }
      return _SQI2[s + 1][m[0]][m[1]][idx];
    }



  private:
    TypeA SquareA(const int &s, const std::vector<unsigned> &m, const std::vector <TypeA> &a, const TypeA &d);
    TypeA SquareA1(const int &s, std::vector<unsigned> &m,
                   const std::vector <TypeA> &a, const std::vector <TypeA> &ma,
                   const TypeA & d, const TypeA & md);
    TypeA SquareB(const int &s, const std::vector<unsigned> &m, const std::vector <TypeA> &a, const TypeA &d);
    TypeA SquareC(const int &s, std::vector<unsigned>& m,
                  const std::vector <TypeA> &a, const std::vector <TypeA> &ma,
                  const TypeA & d, const TypeA & md);

    typedef std::tuple< TypeA, TypeA, TypeA> keydef;
    std::vector<std::vector<std::vector< std::map < keydef, TypeA > > > > _SQImap;

    std::vector<std::vector<std::vector< std::vector<TypeA > > > > _SQI2;
    std::vector<std::vector<std::vector< std::vector<bool > > > > _SQI2init;

    typename std::map < keydef, TypeA >::iterator _itA, _itA1;
    keydef _keyA, _keyA1;

    TypeA _IA, _sum;
    unsigned _cnt, _cntB;
    std::vector <unsigned> _m;
    std::vector <TypeA> _a;
    std::vector <TypeA> _ma;
    std::vector<TypeA> _a2;

};



template <class TypeIO, class TypeA>
TypeIO SQImap<TypeIO, TypeA>::operator()(const int &s, const std::vector<unsigned> &m, const std::vector <TypeIO> &a, const TypeIO &d) {
  _a2.resize(a.size());
  TypeIO d1 = d;
  for(unsigned i = 0; i < a.size(); i++) {
    d1 -= a[i];
    _a2[i] = static_cast<TypeA>(2 * a[i]);
  }
  TypeA d2 = static_cast<TypeA>(d1);
  return pow(2, _a2.size()) * static_cast<TypeIO>(this->sqiA(s, m, _a2, d2));
}

template <class TypeIO, class TypeA>
TypeA SQImap<TypeIO, TypeA>::SquareA(const int &s, const std::vector<unsigned> &m, const std::vector <TypeA> &a, const TypeA &d) {

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

  if(_a.size() > 1) { // all the left a[i] coefficients \ne 0
    _ma.resize(_a.size());
    for(unsigned i = 0; i < _a.size(); i++) { // reorder iterated integral from the smallest to the largest coefficient a[i]
      for(unsigned j = i + 1; j < _a.size(); j++) {
        if(fabs(_a[i]) > fabs(_a[j])) {
          std::swap(_a[i], _a[j]);
          std::swap(_m[i], _m[j]);
        }
      }
      _ma[i] = -_a[i];
    }
    return HCI * this->SquareA1(s, _m, _a, _ma, d, -d);
  }
  else if(_a.size() > 0) {
    return HCI * this->lsi(s, _m[0], std::make_pair(_a[0], d));
  }
  else {
    return -HCI * this->limLi(s, d);
  }
}

template <class TypeIO, class TypeA>
TypeA SQImap<TypeIO, TypeA>::SquareA1(const int &s, std::vector<unsigned> &m,
                                      const std::vector <TypeA> &a, const std::vector <TypeA> &ma,
                                      const TypeA & d, const TypeA & md) {

  _sum = a[0] + a[1] + d;

  switch(s) {
    case -1: // interface integral
      if(_sum <= 0) {
        return SquareB(-1, m, a, d);
      }
      else {
        return SquareB(-1, m, ma, md);
      }
      break;
    default:
      if(_sum <= fabs(a[1])) {
        return SquareB(s, m, a, d);
      }
      else {
        return SquareC(s, m, a, ma, d, md);
      }
  }
}

template <class TypeIO, class TypeA>
TypeA SQImap<TypeIO, TypeA>::SquareB(const int &s, const std::vector<unsigned> &m, const std::vector <TypeA> &a, const TypeA &d) {

  TypeA dpa = d + a[1];
  TypeA mdpa = -dpa;
  unsigned mp1 = m[1] + 1;
  int spmp1 = s + mp1;

  TypeA aI = 1 / a[1];
  TypeA c = aI;

  std::pair <TypeA, TypeA> lsikey(a[0], dpa);

  TypeA SQIb = 0;
  for(unsigned j = 1; j <= m[1]; c *= -aI * (mp1 - j), j++) {
    SQIb += this->lsi(s + j, m[0], lsikey) * c;
  }
  SQIb += this->lsi(spmp1, m[0], lsikey) * c;
  SQIb -= this->lsi(spmp1, m[0], std::make_pair(a[0], d)) * c;
  return SQIb;

}

template <class TypeIO, class TypeA>
TypeA SQImap<TypeIO, TypeA>::SquareC(const int &s, std::vector<unsigned>& m,
                                     const std::vector <TypeA> &a, const std::vector <TypeA> &ma,
                                     const TypeA & d, const TypeA & md) { //alternative formula at s = -1

  TypeA dpa = d + a[1];
  TypeA mdpa = -dpa;
  unsigned mp1 = m[1] + 1;

  std::pair <TypeA, TypeA> lsikey(a[0], dpa);

  TypeA c = 1 / TypeA(mp1);
  TypeA SQIc = 0;
  for(unsigned i = 0; i < s; i++, c *= -a[1] / (mp1 + i)) {
    SQIc += this->lsi(s - i, m[0], lsikey) * c;
  }
  SQIc += this->lsi(0, m[0], lsikey) * c;
  c *= -a[1];

  m[1] += (s + 1);
  //SQIc += this->sqiA1(-1, m, a, ma, d, md) * c ;
  SQIc += this->SquareA1(-1, m, a, ma, d, md) * c ;
  m[1] -= (s + 1);

  return SQIc;

}

#endif







