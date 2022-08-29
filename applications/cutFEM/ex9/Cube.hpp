
#ifndef __femus_cut_fem_CBI_hpp__
#define __femus_cut_fem_CBI_hpp__

#include "Square.hpp"

template <class TypeIO, class TypeA>
class CBImap : public SQImap <TypeA, TypeA> {
  public:
    CBImap(const unsigned &mMax, const unsigned &sMax = 0, const unsigned &ds = 0) : SQImap <TypeA, TypeA> (mMax + 1, sMax, mMax + 1) {

      // std::cout << mMax << " " << sMax << "\n";

      unsigned dim = 3;
      _CBImapA.resize(dim);
      _CBImapA1.resize(dim);

      _cnt.assign(dim, 0);

      for(unsigned d = 0; d < _CBImapA.size(); d++) {
        _CBImapA[d].resize(2u + sMax + (d + 1 != dim) * (mMax + dim - d - 1));
        _CBImapA1[d].resize(2u + sMax + (d + 1 != dim) * (mMax + dim - d - 1));
      }

      _fast = false;
    };

    ~CBImap() {
      clear();
    };

    void printCounter() {
      SQImap <TypeA, TypeA>::printCounter();
      for(unsigned d = 0; d < _cnt.size(); d++) {
        std::cout << "CBI_" << d + 1 << " counters = " << _cnt[d] << std::endl;
      }
    }

    void clear() {
      SQImap <TypeA, TypeA>::clear();
      for(unsigned d = 0; d < _CBImapA.size(); d++) {
        for(unsigned s = 0; s < _CBImapA[d].size(); s++) {
          _CBImapA[d][s].clear();
          _CBImapA1[d][s].clear();
        }
        _cnt[d] = 0;
      }
    };

    TypeIO operator()(const int &s, const std::vector<unsigned> &m, const std::vector<TypeIO> &a, const TypeIO &d);

  protected:

    TypeA cbia1(const int &s, std::vector<unsigned> &m, const TypeA &a0,
                const std::vector <TypeA> &a, const std::vector <TypeA> &ma,
                const TypeA & d, const TypeA & md) {

      return this->CubeA1(s, m, a0, a, ma, d, md);

//       keydef  key = std::make_tuple(std::vector<unsigned>(m.begin(), m.begin() + n + 1),
//                                     std::vector<TypeA>(a.begin(), a.begin() + n + 1), d);
//       typename std::map < keydef, TypeA >::iterator it = _CBImapA1[n][s + 1].find(key);
//
//       if(it == _CBImapA1[n][s + 1].end()) {
//         _cnt[n]++;
//         _IA = this->HyperCubeA1(n, s, m, a, ma, d, md);
//         _CBImapA1[n][s + 1][key] = _IA;
//         return _IA;
//       }
//       else {
//         return it->second;
//       }
    }

  private:
    TypeA cbia(const int &s, const std::vector<unsigned> &m, const std::vector<TypeA> &a, const TypeA &d) {

      return this->CubeA(s, m, a, d);

//       _key = std::make_tuple(m, a, d);
//       _it = _CBImapA[a.size() - 1][s + 1].find(_key);
//       if(_it == _CBImapA[a.size() - 1][s + 1].end()) {
//         _IA = this->HyperCubeA(s, m, a, d);
//         _CBImapA[a.size() - 1][s + 1][_key] = _IA;
//         return _IA;
//       }
//       else {
//         return _it->second;
//       }
    }

    TypeA CubeA(const int &s, const std::vector<unsigned> &m, const std::vector <TypeA> &a, const TypeA &d);
    TypeA CubeA1(const int &s, std::vector<unsigned> &m, const TypeA &a0,
                 const std::vector <TypeA> &a, const std::vector <TypeA> &ma,
                 const TypeA & d, const TypeA & md);
    TypeA CubeB(const int &s, std::vector<unsigned> &m, const TypeA &a0,
                const std::vector <TypeA> &a, const std::vector <TypeA> &ma,
                const TypeA &d, const TypeA &md);
    TypeA CubeC(const int &s, std::vector<unsigned>& m, const TypeA &a0,
                const std::vector <TypeA> &a, const std::vector <TypeA> &ma,
                const TypeA & d, const TypeA & md);

    typedef std::tuple< std::vector <unsigned>, std::vector<TypeA>, TypeA> keydef;
    typedef std::tuple < TypeA, TypeA, TypeA> sqrKeydef;


    std::vector<std::vector< std::map < keydef, TypeA > > > _CBImapA;
    std::vector<std::vector< std::map < keydef, TypeA > > > _CBImapA1;

    typename std::map < keydef, TypeA >::iterator _it;
    keydef _key;
    sqrKeydef _sqrKey;


    TypeA _IA, _sum;
    std::vector <unsigned> _cnt;
    std::vector <unsigned> _m;
    std::vector <TypeA> _a;
    std::vector <TypeA> _ma;
    std::vector<TypeA> _a2;
    bool _fast;

};



template <class TypeIO, class TypeA>
TypeIO CBImap<TypeIO, TypeA>::operator()(const int &s, const std::vector<unsigned> &m, const std::vector <TypeIO> &a, const TypeIO &d) {
  _a2.resize(a.size());
  TypeIO d1 = d;
  for(unsigned i = 0; i < a.size(); i++) {
    d1 -= a[i];
    _a2[i] = static_cast<TypeA>(2 * a[i]);
  }
  TypeA d2 = static_cast<TypeA>(d1);
  return pow(2, _a2.size()) * static_cast<TypeIO>(this->cbia(s, m, _a2, d2));
}

template <class TypeIO, class TypeA>
TypeA CBImap<TypeIO, TypeA>::CubeA(const int &s, const std::vector<unsigned> &m, const std::vector <TypeA> &a, const TypeA &d) {

  _a = a;
  _m = m;

  TypeA CBI = 1;

  for(int i = _a.size() - 1; i >= 0; i--) {
    if(_a[i] == 0) {
      CBI /= TypeA(_m[i] + 1);
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
    if(_a.size() == 3) return CBI * this->cbia1(s, _m, _a[0], _a, _ma, d, -d);
    else return CBI * this->sqiA1(s, _m, _a, _ma, d, -d);
  }
  else if(_a.size() > 0) {
    return CBI * this->lsi(s, _m[0], _a[0], d);
  }
  else {
    return -CBI * this->limLi(s, d);
  }
}

template <class TypeIO, class TypeA>
TypeA CBImap<TypeIO, TypeA>::CubeA1(const int &s, std::vector<unsigned> &m, const TypeA &a0,
                                    const std::vector <TypeA> &a, const std::vector <TypeA> &ma,
                                    const TypeA & d, const TypeA & md) {

  _sum = a[0] + a[1] + a[2] + d;

  switch(s) {
    case -1: // interface integral
      if(_sum <= 0) {
        return this->CubeB(-1, m, a0, a, ma, d, md);
      }
      else {
        return this->CubeB(-1, m, a0, ma, a, md, d);
      }
      break;
    default:
      if(_sum <= fabs(a[2])) {
        return this->CubeB(s, m, a0, a, ma, d, md);
      }
      else {
        return this->CubeC(s, m, a0, a, ma, d, md);
      }
  }
}

template <class TypeIO, class TypeA>
TypeA CBImap<TypeIO, TypeA>::CubeB(const int &s, std::vector<unsigned> &m, const TypeA &a0,
                                   const std::vector <TypeA> &a, const std::vector <TypeA> &ma,
                                   const TypeA &d, const TypeA &md) {

  unsigned shift = (a[0] == a0) ? 0 : 2;  
    
  TypeA dpa = d + a[2];
  TypeA mdpa = -dpa;
  unsigned mp1 = m[2] + 1;
  int spmp1 = s + mp1;

  TypeA aI = 1 / a[2];
  TypeA c = aI; // this is m!/(m+1-j)! 1/a^j for j = 1,...,m + 1

  TypeA CBIb = 0;

  //_sqrKey = std::make_tuple(a[0], a[1], dpa);
  for(unsigned j = 1; j <= m[2]; c *= -aI * (mp1 - j), j++) {
    //CBIb += this->sqiA1(s + j, m, _sqrKey, a, ma, dpa, mdpa) * c;  
    CBIb += this->sqiA2(s + j, m, shift + 1, a, ma, dpa, mdpa) * c;
  }
  //CBIb += this->sqiA1(spmp1, m, _sqrKey, a, ma, dpa, mdpa) * c;
  CBIb += this->sqiA2(spmp1, m, shift + 1, a, ma, dpa, mdpa) * c;
  
  //CBIb -= this->sqiA1(spmp1, m, a, ma, d, md) * c;
  CBIb -= this->sqiA2(spmp1, m, shift, a, ma, d, md) * c;
  return CBIb;

}

template <class TypeIO, class TypeA>
TypeA CBImap<TypeIO, TypeA>::CubeC(const int &s, std::vector<unsigned>& m, const TypeA &a0,
                                   const std::vector <TypeA> &a, const std::vector <TypeA> &ma,
                                   const TypeA & d, const TypeA & md) { //alternative formula at s = -1

 
  TypeA dpa = d + a[2];
  TypeA mdpa = -dpa;
  unsigned mp1 = m[2] + 1;

  TypeA c = 1 / TypeA(mp1);
  TypeA CBIc = 0;
  //_sqrKey = std::make_tuple(a[0], a[1], dpa);
  for(unsigned i = 0; i < s; i++, c *= -a[2] / (mp1 + i)) {
    //CBIc += this->sqiA1(s - i, m, _sqrKey, a, ma, dpa, mdpa) * c;
    CBIc += this->sqiA2(s - i, m, 1, a, ma, dpa, mdpa) * c;
  }
  //CBIc += this->sqiA1(0, m, _sqrKey, a, ma, dpa, mdpa) * c;
  CBIc += this->sqiA2(0, m, 1, a, ma, dpa, mdpa) * c;
  c *= -a[2];

  _fast = true;
  m[2] += (s + 1);
  CBIc += this->cbia1(-1, m, a0, a, ma, d, md) * c ;
  m[2] -= (s + 1);

  return CBIc;

}

#endif




