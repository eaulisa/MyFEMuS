
#ifndef __femus_cut_fem_HCI_hpp__
#define __femus_cut_fem_HCI_hpp__

#include "CutLine.hpp"
#include "CutFem.hpp"

template <class TypeIO, class TypeA>
class HCImap : /*public CutFEMmap<TypeIO, TypeA>,*/ public LSImap <TypeA> {
  public:
    HCImap(const unsigned &dim, const unsigned &mMax, const unsigned &sMax = 0) : LSImap <TypeA> (mMax + dim - 1, sMax, mMax + dim - 1) {

      _HCImapA.resize(dim);
      _HCImapA1.resize(dim);

      _cnt.assign(dim, 0);

      for(unsigned d = 0; d < _HCImapA.size(); d++) {
        _HCImapA[d].resize(2u + sMax + (d + 1 != dim) * (mMax + dim - d - 1));
        _HCImapA1[d].resize(2u + sMax + (d + 1 != dim) * (mMax + dim - d - 1));
      }
      
      _fast = false;
    };

    ~HCImap() {
      clear();
    };

    void printCounter() {
      LSImap<TypeA>::printCounter();
      for(unsigned d = 0; d < _cnt.size(); d++) {
        std::cout << "HCI_" << d + 1 << " counters = " << _cnt[d] << std::endl;
      }
    }

    void clear() {
      LSImap<TypeA>::clear();
      for(unsigned d = 0; d < _HCImapA.size(); d++) {
        for(unsigned s = 0; s < _HCImapA[d].size(); s++) {
          _HCImapA[d][s].clear();
          _HCImapA1[d][s].clear();
        }
        _cnt[d] = 0;
      }
    };

    TypeIO operator()(const int &s, const std::vector<unsigned> &m, const std::vector<TypeIO> &a, const TypeIO &d);

  protected:

    TypeA hcia1(const unsigned & n, const int &s, std::vector<unsigned> &m,
                const std::vector <TypeA> &a, const std::vector <TypeA> &ma,
                const TypeA & d, const TypeA & md) {
     // if(n>=1) std::cout << n<<" "<< s<<" " << m[0]<<" "<< m[1]<<" "<<m[2]<<std::endl;
      return this->HyperCubeA1(n, s, m, a, ma, d, md);

//       keydef  key = std::make_tuple(std::vector<unsigned>(m.begin(), m.begin() + n + 1),
//                                     std::vector<TypeA>(a.begin(), a.begin() + n + 1), d);
//       typename std::map < keydef, TypeA >::iterator it = _HCImapA1[n][s + 1].find(key);
// 
//       if(it == _HCImapA1[n][s + 1].end()) {
//         _cnt[n]++;
//         _IA = this->HyperCubeA1(n, s, m, a, ma, d, md);
//         _HCImapA1[n][s + 1][key] = _IA;
//         return _IA;
//       }
//       else {
//         return it->second;
//       }
    }

  private:
    TypeA hcia(const int &s, const std::vector<unsigned> &m, const std::vector<TypeA> &a, const TypeA &d) {
      
      return this->HyperCubeA(s, m, a, d);

//       _key = std::make_tuple(m, a, d);
//       _it = _HCImapA[a.size() - 1][s + 1].find(_key);
//       if(_it == _HCImapA[a.size() - 1][s + 1].end()) {
//         _IA = this->HyperCubeA(s, m, a, d);
//         _HCImapA[a.size() - 1][s + 1][_key] = _IA;
//         return _IA;
//       }
//       else {
//         return _it->second;
//       }
    }

    TypeA HyperCubeA(const int &s, const std::vector<unsigned> &m, const std::vector <TypeA> &a, const TypeA &d);
    TypeA HyperCubeA1(const unsigned & n, const int &s, std::vector<unsigned> &m,
                      const std::vector <TypeA> &a, const std::vector <TypeA> &ma,
                      const TypeA & d, const TypeA & md);
    TypeA HyperCubeB(const unsigned &n, const int &s, std::vector<unsigned> &m,
                     const std::vector <TypeA> &a, const std::vector <TypeA> &ma,
                     const TypeA &d, const TypeA &md);
    TypeA HyperCubeC(const unsigned & n, const int &s, std::vector<unsigned>& m,
                     const std::vector <TypeA> &a, const std::vector <TypeA> &ma,
                     const TypeA & d, const TypeA & md);

    typedef std::tuple< std::vector <unsigned>, std::vector<TypeA>, TypeA> keydef;

    std::vector<std::vector< std::map < keydef, TypeA > > > _HCImapA;
    std::vector<std::vector< std::map < keydef, TypeA > > > _HCImapA1;

    typename std::map < keydef, TypeA >::iterator _it;
    keydef _key;

    TypeA _IA,_sum;
    std::vector <unsigned> _cnt;
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
  return/* pow(2, _a2.size()) **/ static_cast<TypeIO>(this->hcia(s, m, _a2, d2));
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
    // all the left a[i] coefficients \ne 0
    for(unsigned i = 0; i < _a.size() - 1; i++) { // reorder iterated integral from the smallest to the largest coefficient a[i]
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
TypeA HCImap<TypeIO, TypeA>::HyperCubeA1(const unsigned & n, const int &s, std::vector<unsigned> &m,
                                         const std::vector <TypeA> &a, const std::vector <TypeA> &ma,
                                         const TypeA & d, const TypeA & md) {

  if(_fast) goto fast;

  if(n == 0)  {
    return this->lsi(s, m[0], a[0], d);
  }

  _sum = d;
  for(unsigned i = 0; i <= n; i++) _sum += a[i];

  switch(s) {
    case -1: // interface integral
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
TypeA HCImap<TypeIO, TypeA>::HyperCubeB(const unsigned &n, const int &s, std::vector<unsigned> &m,
                                        const std::vector <TypeA> &a, const std::vector <TypeA> &ma,
                                        const TypeA &d, const TypeA &md) {

  TypeA dpa = d + a[n];
  TypeA mdpa = -dpa;
  unsigned mp1 = m[n] + 1;
  unsigned nm1 = n - 1;
  int spmp1 = s + mp1;

  TypeA aI = 1 / a[n];
  TypeA c = aI; // this is m!/(m+1-j)! 1/a^j for j = 1,...,m + 1

  TypeA HCIb = 0;
  for(unsigned j = 1; j <= m[n]; c *= -aI * (mp1 - j), j++) {
    //if(n==2) std::cout << m[n] << " " << j<< " AA ";  
    HCIb += this->hcia1(nm1, s + j, m, a, ma, dpa, mdpa) * c;
  }
  HCIb += this->hcia1(nm1, spmp1, m, a, ma, dpa, mdpa) * c;
  HCIb -= this->hcia1(nm1, spmp1, m, a, ma, d, md) * c;
  return HCIb;

}

template <class TypeIO, class TypeA>
TypeA HCImap<TypeIO, TypeA>::HyperCubeC(const unsigned & n, const int &s, std::vector<unsigned>& m,
                                        const std::vector <TypeA> &a, const std::vector <TypeA> &ma,
                                        const TypeA & d, const TypeA & md) { //alternative formula at s = -1

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

#endif




