
#ifndef __femus_cut_fem_HCI_hpp__
#define __femus_cut_fem_HCI_hpp__

template<typename T>
std::vector<T> slice(const std::vector<T> &v, const unsigned &m, const unsigned &n) {
  std::vector<T> vec(n - m + 1);
  std::copy(v.begin() + m, v.begin() + n + 1, vec.begin());
  return vec;
}

template <class TypeIO, class TypeA>
class HCImap : public LSImap <TypeA> {
  public:
    HCImap(const unsigned &dim, const unsigned &mMax, const unsigned &sMax = 0) : LSImap <TypeA> (mMax + dim - 1, sMax, mMax + dim - 1) {

      _HCImapA.resize(dim);
      _HCImapA1.resize(dim);

      _cnt.assign(dim, 0);

      for(unsigned d = 0; d < _HCImapA.size(); d++) {
        _HCImapA[d].resize(2u + sMax + (d + 1 != dim) * (mMax + dim - d - 1));
        _HCImapA1[d].resize(2u + sMax + (d + 1 != dim) * (mMax + dim - d - 1));
      }
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

      typename std::map < std::tuple< std::vector <unsigned>, std::vector<TypeA>, TypeA>, TypeA >::iterator it;
      std::tuple< std::vector <unsigned>, std::vector<TypeA>, TypeA>  key;

      key = std::make_tuple(slice(m, 0, n), slice(a, 0, n), d);
      it = _HCImapA1[n][s + 1].find(key);
      if(it == _HCImapA1[n][s + 1].end()) {
        _cnt[n]++;
        _I1 = HyperCubeA1(n, s, m, a, ma, d, md);
        _HCImapA1[n][s + 1][key] = _I1;
        return _I1;
      }
      else {
        return it->second;
      }
    }

  private:

    TypeA hcia(const int &s, const std::vector<unsigned> &m, const std::vector<TypeA> &a, const TypeA &d) {
      _key = std::make_tuple(m, a, d);
      _it = _HCImapA[a.size() - 1][s + 1].find(_key);
      if(_it == _HCImapA[a.size() - 1][s + 1].end()) {
        _I1 = HyperCubeA(s, m, a, d);
        _HCImapA[a.size() - 1][s + 1][_key] = _I1;
        return _I1;
      }
      else {
        return _it->second;
      }
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


    std::vector<std::vector< std::map < std::tuple< std::vector <unsigned>, std::vector<TypeA>, TypeA>, TypeA > > > _HCImapA;
    std::vector<std::vector< std::map < std::tuple< std::vector <unsigned>, std::vector<TypeA>, TypeA>, TypeA > > > _HCImapA1;
    typename std::map < std::tuple< std::vector <unsigned>, std::vector<TypeA>, TypeA>, TypeA >::iterator _it;
    std::tuple< std::vector <unsigned>, std::vector<TypeA>, TypeA>  _key;

    TypeA _I1;
    std::vector <unsigned> _cnt;
    std::vector <unsigned> _m;
    std::vector <TypeA> _a;
    std::vector <TypeA> _ma;
    std::vector<TypeA> _a2;
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
  return static_cast<TypeIO>(pow(TypeIO(2), _a2.size()) * this->hcia(s, m, _a2, d2));
}

template <class TypeIO, class TypeA>
TypeA HCImap<TypeIO, TypeA>::HyperCubeA(const int &s, const std::vector<unsigned> &m, const std::vector <TypeA> &a, const TypeA &d) {

  _a = a;
  _m = m;

  TypeA HCI = 1;

  for(int i = _a.size() - 1; i >= 0; i--) {
    if(_a[i] == 0) {
      HCI *= 1 / TypeA(_m[i] + 1);
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

    return HCI * this->hcia1(-1 + _a.size(), s, _m, _a, _ma, d, -d);
  }
  else {
    return -HCI * this->limLi(s, d);
  }
}

template <class TypeIO, class TypeA>
TypeA HCImap<TypeIO, TypeA>::HyperCubeA1(const unsigned & n, const int &s, std::vector<unsigned> &m,
                                         const std::vector <TypeA> &a, const std::vector <TypeA> &ma,
                                         const TypeA & d, const TypeA & md) {

  if(n == 0)  {
    return this->lsi(s, m[0], a[0], d);
  }

  TypeA sum = d;
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

template <class TypeIO, class TypeA>
TypeA HCImap<TypeIO, TypeA>::HyperCubeB(const unsigned &n, const int &s, std::vector<unsigned> &m,
                                        const std::vector <TypeA> &a, const std::vector <TypeA> &ma,
                                        const TypeA &d, const TypeA &md) {

  TypeA aI = 1 / a[n];
  TypeA c = aI; // this is m!/(m+1-j)! 1/a^j for j = 1,...,m + 1

  TypeA I1 = 0;
  for(int j = 1; j <= m[n] + 1;  c *= -aI * (m[n] + 1 - j), j++) {
    I1 += this->hcia1(n - 1, s + j, m, a, ma, d + a[n], -d - a[n]) * c;
  }
  TypeA I2 = this->hcia1(n - 1, s + m[n] + 1, m, a, ma, d, md) * factorial<TypeA>(m[n]) * pow(-aI, m[n] + 1);
  return I1 + I2;

}

template <class TypeIO, class TypeA>
TypeA HCImap<TypeIO, TypeA>::HyperCubeC(const unsigned & n, const int &s, std::vector<unsigned>& m,
                                        const std::vector <TypeA> &a, const std::vector <TypeA> &ma,
                                        const TypeA & d, const TypeA & md) { //alternative formula at s = -1

  TypeA an = a[n];
  unsigned mn = m[n];

  m[n] += (s + 1);
  TypeA I1 = this->hcia1(n, -1, m, a, ma, d, md) * pow(-an, s + 1) * factorial<TypeA>(mn) / factorial<TypeA>(mn + s + 1) ;
  m[n] -= (s + 1);

  TypeA I2 = 0;
  TypeA c1 = 1 / TypeA(mn + 1);
  for(unsigned i = 0; i <= s; c1 *= -an / (mn + 2 + i), i++) {
    I2 += this->hcia1(n - 1, s - i, m, a, ma, d + an, -(d + an)) * c1;
  }
  return I1 + I2;

}

#endif

