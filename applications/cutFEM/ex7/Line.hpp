#ifndef __femus_cut_fem_LSI_hpp__
#define __femus_cut_fem_LSI_hpp__

#include <iostream>
#include <boost/math/special_functions/factorials.hpp>

#include <map>

using boost::math::factorial;

template <class Type>
class LimLimap {
  public:

    LimLimap(const unsigned &sMax) {
      _LimLimap.resize(2u + sMax);
    };
    ~LimLimap() { _cnt = 0;};

    Type LimLi(const int &s, const Type &x);

    void clear() {
      for(unsigned s = 0; s < _LimLimap.size(); s++) {
        _LimLimap[s].clear();
        _cnt = 0;
      }
    };
    void printCounter() {
      std::cout << "limLi counter = " << _cnt << std::endl;
    }

    Type limLi(const int &s, const Type &d) {
      _it = _LimLimap[s + 1].find(d);
      if(_it == _LimLimap[s + 1].end()) {
        _cnt++;
        _LimLimap[s + 1][d] = this->LimLi(s, d);
        //std::cout << "n0 s = "<<s<<" d = "<< d <<std::endl;
      }
      return _LimLimap[s + 1][d];
    }

    Type operator()(const int &s, const Type &d) {
      return this->limLi(s, d);
    }

  private:
    std::vector< std::map < Type, Type > > _LimLimap;
    typename std::map< Type, Type >::iterator _it;

    unsigned _cnt;
};


template <class Type>
Type LimLimap<Type>::LimLi(const int &s, const Type &x) {
  if(x < 0) return Type(0);
  else if(s != 0) return -pow(x, s) / factorial<Type>(s);
  else if(x > 0) return Type(-1);
  else return Type(-0.5);
}

/////////////////////////////////////////////////////////////////////////

template <class Type>
class LSImap : public LimLimap <Type> {
  public:

    LSImap(const unsigned &mMax, const unsigned &sMax = 0, const unsigned &ds = 0) : LimLimap <Type> (sMax + mMax + 1u) {
      _LSImap.resize(2u + sMax + ds);
      unsigned max = 2u + mMax + sMax;

      for(unsigned s = 0; s < _LSImap.size(); s++) _LSImap[s].resize(max);
      _cnt = 0;
    };
    ~LSImap() {};

    Type LSIm1(const int &m, const Type &a, Type d);
    Type LSI(const int &s, const unsigned &m, const Type &a, Type d);

    void clear() {
      LimLimap<Type>::clear();
      for(unsigned s = 0; s < _LSImap.size(); s++) {
        for(unsigned i = 0; i < _LSImap[s].size(); i++) {
          _LSImap[s][i].clear();
        }
      }
      _cnt = 0;
    };

    void printCounter() {
      LimLimap<Type>::printCounter();
      std::cout << "LSI counter = " << _cnt << std::endl;
    }

    Type lsi(const int &s, const unsigned &m, const Type &a, const Type &d) {
      _index = std::make_pair(a, d);
      _it = _LSImap[s + 1][m].find(_index);
      if(_it == _LSImap[s + 1][m].end()) {
        _cnt++;
        _LSImap[s + 1][m][_index] = LSI(s, m, a, d);
        //std::cout << "n1 s = " << s << " m = " << m << " a = " << a << " d = " << d << std::endl;
      }
      return _LSImap[s + 1][m][_index];
    }

    Type operator()(const int &s, const unsigned &m, const Type &a, const Type &d) {
      return this->lsi(s, m, a, d);
    }

  private:
    std::vector<std::vector<std::map < std::pair<Type, Type>, Type > > > _LSImap;
    typename std::map< std::pair<Type, Type>, Type >::iterator _it;
    std::pair<Type, Type> _index;

    unsigned _cnt;
};


template <class Type>
Type LSImap<Type>::LSIm1(const int &m, const Type &a, Type d) {

  if(a == 0) {
    std::cout << "Something is wrong! The function LSIm1 can not be called with a = 0" << std::endl;
    abort();
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
Type LSImap<Type>::LSI(const int &s, const unsigned &m, const Type &a, Type d) {

  switch(s) {
    case -1:
      return LSIm1(m, a, d);
      break;

    default:

      if(a == 0) {
        return -this->limLi(s, d) / Type(m + 1);
      }

      Type INT(0);
      Type x(a + d);

      if(x < 0 || d < 0) { // in all these cases no-significant digit cancellation occurs
        if(x >= 0) {
          Type g =  1 / (-a);
          for(unsigned i = 1; i <= m + 1; i++) {
            INT += g * this->limLi(s + i, x);
            g *= (m + 1 - i) / (-a);
          }
        }
        else if(d >= 0.) {
          INT += factorial<Type>(m) * this->limLi(s + m + 1, d) / (a * pow(-a, m));
        }
      }
      else { //alternative formula to avoid significant digit cancellations when s>1, and (a+d) and d are non-negative and d >> a
        for(int i = 1; i <= s; i++) {
          INT -= pow(-a, s - i) / factorial<Type>(m + 1 + s - i) * this->limLi(i, x) ;
        }
        INT += pow(-a, s) / factorial<Type>(m + 1 + s);
        INT *= factorial<Type>(m);
      }
      return INT;
  }
}




#endif
