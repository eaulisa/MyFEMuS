#ifndef __femus_cut_fem_LSI_hpp__
#define __femus_cut_fem_LSI_hpp__

#include <iostream>
#include <boost/math/special_functions/factorials.hpp>

#include <map>

using boost::math::factorial;

template <class Type>
class LSImap {
  public:

    LSImap() {};
    ~LSImap() {};

    Type LimLi(const int &s, const Type &x);
    Type LSIm1(const int &m, const Type &a, Type d);
    Type LSI0(const int &m, const Type &a, Type d);
    Type LSI(const int &s, const unsigned &m, const Type &a, Type d);
    
    void clear() {
      _LSImap.clear();
    };

    Type operator()(const int &s, const unsigned &m, const Type &a, const Type &d) {
      _index = std::make_pair(std::make_pair(s, m), std::make_pair(a, d));
      _it = _LSImap.find(_index);
      if(_it == _LSImap.end()) {
        std::cout << "new ";
        _LSImap[_index] = LSI(s, m, a, d);
      }
      return _LSImap[_index];
    }

  private:
    std::map< std::pair< std::pair<int, unsigned>, std::pair<Type, Type> >, Type > _LSImap;
    typename std::map< std::pair< std::pair<int, unsigned>, std::pair<Type, Type> >, Type >::iterator _it;
    std::pair< std::pair<int, unsigned>, std::pair<Type, Type> > _index;
};


template <class Type>
Type LSImap<Type>::LimLi(const int &s, const Type &x) {

  if(x < 0) return Type(0);
  else if(s != 0) return -pow(x, s) / factorial<Type>(s);
  else if(x > 0) return Type(-1);
  else return Type(-0.5);
}


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
Type LSImap<Type>::LSI0(const int &m, const Type &a, Type d) {

  if(a == 0) {
    return -LimLi(0, d) / Type(m + 1);
  }

  d = -d / a;


  if(a > 0) {
    if(d <= 0) return 1 / Type(m + 1);
    else if(d < 1) return (Type(1) - pow(d, m + 1)) / Type(m + 1);
    else return Type(0);
  }
  else if(a < 0) {
    if(d >= 1) return 1 / Type(m + 1);
    else if(d > 0) return pow(d, m + 1) / Type(m + 1);
    else return 0;
  }
}

template <class Type>
Type LSImap<Type>::LSI(const int &s, const unsigned &m, const Type &a, Type d) {

  Type I1;

  switch(s) {
    case -1:
      return LSIm1(m, a, d);
      break;

//     case 0:
//       return LSI0(m, a, d);
//       break;

    default:

      if(a == 0) {
        return -LimLi(s, d) / Type(m + 1);
      }

      Type INT(0);
      Type x(a + d);

      if(x < 0 || d < 0) { // in all these cases no-significant digit cancellation occurs
        if(x >= 0) {
          Type g =  1 / (-a);
          for(unsigned i = 1; i <= m + 1; i++) {
            INT += g * LimLi(s + i, x);
            g *= (m + 1 - i) / (-a);
          }
        }
        else if(d >= 0.) {
          INT += factorial<Type>(m) * LimLi(s + m + 1, d) / (a * pow(-a, m));
        }
      }
      else { //alternative formula to avoid significant digit cancellations when s>1, and (a+d) and d are non-negative and d >> a
        for(int i = 1; i <= s; i++) {
          INT -= pow(-a, s - i) / factorial<Type>(m + 1 + s - i) * LimLi(i, x) ;
        }
        INT += pow(-a, s) / factorial<Type>(m + 1 + s);
        INT *= factorial<Type>(m);
      }
      return INT;
  }
}




#endif
