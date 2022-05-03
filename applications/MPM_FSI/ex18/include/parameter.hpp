class parameter {
  public:
    typedef bool (*BoundaryFunc)(const std::vector < double >& x, const char name[], double &value, const int FaceName, const double time);
    typedef double (*TimeFunc)(const double time);

    parameter(bool weakP, double rhoInf, std::vector < double > gravity,
              double GAMMA, double gammacF, double gammacS, double gammap,
              bool NeoHookean, bool plainStress, double rhos, double rhof, double nu, double E, double muf,
              std::string mMesh, double mScale, unsigned mUniform, unsigned mAdaptive,
              std::string bMesh, double bScale, int deltaUniform,
              BoundaryFunc bdcFunction, TimeFunc timeFunction) {
//       _weakP = weakP;

      _af = rhoInf / (rhoInf + 1.);
      _am = (2. * rhoInf - 1.) / (rhoInf + 1.);

      _beta = 0.25 * (1. + _af - _am) * (1. + _af - _am);
      _gamma = 0.5 + _af - _am;

      _theta = 1. - _af;
      _gravity = gravity;

      _GAMMA = GAMMA;
      _gammacF = gammacF;
      _gammacS = gammacS;
      _gammap = gammap;
      _gammau = 0.05 * _gammacF;

      _NeoHookean = NeoHookean;
      _plainStress = plainStress;
      _rhos = rhos;
      _rhof = rhof;
      _nu = nu;
      _E = E;
      _muf = muf;

      _mMesh = mMesh;
      _mScale = mScale;
      _mUniform = mUniform;
      _mAdaptive = mAdaptive;

      _bMesh = bMesh;
      _bScale = bScale;
      _bUniform = _mUniform + deltaUniform;

      _bdcFunction = bdcFunction;
      _timeFunction = timeFunction;

    }

  public:
//     bool _weakP;
    double _theta, _af, _am;
    double _beta, _gamma;
    std::vector < double > _gravity;

    double _GAMMA;
    double _gammacF;
    double _gammacS;
    double _gammap;
    double _gammau;

    bool _NeoHookean;
    bool _plainStress;
    double _rhos;
    double _rhof;
    double _nu;
    double _E;
    double _muf;

    std::string _mMesh;
    double _mScale;
    unsigned _mUniform;
    unsigned _mAdaptive;
    double _bScale;
    std::string _bMesh;
    unsigned _bUniform;

    BoundaryFunc _bdcFunction;
    TimeFunc _timeFunction;

};

unsigned vectorToUint(const std::vector<unsigned> &num) {
  unsigned n = 0;
  unsigned N = num.size();
  for(unsigned i = 0; i < N; i++) {
    n += num[i] * pow(10, N - i - 1);
  }
  return n;
}

unsigned mapToUint(const std::map<unsigned, bool> &num) {
  std::map<unsigned, bool>::const_iterator it;
  unsigned n = 0;
  unsigned fac = 1;
  for(it = num.begin(); it != num.end(); it++, fac *= 10) {
    n += it->first * fac;
  }
  return n;
}

std::vector<unsigned> uintToVector(unsigned n) {
  std::vector<unsigned> vec;
  while(n != 0) {
    vec.push_back(n % 10);
    n /= 10;
  }
  reverse(vec.begin(), vec.end());
  return vec;
}

std::map<unsigned, bool> uintToMap(unsigned n) {
  std::map<unsigned, bool> num;
  while(n != 0) {
    num[n % 10] = true;
    n /= 10;
  }
  return num;
}

