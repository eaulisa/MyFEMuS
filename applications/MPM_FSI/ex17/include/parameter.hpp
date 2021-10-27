class parameter {
  public:
    typedef bool (*BoundaryFunc) (const std::vector < double >& x, const char name[], double &value, const int FaceName, const double time);  
    typedef double (*TimeFunc) (const double time);
    
    parameter(bool weakP, double theta, std::vector < double > gravity,
              double GAMMA, double gammacF, double gammacS, double gammap,
              bool NeoHookean, double rhos, double rhof, double nu, double E, double muf,
              std::string mMesh, double mScale, unsigned mUniform, unsigned mAdaptive,
              std::string bMesh, double bScale, unsigned deltaUniform, 
              BoundaryFunc bdcFunction, TimeFunc timeFunction) {
      _weakP = weakP;
      _theta = theta;
      _af = 1 - theta;
      _am = _af - 0.1;
      _beta = 0.25 + 0.5 * (_af - _am);
      _gamma = 0.5 + (_af - _am);
      _gravity = gravity;

      _GAMMA = GAMMA;
      _gammacF = gammacF;
      _gammacS = gammacS;
      _gammap = gammap;
      _gammau = 0.05 * _gammacF;

      _NeoHookean = NeoHookean;
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
    bool _weakP;
    double _theta, _af, _am, _beta, _gamma;
    std::vector < double > _gravity;

    double _GAMMA;
    double _gammacF;
    double _gammacS;
    double _gammap;
    double _gammau;

    bool _NeoHookean;
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
