#ifndef __femus_cut_fem_hpp__
#define __femus_cut_fem_hpp__

template <class TypeA>
class CutFEMmap {
  public:
    virtual void clear() = 0;
    virtual TypeA operator()(const int &s, const std::vector<unsigned> &m, const std::vector<TypeA> &a, const TypeA &d) = 0;
    virtual ~CutFEMmap(){}; 
};



#endif
