#ifndef __femus_cut_fem_hpp__
#define __femus_cut_fem_hpp__

template <class TypeIO, class TypeA>
class cutFEMmap {
  public:
    virtual void clear() = 0;
    virtual TypeIO operator()(const int &s, const std::vector<unsigned> &m, const std::vector<TypeIO> &a, const TypeIO &d) = 0;

};



#endif
