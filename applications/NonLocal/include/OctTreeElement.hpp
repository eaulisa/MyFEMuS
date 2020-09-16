#ifndef __femus_OctTreeElement_hpp__
#define __femus_OctTreeElement_hpp__

class OctTreeElement {
    //specific data members
    //std::vector<std::vector < double > >  _xv;
    std::vector<std::vector < double > >  _xi;
    std::vector< OctTreeElement >  _child;
    bool _haveChilderen = false;
    //std::vector < std::vector < double> > _xg[2];
    std::vector < std::vector < double > > _phi;
    //std::vector < double> _weight[2];
    bool _femQuantitiesAreBuilt;

    std::vector<std::vector < std::vector < std::pair < unsigned, double> > > > _fakePMatrix;

    //share data members
    static unsigned _counter;
    static unsigned _counterFEM;
    std::vector<std::vector < std::vector < std::pair < unsigned, double> > > > &_PMatrix  = _fakePMatrix;
    const elem_type *_fem;

  public:

    const unsigned GetNumberOfChildren() {
      return _PMatrix.size();
    }

    const std::vector < std::vector < double> > & GetGaussShapeFunctions() {
      if(!_femQuantitiesAreBuilt) {
        BuildFemQuantities();
      }
      return _phi;
    }

    void Init(const std::vector<std::vector < double > >  &xi,
              const std::vector<std::vector < std::vector < std::pair < unsigned, double> > > > &PMatrix,
              const elem_type *fem, const unsigned &prelocateMemoryLevelMax) {

      this->Clear();
      _xi = xi;
      _PMatrix = PMatrix;

      _counter = 1;
      _counterFEM = 0;

      _haveChilderen = false;
      _fem = fem;
      _femQuantitiesAreBuilt = false;

      if(prelocateMemoryLevelMax > 0) {
        unsigned ng = _fem->GetGaussPointNumber();
        unsigned numberOfChildren = _PMatrix.size();
        this->PrelocateMemory(0, prelocateMemoryLevelMax, ng, numberOfChildren);
      }
    }

  private:

    void PrelocateMemory(const unsigned &level, const unsigned &prelocateMemoryLevelMax,
                         const unsigned &ng, const unsigned &numberOfChildren) {

      _phi.resize(ng);

      for(unsigned ig = 0; ig < ng; ig++) {

        _phi[ig].resize(_xi[0].size());
      }

      if(prelocateMemoryLevelMax > level + 1) {
        _child.resize(numberOfChildren);
        for(unsigned i = 0; i < numberOfChildren; i++) {

          _child[i]._xi.resize(_xi.size());
          for(unsigned j = 0; j < _xi.size(); j++) {

            _child[i]._xi[j].resize(_xi[j].size());
          }
          _child[i]._haveChilderen = false;
          _child[i]._femQuantitiesAreBuilt = false;

          _child[i].PrelocateMemory(level + 1, prelocateMemoryLevelMax, ng, numberOfChildren);
        }
      }
    }

    void BuildFemQuantities() {
      _counterFEM++;
      unsigned dim = _xi.size();
      unsigned numberOfNodes = _xi[0].size();

      std::vector < double>  xiFg(dim); //local coordinates of the Father element evaluated in Child element
      const double *phiC; //Shape functions of the Child element

      unsigned ng = _fem->GetGaussPointNumber();

      _phi.resize(ng);


      for(unsigned ig = 0; ig < ng; ig++) {
        phiC = _fem->GetPhi(ig);
        xiFg.assign(dim, 0.);
        for(unsigned k = 0; k < dim; k++) {
          for(unsigned j = 0; j < numberOfNodes; j++) {
            xiFg[k] += _xi[k][j] * phiC[j];
          }
        }
        _phi[ig].resize(numberOfNodes);
        _fem->GetPhi(_phi[ig], xiFg);
      }
      _femQuantitiesAreBuilt  = true;
    }

  public:
    void Clear() {
      if(_haveChilderen) {
        for(unsigned i = 0; i < _child.size(); i++) {
          _child[i].Clear();
        }
      }
      _haveChilderen = false;
      _femQuantitiesAreBuilt = false;
    }

    OctTreeElement *GetElement(const std::vector<unsigned> &genealogy) {
      return this->GetElement(genealogy, 0);
    }

  private:
    OctTreeElement *GetElement(const std::vector<unsigned> &genealogy, unsigned level) {
      if(level == genealogy.size()) {
        return this;
      }
      else {
        if(this->_haveChilderen == false) {
          this->GenerateChildren();
        }
        return this->_child[ genealogy[level] ].GetElement(genealogy, level + 1);
      }
    }

    void GenerateChildren() {

      _counter += _PMatrix.size();
      _child.resize(_PMatrix.size());

      for(unsigned i = 0; i < _child.size(); i++) {
        _child[i]._xi.resize(_xi.size());
        for(unsigned j = 0; j < _xi.size(); j++) {
          _child[i]._xi[j].resize(_xi[j].size());
        }
        _child[i]._PMatrix = _PMatrix;
        _child[i]._haveChilderen = false;
        _child[i]._fem = _fem;
      }

      std::vector< OctTreeElement >::iterator xCi;
      std::vector < std::vector<double>>::iterator xCik;
      std::vector<double>::iterator xCikj;

      std::vector < std::vector<double>>::const_iterator xFk;

      std::vector <std::vector <std::vector < std::pair<unsigned, double>>>>::const_iterator Pi;
      std::vector <std::vector < std::pair<unsigned, double>>>::const_iterator Pij;
      std::vector < std::pair<unsigned, double>>::const_iterator Pijl;

      for(Pi = _PMatrix.begin(), xCi = _child.begin(); Pi != _PMatrix.end(); xCi++, Pi++) {
        for(xCik = (*xCi)._xi.begin(), xFk = _xi.begin(); xCik != (*xCi)._xi.end(); xCik++, xFk++) {
          for(Pij = (*Pi).begin(), xCikj = (*xCik).begin(); Pij != (*Pi).end(); Pij++, xCikj++) {
            *xCikj = 0.;
            for(Pijl = (*Pij).begin(); Pijl != (*Pij).end(); Pijl++) {
              *xCikj += Pijl->second * (*xFk)[Pijl->first];
            }
          }
        }
      }
      _haveChilderen = true;
    }

    void PrintElement(std::ofstream & fout) const {
      for(unsigned j = 0; j < 4; j++) {
        fout << _xi[0][j] << " " << _xi[1][j] << std::endl;
      }
      fout << _xi[0][0] << " " << _xi[1][0] << std::endl << std::endl;
      if(_haveChilderen) {
        for(unsigned i = 0; i < _child.size(); i++) {
          _child[i].PrintElement(fout);
        }
      }
    }

  public:
    void PrintElement(std::string filename, bool clearFile = true) const {
      std::ofstream fout;
      if(clearFile) {
        fout.open(filename);
      }
      else {
        fout.open(filename, std::ios::app);
      }
      PrintElement(fout);
      fout.close();
    }

    void PrintCounter() const {
      std::cout << "The total number of elements formed is " << _counter << std::endl;
      std::cout << "The total number of FEM formed is " << _counterFEM << std::endl;
    }

};

unsigned OctTreeElement::_counter = 0;
unsigned OctTreeElement::_counterFEM = 0;

#endif
