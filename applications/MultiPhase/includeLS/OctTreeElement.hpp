#ifndef __femus_OctTreeElement_hpp__
#define __femus_OctTreeElement_hpp__

class OctTreeElement {
    //specific data members
    std::vector<std::vector < double > > _xi;
    std::vector< OctTreeElement > _child;
    bool _haveChilderen = false;
    std::vector < std::vector < double > > _phi;

    std::vector<std::vector < std::vector < std::pair < unsigned, double> > > > _fakePMatrix;

    //share data members
    std::vector<std::vector < std::vector < std::pair < unsigned, double> > > > &_PMatrix  = _fakePMatrix;
    const elem_type *_fem;

  public:

    const unsigned GetNumberOfChildren() const {
      return _PMatrix.size();
    }

    const std::vector < std::vector < double> > & GetGaussShapeFunctions() const {
      return _phi;
    }

    void Init(const std::vector<std::vector < double > >  &xi,
              const std::vector<std::vector < std::vector < std::pair < unsigned, double> > > > &PMatrix,
              const elem_type *fem, const unsigned &levelMax) {

      this->Clear();
      _xi = xi;
      _PMatrix = PMatrix;

      _haveChilderen = false;
      _fem = fem;

      if(levelMax > 0) {
        unsigned ng = _fem->GetGaussPointNumber();
        unsigned numberOfChildren = _PMatrix.size();
        this->RecursiveBuild(0, levelMax, ng, numberOfChildren);
      }
      else {
        std::cout << "Error! Level max should be greater than 0" << std::endl;
        abort;
      }
    }

  private:

    void RecursiveBuild(const unsigned &level, const unsigned &levelMax,
                        const unsigned &ng, const unsigned &numberOfChildren) {

      BuildFemQuantities();
      if(levelMax > level + 1) {
        GenerateChildren();
        for(unsigned i = 0; i < numberOfChildren; i++) {
          _child[i]._haveChilderen = false;
          _child[i].RecursiveBuild(level + 1, levelMax, ng, numberOfChildren);
        }
      }
    }

    void BuildFemQuantities() {
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
    }

  public:
    void Clear() {
      if(_haveChilderen) {
        for(unsigned i = 0; i < _child.size(); i++) {
          _child[i].Clear();
        }
      }
      _haveChilderen = false;
    }

    const OctTreeElement *GetElement(const std::vector<unsigned> &genealogy) const {
      return this->GetElement(genealogy, 0);
    }

  private:
    const OctTreeElement *GetElement(const std::vector<unsigned> &genealogy, unsigned level) const {
      if(level == genealogy.size()) {
        return this;
      }
      else {
        if(this->_haveChilderen == false) {
          std::cout << "Error! This level child is not available" << std::endl;
          abort;
        }
        return this->_child[ genealogy[level] ].GetElement(genealogy, level + 1);
      }
    }

    void GenerateChildren() {

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
};


#endif
