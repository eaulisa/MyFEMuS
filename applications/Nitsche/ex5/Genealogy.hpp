#ifndef __femus_GenealogyTree_hpp__
#define __femus_GenealogyTree_hpp__

class GenealogyTree {
    std::vector<std::vector < double > >  _xv;
    std::vector<std::vector < double > >  _xi;
    std::vector< GenealogyTree >  _child;
    bool _haveChilderen;

    std::vector<std::vector < std::vector < std::pair < unsigned, double> > > > _fakePMatrix;
    std::vector<std::vector < std::vector < std::pair < unsigned, double> > > > &_PMatrix = _fakePMatrix;
public:
    void Init(const std::vector<std::vector < double > >  &xv,
              const std::vector<std::vector < double > >  &xi,
              const std::vector<std::vector < std::vector < std::pair < unsigned, double> > > > &PMatrix) {

      this->Clear();
      _xv = xv;
      _xi = xi;
      _PMatrix = PMatrix;
      _child.resize(0);

    }

    void Clear() {
      if(_haveChilderen) {
        for(unsigned k = 0; k < _child.size(); k++) {
          _child[k].Clear();
        }
      }
      _haveChilderen = false;
    }

    GenealogyTree &GetElementQuantities(const std::vector<unsigned> &childPosition, unsigned level = 0) {
      if(level == childPosition.size()) {
        return *this;
      }
      else {
        if(_haveChilderen == false) {
          GenerateChildren(_PMatrix);
        }
        return _child[ childPosition[level] ].GetElementQuantities(childPosition, level + 1);
      }
    }

    void GenerateChildren(const std::vector<std::vector < std::vector < std::pair < unsigned, double> > > > &PMatrix) {

      _PMatrix = PMatrix;
      _child.resize(PMatrix.size());

      for(unsigned i = 0; i < _child.size(); i++) {
        _child[i]._xv.resize(_xv.size());
        _child[i]._xi.resize(_xi.size());
        for(unsigned j = 0; j < _xv.size(); j++) {
          _child[i]._xv[j].resize(_xv[j].size());
          _child[i]._xi[j].resize(_xi[j].size());
        }
        _child[i]._haveChilderen = false;
      }

      std::vector< GenealogyTree >::iterator xCi;
      std::vector < std::vector<double>>::iterator xCik;
      std::vector<double>::iterator xCikj;

      std::vector < std::vector<double>>::const_iterator xFk;

      std::vector <std::vector <std::vector < std::pair<unsigned, double>>>>::const_iterator Pi;
      std::vector <std::vector < std::pair<unsigned, double>>>::const_iterator Pij;
      std::vector < std::pair<unsigned, double>>::const_iterator Pijl;

      for(Pi = PMatrix.begin(), xCi = _child.begin(); Pi != PMatrix.end(); xCi++, Pi++) {
        for(xCik = (*xCi)._xv.begin(), xFk = _xv.begin(); xCik != (*xCi)._xv.end(); xCik++, xFk++) {
          for(Pij = (*Pi).begin(), xCikj = (*xCik).begin(); Pij != (*Pi).end(); Pij++, xCikj++) {
            *xCikj = 0.;
            for(Pijl = (*Pij).begin(); Pijl != (*Pij).end(); Pijl++) {
              *xCikj += Pijl->second * (*xFk)[Pijl->first];
            }
          }
        }
      }

      for(Pi = PMatrix.begin(), xCi = _child.begin(); Pi != PMatrix.end(); xCi++, Pi++) {
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

    void Print(std::string filename, bool clearFile = true) {
      std::ofstream fout;
      if(clearFile) {
        fout.open(filename);
      }
      else {
        fout.open(filename, std::ios::app);
      }
      for(unsigned j = 0; j < 4; j++) {
        fout << _xv[0][j] << " " << _xv[1][j] << std::endl;
      }
      fout << _xv[0][0] << " " << _xv[1][0] << std::endl << std::endl;
      fout.close();
      if(_haveChilderen) {
        this->Print(filename, false);
      }
    }




};
#endif
