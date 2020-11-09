#ifndef __femus_RefineElement_hpp__
#define __femus_RefineElement_hpp__

#include "OctTreeElement.hpp"

class RefineElement {
  public:
    RefineElement(unsigned const &lmax, const char* geom_elem, const char* fe_order, const char* order_gauss_fem1,
                  const char* order_gauss_fem2, const char* gauss_type = "legendre");
    ~RefineElement();
    const std::vector<std::vector < std::vector < std::pair < unsigned, double> > > > & GetProlongationMatrix();

    void BuildElement1Prolongation(const unsigned &level, const unsigned &i);
    void SetConstants(const double &eps);

    double GetSmoothStepFunction(const double &dg1) const {
      if(dg1 < - _eps)
        return 0.;
      else if(dg1 < _eps) {
        double dg2 = dg1 * dg1;
        return (_a0 + dg1 * (_a1 + dg2 * (_a3 + dg2 * (_a5 + dg2 * (_a7 + dg2 * _a9)))));
      }
      else
        return 1.;
    };

    const unsigned &GetMaxLevel() const {
      return _lmax;
    }
        
    const elem_type &GetFem1() const {
      return *_finiteElement1;
    }
    const elem_type &GetFem2() const {
      return *_finiteElement2;
    }

    const elem_type *GetFem2Pointer() const {
      return _finiteElement2;
    }

    const unsigned &GetNumberOfNodes() const {
      return _numberOfNodes;
    }
    const unsigned &GetNumberOfLinearNodes() const {
      return _numberOfLinearNodes;
    }
    const unsigned &GetNumberOfChildren() const {
      return _numberOfChildren;
    }
    const unsigned &GetDimension() const {
      return _dim;
    }

    const OctTreeElement& GetOctTreeElement1() const {
      return _octTreeElement1;
    }

    void InitElement1(std::vector<std::vector<double>> &xv) {
      _xv1l.resize(_lmax);
      _xv1l[0].resize(1);
      _xv1l[0][0] = xv;
      for(unsigned l = 1; l < _lmax; l++) {
        _xv1l[l].resize(_numberOfChildren);
        for(unsigned i = 0; i < _numberOfChildren; i++) {
          _xv1l[l][i].resize(_dim);
          for(unsigned k = 0; k < _dim; k++) {
            _xv1l[l][i][k].resize(_numberOfNodes);
          }
        }
      }
    };

    const std::vector<std::vector<double>> & GetElement1NodeCoordinates(const unsigned &level, const unsigned &i)const {
      return _xv1l[level][i];
    }
    
    const unsigned GetLevelSize(const unsigned &level)const {
      return _xv1l[level].size();
    }
    

    const double &GetEps() const {
      return _eps;
    }

  private:
    unsigned _lmax;  
    unsigned _dim;
    unsigned _numberOfChildren;
    unsigned _numberOfNodes;
    unsigned _numberOfLinearNodes;
    const elem_type *_finiteElement1;
    const elem_type *_finiteElement2;
    const elem_type *_finiteElementLinear;

    OctTreeElement _octTreeElement1;
    std::vector<std::vector < std::vector < std::pair < unsigned, double> > > > _PMatrix;
    void BuildPMat();
    basis* _basis;
    std::vector< std::vector < std::vector < std::vector <double> > > > _xv1l;/*
    std::vector< std::vector < std::vector < std::vector <double> > > > _xi1l;*/

    double _a0, _a1, _a3, _a5, _a7, _a9, _eps;

};


RefineElement::RefineElement(unsigned const &lmax, const char* geom_elem, const char* fe_order,
                             const char* order_gauss_fem1, const char* order_gauss_fem2, const char* gauss_type) {
  _lmax = lmax;  
    
  if(!strcmp(geom_elem, "line")) {
    _numberOfChildren = 2;

    _finiteElement1 = new const elem_type_1D(geom_elem, fe_order, order_gauss_fem1, gauss_type);
    _finiteElement2 = new const elem_type_1D(geom_elem, fe_order, order_gauss_fem2, gauss_type);;
    _finiteElementLinear = new const elem_type_1D(geom_elem, "linear", "zero", gauss_type);
  }
  else if(!strcmp(geom_elem, "quad") || !strcmp(geom_elem, "tri")) {
    _numberOfChildren = 4;
    _finiteElement1 = new const elem_type_2D(geom_elem, fe_order, order_gauss_fem1, gauss_type);
    _finiteElement2 = new const elem_type_2D(geom_elem, fe_order, order_gauss_fem2, gauss_type);
    _finiteElementLinear = new const elem_type_2D(geom_elem, "linear", "zero", gauss_type);
  }
  else if(!strcmp(geom_elem, "hex") || !strcmp(geom_elem, "wedge") || !strcmp(geom_elem, "tet")) {
    _numberOfChildren = 8;
    _finiteElement1 = new const elem_type_3D(geom_elem, fe_order, order_gauss_fem1, gauss_type);
    _finiteElement2 = new const elem_type_3D(geom_elem, fe_order, order_gauss_fem2, gauss_type);
    _finiteElementLinear = new const elem_type_3D(geom_elem, "linear", "zero", gauss_type);
  }

  _dim = _finiteElement1->GetDim();
  _numberOfNodes = _finiteElement1->GetNDofs();

  _basis = _finiteElement1->GetBasis();

  _numberOfLinearNodes = _finiteElementLinear->GetNDofs();

  BuildPMat();

  _xv1l.resize(1);
  _xv1l[0].resize(1);
  _xv1l[0][0].resize(_dim);
  for(unsigned k = 0; k < _dim; k++) {
    _xv1l[0][0][k].resize(_numberOfNodes);
    for(unsigned i = 0; i < _numberOfNodes; i++) {
      _xv1l[0][0][k][i] =  *(_basis->GetXcoarse(i) + k);
    }
  }

  _octTreeElement1.Init(_xv1l[0][0], _PMatrix, _finiteElement1, _lmax);

  delete _finiteElementLinear;

}

RefineElement::~RefineElement() {
  delete _finiteElement1;
  delete _finiteElement2;
}

const std::vector<std::vector < std::vector < std::pair < unsigned, double> > > > & RefineElement::GetProlongationMatrix() {
  return _PMatrix;
}

void RefineElement::BuildPMat() {

  std::vector<std::vector < std::vector < std::pair < unsigned, double> > > > PMatrix;

  PMatrix.resize(_numberOfChildren);
  for(unsigned i = 0; i < _numberOfChildren; i++) {
    PMatrix[i].resize(_numberOfNodes);
  }

  std::vector< std::vector< double > > xvLinearChild(_dim);
  for(unsigned k = 0; k < _dim; k++) {
    xvLinearChild[k].resize(_numberOfLinearNodes);
  }

  std::vector <double> phiLinear(_numberOfLinearNodes);
  std::vector <double> phi(_numberOfNodes);
  std::vector< double > xiChild(_dim);

  for(unsigned i = 0; i < _numberOfChildren; i++) {
    for(unsigned j = 0; j < _numberOfLinearNodes; j++) {
      for(int k = 0; k < _dim; k++)  xvLinearChild[k][j] = *(_basis->GetXcoarse(_basis->GetFine2CoarseVertexMapping(i, j)) + k);
    }


    for(unsigned j = 0; j < _numberOfNodes; j++) {

      std::vector<double> xij(_dim);
      for(int k = 0; k < _dim; k++)  xij[k] = * (_basis->GetXcoarse(j) + k);

      _finiteElementLinear->GetPhi(phiLinear, xij);
      xiChild.assign(_dim, 0.);
      for(unsigned jj = 0; jj < _numberOfLinearNodes; jj++) {
        for(unsigned k = 0; k < _dim; k++) {
          xiChild[k] += phiLinear[jj] * xvLinearChild[k][jj];
        }
      }

      PMatrix[i][j].resize(_numberOfNodes);
      unsigned cnt = 0;
      _finiteElement1->GetPhi(phi, xiChild);
      for(unsigned jj = 0; jj < _numberOfNodes; jj++) {
        if(fabs(phi[jj]) > 1.0e-10) {
          PMatrix[i][j][cnt].first = jj;
          PMatrix[i][j][cnt].second = phi[jj];
          cnt++;
        }
      }
      PMatrix[i][j].resize(cnt);
    }
  }
  _PMatrix = PMatrix;
}


void RefineElement::BuildElement1Prolongation(const unsigned &level, const unsigned &i) {

  std::vector< std::vector < std::vector <double> > >::iterator xCi;
  std::vector < std::vector<double>>::iterator xCik;
  std::vector<double>::iterator xCikj;

  std::vector < std::vector<double>>::const_iterator xFk;

  std::vector <std::vector <std::vector < std::pair<unsigned, double>>>>::const_iterator Pi;
  std::vector <std::vector < std::pair<unsigned, double>>>::const_iterator Pij;
  std::vector < std::pair<unsigned, double>>::const_iterator Pijl;

  for(Pi = _PMatrix.begin(), xCi = _xv1l[level + 1].begin(); Pi != _PMatrix.end(); xCi++, Pi++) {
    for(xCik = (*xCi).begin(), xFk = _xv1l[level][i].begin(); xCik != (*xCi).end(); xCik++, xFk++) {
      for(Pij = (*Pi).begin(), xCikj = (*xCik).begin(); Pij != (*Pi).end(); Pij++, xCikj++) {
        *xCikj = 0.;
        for(Pijl = (*Pij).begin(); Pijl != (*Pij).end(); Pijl++) {
          *xCikj += Pijl->second * (*xFk)[Pijl->first];
        }
      }
    }
  }
}


void RefineElement::SetConstants(const double &eps) {
  _eps = eps;
  _a0 = 0.5; // 128./256.;
  _a1 = pow(eps, -1.) * 1.23046875; // 315/256.;
  _a3 = -pow(eps, -3.) * 1.640625; //420./256.;
  _a5 = pow(eps, -5.) * 1.4765625; // 378./256.;
  _a7 = -pow(eps, -7.) * 0.703125; // 180./256.;
  _a9 = pow(eps, -9.) * 0.13671875; // 35./256.;
}


#endif