#ifndef __femus_RefineElement_hpp__
#define __femus_RefineElement_hpp__

class RefineElement {
  public:
    RefineElement(const char* geom_elem, const char* fe_order, const char* order_gauss_coarse,
                  const char* order_gauss_medium, const char* order_gauss_fine, const char* gauss_type = "legendre");
    ~RefineElement();
    const std::vector<std::vector < std::vector < std::pair < unsigned, double> > > > & GetProlongationMatrix();

    void BuildElementProlongation(const unsigned &level, const unsigned &i);

    void SetConstants(const double &eps);

    double GetSmoothStepFunction(const double &dg1);

    const elem_type &GetFEMCoarse() const {
      return *_finiteElementCoarse;
    }
    const elem_type &GetFEMMedium() const {
      return *_finiteElementMedium;
    }
    const elem_type &GetFEMFine() const {
      return *_finiteElementFine;
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

    void InitElement(std::vector<std::vector<double>> &xv, const unsigned &lMax) {

      _xvl.resize(lMax);
      _xil.resize(lMax);
      for(unsigned l = 0; l < lMax; l++) {
        _xvl[l].resize(_numberOfChildren);
        _xil[l].resize(_numberOfChildren);
        for(unsigned i = 0; i < _numberOfChildren; i++) {
          _xvl[l][i].resize(_dim);
          _xil[l][i].resize(_dim);
          for(unsigned k = 0; k < _dim; k++) {
            _xvl[l][i][k].resize(_numberOfNodes);
            _xil[l][i][k].resize(_numberOfNodes);
          }
        }
      }

      _xvl[0][0] = xv;
      for(unsigned k = 0; k < _dim; k++) {
        for(unsigned i = 0; i < _numberOfNodes; i++) {
          _xil[0][0][k][i] = *(_basis->GetXcoarse(i) + k);
        }
      }
    };

    const std::vector<std::vector<double>> & GetNodeCoordinates(const unsigned &level, const unsigned &i)const {
      return _xvl[level][i];
    }

    const std::vector<std::vector<double>> & GetNodeLocalCoordinates(const unsigned &level, const unsigned &i)const {
      return _xil[level][i];
    }

    const double &GetEps() {
      return _eps;
    }

  private:
    unsigned _dim;
    unsigned _numberOfChildren;
    unsigned _numberOfNodes;
    unsigned _numberOfLinearNodes;
    const elem_type *_finiteElementCoarse;
    const elem_type *_finiteElementMedium;
    const elem_type *_finiteElementFine;
    const elem_type *_finiteElementLinear;
    std::vector<std::vector < std::vector < std::pair < unsigned, double> > > > _PMatrix;
    void BuildPMat();
    basis* _basis;
    std::vector< std::vector < std::vector < std::vector <double> > > > _xvl;
    std::vector< std::vector < std::vector < std::vector <double> > > > _xil;
    double _a0, _a1, _a3, _a5, _a7, _a9, _eps;

};


RefineElement::RefineElement(const char* geom_elem, const char* fe_order, const char* order_gauss_coarse,
                             const char* order_gauss_medium, const char* order_gauss_fine, const char* gauss_type) {
  if(!strcmp(geom_elem, "line")) {
    _numberOfChildren = 2;

    _finiteElementCoarse = new const elem_type_1D(geom_elem, fe_order, order_gauss_coarse, gauss_type);
    _finiteElementMedium = new const elem_type_1D(geom_elem, fe_order, order_gauss_medium, gauss_type);
    _finiteElementFine = new const elem_type_1D(geom_elem, fe_order, order_gauss_fine, gauss_type);
    _finiteElementLinear = new const elem_type_1D(geom_elem, "linear", "zero", gauss_type);
  }
  else if(!strcmp(geom_elem, "quad") || !strcmp(geom_elem, "tri")) {
    _numberOfChildren = 4;
    _finiteElementCoarse = new const elem_type_2D(geom_elem, fe_order, order_gauss_coarse, gauss_type);
    _finiteElementMedium = new const elem_type_2D(geom_elem, fe_order, order_gauss_medium, gauss_type);
    _finiteElementFine = new const elem_type_2D(geom_elem, fe_order, order_gauss_fine, gauss_type);
    _finiteElementLinear = new const elem_type_2D(geom_elem, "linear", "zero", gauss_type);
  }
  else if(!strcmp(geom_elem, "hex") || !strcmp(geom_elem, "wedge") || !strcmp(geom_elem, "tet")) {
    _numberOfChildren = 8;
    _finiteElementCoarse = new const elem_type_3D(geom_elem, fe_order, order_gauss_coarse, gauss_type);
    _finiteElementMedium = new const elem_type_3D(geom_elem, fe_order, order_gauss_medium, gauss_type);
    _finiteElementFine = new const elem_type_3D(geom_elem, fe_order, order_gauss_fine, gauss_type);
    _finiteElementLinear = new const elem_type_3D(geom_elem, "linear", "zero", gauss_type);
  }

  _dim = _finiteElementFine->GetDim();
  _numberOfNodes = _finiteElementFine->GetNDofs();

  _basis = _finiteElementFine->GetBasis();

  _numberOfLinearNodes = _finiteElementLinear->GetNDofs();

  BuildPMat();
  delete _finiteElementLinear;

}

RefineElement::~RefineElement() {
  delete _finiteElementCoarse;  
  delete _finiteElementMedium;
  delete _finiteElementFine;
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
        xiChild[0] += phiLinear[jj] * xvLinearChild[0][jj];
        xiChild[1] += phiLinear[jj] * xvLinearChild[1][jj];
      }

      PMatrix[i][j].resize(_numberOfNodes);
      unsigned cnt = 0;
      _finiteElementFine->GetPhi(phi, xiChild);
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
//   std::cout.precision(16);
//   for(unsigned i = 0; i < _numberOfChildren; i++) {
//     for(unsigned j = 0; j < _numberOfNodes; j++) {
//       double sum = 0.;
//       for(unsigned k = 0; k < PMatrix[i][j].size(); k++) {
//         sum += PMatrix[i][j][k].second;
//         std::cout <<  PMatrix[i][j][k].first << " " << PMatrix[i][j][k].second << "\t";
//       }
//       std::cout << sum << " "<< std::endl;
//     }
//     std::cout << std::endl;
//   }
//   std::cout << std::endl;
  _PMatrix = PMatrix;
}

void RefineElement::BuildElementProlongation(const unsigned &level, const unsigned &i) {

  std::vector< std::vector < std::vector <double> > >::iterator xCi;
  std::vector < std::vector<double>>::iterator xCik;
  std::vector<double>::iterator xCikj;

  std::vector < std::vector<double>>::const_iterator xFk;

  std::vector <std::vector <std::vector < std::pair<unsigned, double>>>>::const_iterator Pi;
  std::vector <std::vector < std::pair<unsigned, double>>>::const_iterator Pij;
  std::vector < std::pair<unsigned, double>>::const_iterator Pijl;

  for(Pi = _PMatrix.begin(), xCi = _xvl[level + 1].begin(); Pi != _PMatrix.end(); xCi++, Pi++) {
    for(xCik = (*xCi).begin(), xFk = _xvl[level][i].begin(); xCik != (*xCi).end(); xCik++, xFk++) {
      for(Pij = (*Pi).begin(), xCikj = (*xCik).begin(); Pij != (*Pi).end(); Pij++, xCikj++) {
        *xCikj = 0.;
        for(Pijl = (*Pij).begin(); Pijl != (*Pij).end(); Pijl++) {
          *xCikj += Pijl->second * (*xFk)[Pijl->first];
        }
      }
    }
  }

  for(Pi = _PMatrix.begin(), xCi = _xil[level + 1].begin(); Pi != _PMatrix.end(); xCi++, Pi++) {
    for(xCik = (*xCi).begin(), xFk = _xil[level][i].begin(); xCik != (*xCi).end(); xCik++, xFk++) {
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




double RefineElement::GetSmoothStepFunction(const double &dg1) {
  double dg2 = dg1 * dg1;

  if(dg1 < - _eps)
    return 0.;
  else if(dg1 < _eps) {
    return (_a0 + dg1 * (_a1 + dg2 * (_a3 + dg2 * (_a5 + dg2 * (_a7 + dg2 * _a9)))));
  }
  else
    return 1.;
}

#endif
