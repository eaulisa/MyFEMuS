#ifndef __femus_RefineElement_hpp__
#define __femus_RefineElement_hpp__

class RefineElement {
  public:
    RefineElement(const char* geom_elem, const char* fe_order, const char* order_gauss_coarse,
                  const char* order_gauss_medium, const char* order_gauss_fine, const char* gauss_type = "legendre");
    ~RefineElement();
    const std::vector<std::vector < std::vector < std::pair < unsigned, double> > > > & GetProlongationMatrix();

    void BuildElement1Prolongation(const unsigned &level, const unsigned &i);
    void BuildElement2Prolongation(const unsigned &level, const unsigned &i);

    void SetConstants(const double &eps, const double &eps0);

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


    void InitElement1(std::vector<std::vector<double>> &xv, const unsigned &lMax) {

      _xv1l.resize(lMax);
      _xi1l.resize(lMax);
      for(unsigned l = 0; l < lMax; l++) {
        _xv1l[l].resize(_numberOfChildren);
        _xi1l[l].resize(_numberOfChildren);
        for(unsigned i = 0; i < _numberOfChildren; i++) {
          _xv1l[l][i].resize(_dim);
          _xi1l[l][i].resize(_dim);
          for(unsigned k = 0; k < _dim; k++) {
            _xv1l[l][i][k].resize(_numberOfNodes);
            _xi1l[l][i][k].resize(_numberOfNodes);
          }
        }
      }

      _xv1l[0][0] = xv;
      for(unsigned k = 0; k < _dim; k++) {
        for(unsigned i = 0; i < _numberOfNodes; i++) {
          _xi1l[0][0][k][i] = *(_basis->GetXcoarse(i) + k);
        }
      }
    };



    void InitElement2(std::vector<std::vector<double>> &xv, const unsigned &lMax) {

      _xv2l.resize(lMax);
      _xi2l.resize(lMax);
      for(unsigned l = 0; l < lMax; l++) {
        _xv2l[l].resize(_numberOfChildren);
        _xi2l[l].resize(_numberOfChildren);
        for(unsigned i = 0; i < _numberOfChildren; i++) {
          _xv2l[l][i].resize(_dim);
          _xi2l[l][i].resize(_dim);
          for(unsigned k = 0; k < _dim; k++) {
            _xv2l[l][i][k].resize(_numberOfNodes);
            _xi2l[l][i][k].resize(_numberOfNodes);
          }
        }
      }

      _xv2l[0][0] = xv;
      for(unsigned k = 0; k < _dim; k++) {
        for(unsigned i = 0; i < _numberOfNodes; i++) {
          _xi2l[0][0][k][i] = *(_basis->GetXcoarse(i) + k);
        }
      }
    };


    const std::vector<std::vector<double>> & GetElement1NodeCoordinates(const unsigned &level, const unsigned &i)const {
      return _xv1l[level][i];
    }

    const std::vector<std::vector<double>> & GetElement1LocalCoordinates(const unsigned &level, const unsigned &i)const {
      return _xi1l[level][i];
    }



    const std::vector<std::vector<double>> & GetElement2NodeCoordinates(const unsigned &level, const unsigned &i)const {
      return _xv2l[level][i];
    }

    const std::vector<std::vector<double>> & GetElement2LocalCoordinates(const unsigned &level, const unsigned &i)const {
      return _xi2l[level][i];
    }








    const double &GetEps() {
      return _eps;
    }
    const double &GetEps0() {
      return _eps0;
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
    std::vector< std::vector < std::vector < std::vector <double> > > > _xv2l;
    std::vector< std::vector < std::vector < std::vector <double> > > > _xi2l;

    std::vector< std::vector < std::vector < std::vector <double> > > > _xv1l;
    std::vector< std::vector < std::vector < std::vector <double> > > > _xi1l;

    double _a0, _a1, _a3, _a5, _a7, _a9, _eps, _eps0;

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

  for(Pi = _PMatrix.begin(), xCi = _xi1l[level + 1].begin(); Pi != _PMatrix.end(); xCi++, Pi++) {
    for(xCik = (*xCi).begin(), xFk = _xi1l[level][i].begin(); xCik != (*xCi).end(); xCik++, xFk++) {
      for(Pij = (*Pi).begin(), xCikj = (*xCik).begin(); Pij != (*Pi).end(); Pij++, xCikj++) {
        *xCikj = 0.;
        for(Pijl = (*Pij).begin(); Pijl != (*Pij).end(); Pijl++) {
          *xCikj += Pijl->second * (*xFk)[Pijl->first];
        }
      }
    }
  }
}







void RefineElement::BuildElement2Prolongation(const unsigned &level, const unsigned &i) {

  std::vector< std::vector < std::vector <double> > >::iterator xCi;
  std::vector < std::vector<double>>::iterator xCik;
  std::vector<double>::iterator xCikj;

  std::vector < std::vector<double>>::const_iterator xFk;

  std::vector <std::vector <std::vector < std::pair<unsigned, double>>>>::const_iterator Pi;
  std::vector <std::vector < std::pair<unsigned, double>>>::const_iterator Pij;
  std::vector < std::pair<unsigned, double>>::const_iterator Pijl;

  for(Pi = _PMatrix.begin(), xCi = _xv2l[level + 1].begin(); Pi != _PMatrix.end(); xCi++, Pi++) {
    for(xCik = (*xCi).begin(), xFk = _xv2l[level][i].begin(); xCik != (*xCi).end(); xCik++, xFk++) {
      for(Pij = (*Pi).begin(), xCikj = (*xCik).begin(); Pij != (*Pi).end(); Pij++, xCikj++) {
        *xCikj = 0.;
        for(Pijl = (*Pij).begin(); Pijl != (*Pij).end(); Pijl++) {
          *xCikj += Pijl->second * (*xFk)[Pijl->first];
        }
      }
    }
  }

  for(Pi = _PMatrix.begin(), xCi = _xi2l[level + 1].begin(); Pi != _PMatrix.end(); xCi++, Pi++) {
    for(xCik = (*xCi).begin(), xFk = _xi2l[level][i].begin(); xCik != (*xCi).end(); xCik++, xFk++) {
      for(Pij = (*Pi).begin(), xCikj = (*xCik).begin(); Pij != (*Pi).end(); Pij++, xCikj++) {
        *xCikj = 0.;
        for(Pijl = (*Pij).begin(); Pijl != (*Pij).end(); Pijl++) {
          *xCikj += Pijl->second * (*xFk)[Pijl->first];
        }
      }
    }
  }
}

void RefineElement::SetConstants(const double &eps, const double &eps0) {
  _eps = eps;
  _eps0 = eps0;
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
