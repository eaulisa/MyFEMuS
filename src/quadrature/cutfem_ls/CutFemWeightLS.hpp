#ifndef FEMUS_CUT_FEM_WEIGHT_LS_HPP
#define FEMUS_CUT_FEM_WEIGHT_LS_HPP

#include "CutFemWeight.hpp"
#include "FELevelSetMoments.hpp"
#include "Hex27LevelSetMoments.hpp"
#include "Hex27DirectionSelection.hpp"
#include "AdaptiveLineSliceMoments.hpp"
#include "Quad9DirectionSelection.hpp"
#include "Tri7DirectionSelection.hpp"
#include "Tet15LevelSet.hpp"
#include "Tet15LevelSetMoments.hpp"
#include "SquareParabolaMoment.hpp"
#include <stdexcept>

template <class TypeIO, class TypeA>
class CutFemWeightLS : private CutFemWeight<TypeIO, TypeA> {
  public:
    CutFemWeightLS(const GeomElType &geomElemType, const unsigned &qM, const std::string &gaussType) : CutFemWeight<TypeIO, TypeA>(geomElemType, qM, gaussType) {
      _fPair.resize(this->_L, 2);
    }

    void build() {
      UpdateQuadratureRule(this->_qM);
    }

    void clear() {
      CutFemWeight<TypeIO, TypeA>::clear();
      _nodalMomentCacheValid = false;
      _nodalWeightCacheValid = false;
    }

    void operator()(const int &s, const std::vector<TypeA> &a, const TypeA &d, const bool &xSquared, std::vector<TypeIO> &weightCF);
    void operator()(const int &s, const std::vector<TypeA> &aX, const TypeA &dX, const std::vector<TypeA> &aY, const TypeA &dY, const bool &xSquared, std::vector<TypeIO> &weightCF);

    void GetWeightFromMoments(const std::vector<std::vector<std::vector<double>>>& moments, std::vector<TypeIO>& weightCF) {
      if(this->_geomElemType != HEX || moments.size() < this->_qM+1U) throw std::invalid_argument("Expected unit-cube HEX moments.");
      unsigned count = 0;
      for(unsigned q = 0; q <= this->_qM; ++q) {
        for(int ii = q; ii >= 0; --ii) {
          for(int jj = q-static_cast<unsigned>(ii); jj >= 0; --jj) {
            unsigned i = static_cast<unsigned>(ii), j = static_cast<unsigned>(jj), k = q-i-j;
            if(moments[i].size() <= j || moments[i][j].size() <= k || !std::isfinite(moments[i][j][k])) throw std::invalid_argument("Invalid HEX moment table.");
            this->_f[count++] = static_cast<TypeA>(moments[i][j][k]);
          }
        }
      }
      ++this->_cntCall;
      this->EvaluateMomentWeights(weightCF);
    }

    void GetWeightFromMoments(const std::vector<std::vector<double>>& moments, std::vector<TypeIO>& weightCF) {
      if (this->_dim != 2 || (this->_geomElemType != QUAD && this->_geomElemType != TRI) || moments.size() < this->_qM + 1) {
        throw std::invalid_argument("Expected canonical QUAD or TRI moments.");
      }
      unsigned count = 0;
      for (unsigned q = 0; q <= this->_qM; ++q) {
        for (unsigned j = 0; j <= q; ++j) {
          unsigned i = q - j;
          if (moments[i].size() <= j || !std::isfinite(moments[i][j])) {
            throw std::invalid_argument("Invalid canonical moment table.");
          }
          this->_f[count++] = static_cast<TypeA>(moments[i][j]);
        }
      }
      ++this->_cntCall;
      this->EvaluateMomentWeights(weightCF);
    }

    void GetWeightFromMoments(const std::vector<double>& moments, std::vector<TypeIO>& weightCF) {
      if (this->_geomElemType != LINE || moments.size() < this->_qM + 1) {
        throw std::invalid_argument("Expected canonical LINE moments.");
      }
      for (unsigned i = 0; i <= this->_qM; ++i) {
        if (!std::isfinite(moments[i])) {
          throw std::invalid_argument("Invalid canonical LINE moment.");
        }
        this->_f[i] = static_cast<TypeA>(moments[i]);
      }
      ++this->_cntCall;
      this->EvaluateMomentWeights(weightCF);
    }

    void GetWeightFromFE(const int& s, const std::vector<std::vector<double>>& coefficient, std::vector<TypeIO>& weightCF, double absoluteTolerance = 1.0e-12, double relativeTolerance = 1.0e-12) {
      if (this->_geomElemType != QUAD || (s != 0 && s != -1)) {
        throw std::invalid_argument("QUAD9 FE weights require QUAD and s=0 or s=-1.");
      }
      if (coefficient.size() != 3) {
        throw std::invalid_argument("Expected a 3 by 3 QUAD9 coefficient matrix.");
      }
      bool planar = true;
      for (unsigned i = 0; i < 3; ++i) {
        if (coefficient[i].size() != 3) {
          throw std::invalid_argument("Expected a 3 by 3 QUAD9 coefficient matrix.");
        }
        for (unsigned j = 0; j < 3; ++j) {
          if (!std::isfinite(coefficient[i][j])) {
            throw std::invalid_argument("A QUAD9 coefficient is not finite.");
          }
          if (i + j > 1 && coefficient[i][j] != 0.0) {
            planar = false;
          }
        }
      }
      if (planar && (coefficient[1][0] != 0.0 || coefficient[0][1] != 0.0)) {
        std::vector<TypeIO> normal = {static_cast<TypeIO>(coefficient[1][0]), static_cast<TypeIO>(coefficient[0][1])};
        CutFemWeight<TypeIO, TypeA>::operator()(s, normal, static_cast<TypeIO>(coefficient[0][0]), weightCF);
        return;
      }
      std::vector<std::vector<double>> positive;
      std::vector<std::vector<double>> delta;
      FELevelSetMoments::IntegrateQuad9CanonicalSelected(coefficient, this->_qM, positive, delta, absoluteTolerance, relativeTolerance, -1.0, 1.0, -1.0, 1.0, s == -1, false, _quad9DirectionValidation, _quad9IntegrationDirection, nullptr, true);
      GetWeightFromMoments((s == 0) ? positive : delta, weightCF);
    }

    void GetWeightFromFE(const int& s, const std::vector<double>& coefficient, std::vector<TypeIO>& weightCF, double absoluteTolerance = 1.0e-12, double relativeTolerance = 1.0e-12) {
      if (s != 0 && s != -1) {
        throw std::invalid_argument("FE weights require s=0 or s=-1.");
      }
      if (this->_geomElemType == TRI) {
        if (coefficient.size() != 7) {
          throw std::invalid_argument("Expected seven TRI7 coefficients.");
        }
        for (unsigned k = 0; k < 7; ++k) {
          if (!std::isfinite(coefficient[k])) {
            throw std::invalid_argument("A TRI7 coefficient is not finite.");
          }
        }
        if (coefficient[3] == 0.0 && coefficient[4] == 0.0 && coefficient[5] == 0.0 && coefficient[6] == 0.0 && (coefficient[1] != 0.0 || coefficient[2] != 0.0)) {
          std::vector<TypeIO> normal = {static_cast<TypeIO>(coefficient[1]), static_cast<TypeIO>(coefficient[2])};
          CutFemWeight<TypeIO, TypeA>::operator()(s, normal, static_cast<TypeIO>(coefficient[0]), weightCF);
          return;
        }
        std::vector<std::vector<double>> positive;
        std::vector<std::vector<double>> delta;
        FELevelSetMoments::IntegrateTri7CanonicalSelected(coefficient, this->_qM, positive, delta, absoluteTolerance, relativeTolerance, s == -1, false, _tri7DirectionValidation, _tri7IntegrationDirection, nullptr, true);
        GetWeightFromMoments((s == 0) ? positive : delta, weightCF);
      }
      else if (this->_geomElemType == LINE && coefficient.size() == 3) {
        for (unsigned k = 0; k < 3; ++k) {
          if (!std::isfinite(coefficient[k])) {
            throw std::invalid_argument("A LINE3 coefficient is not finite.");
          }
        }
        if (coefficient[0] == 0.0 && coefficient[1] != 0.0) {
          std::vector<TypeIO> normal(1, static_cast<TypeIO>(coefficient[1]));
          CutFemWeight<TypeIO, TypeA>::operator()(s, normal, static_cast<TypeIO>(coefficient[2]), weightCF);
          return;
        }
        std::vector<double> positive;
        std::vector<double> delta;
        FELevelSetMoments::Integrate1DCanonical(coefficient[0], coefficient[1], coefficient[2], this->_qM, positive, delta, -1.0, 1.0, s == -1);
        GetWeightFromMoments((s == 0) ? positive : delta, weightCF);
      }
      else {
        throw std::invalid_argument("Expected TRI7 coefficients or LINE3 coefficients {a,b,c}.");
      }
    }

    void GetWeightFromCanonicalFE(const int& s, const std::vector<std::vector<double>>& coefficient, std::vector<TypeIO>& weightCF, double absoluteTolerance = 1.0e-12, double relativeTolerance = 1.0e-12) {
      if (this->_geomElemType != QUAD || (s != 0 && s != -1)) {
        throw std::invalid_argument("Canonical QUAD9 weights require QUAD and s=0 or s=-1.");
      }
      if (coefficient.size() != 3) {
        throw std::invalid_argument("Expected a 3 by 3 unit-square QUAD9 coefficient matrix.");
      }
      bool planar = true;
      for (unsigned i = 0; i < 3; ++i) {
        if (coefficient[i].size() != 3) {
          throw std::invalid_argument("Expected a 3 by 3 unit-square QUAD9 coefficient matrix.");
        }
        for (unsigned j = 0; j < 3; ++j) {
          if (!std::isfinite(coefficient[i][j])) {
            throw std::invalid_argument("A QUAD9 coefficient is not finite.");
          }
          if (i + j > 1 && coefficient[i][j] != 0.0) {
            planar = false;
          }
        }
      }
      if (planar && (coefficient[1][0] != 0.0 || coefficient[0][1] != 0.0)) {
        std::vector<TypeIO> normal = {static_cast<TypeIO>(0.5L * coefficient[1][0]), static_cast<TypeIO>(0.5L * coefficient[0][1])};
        CutFemWeight<TypeIO, TypeA>::operator()(s, normal, static_cast<TypeIO>(static_cast<long double>(coefficient[0][0]) + 0.5L * coefficient[1][0] + 0.5L * coefficient[0][1]), weightCF);
        return;
      }
      std::vector<std::vector<double>> positive;
      std::vector<std::vector<double>> delta;
      FELevelSetMoments::IntegrateQuad9UnitSelected(coefficient, this->_qM, positive, delta, absoluteTolerance, relativeTolerance, s == -1, false, 24, _quad9IntegrationDirection, 200000, nullptr, _quad9DirectionValidation, nullptr, true);
      GetWeightFromMoments((s == 0) ? positive : delta, weightCF);
    }

    void GetWeightFromCanonicalFE(const int& s, const std::vector<double>& coefficient, std::vector<TypeIO>& weightCF, double absoluteTolerance = 1.0e-12, double relativeTolerance = 1.0e-12) {
      if (s != 0 && s != -1) {
        throw std::invalid_argument("Canonical FE weights require s=0 or s=-1.");
      }
      if (this->_geomElemType == TRI) {
        if (coefficient.size() != 7) {
          throw std::invalid_argument("Expected reflected TRI7 coefficients of 1,x,y,x*x,x*y,y*y,x*y*(x-y).");
        }
        for (unsigned k = 0; k < 7; ++k) {
          if (!std::isfinite(coefficient[k])) {
            throw std::invalid_argument("A TRI7 coefficient is not finite.");
          }
        }
        if (coefficient[3] == 0.0 && coefficient[4] == 0.0 && coefficient[5] == 0.0 && coefficient[6] == 0.0 && (coefficient[1] != 0.0 || coefficient[2] != 0.0)) {
          std::vector<TypeIO> normal = {static_cast<TypeIO>(-coefficient[1]), static_cast<TypeIO>(coefficient[2])};
          CutFemWeight<TypeIO, TypeA>::operator()(s, normal, static_cast<TypeIO>(static_cast<long double>(coefficient[0]) + coefficient[1]), weightCF);
          return;
        }
        std::vector<std::vector<double>> positive;
        std::vector<std::vector<double>> delta;
        FELevelSetMoments::IntegrateTri7ReflectedSelected(coefficient, this->_qM, positive, delta, absoluteTolerance, relativeTolerance, s == -1, false, 24, _tri7IntegrationDirection, 200000, nullptr, _tri7DirectionValidation, nullptr, true);
        GetWeightFromMoments((s == 0) ? positive : delta, weightCF);
      }
      else if (this->_geomElemType == LINE && coefficient.size() == 3) {
        for (unsigned k = 0; k < 3; ++k) {
          if (!std::isfinite(coefficient[k])) {
            throw std::invalid_argument("A LINE3 coefficient is not finite.");
          }
        }
        if (coefficient[0] == 0.0 && coefficient[1] != 0.0) {
          std::vector<TypeIO> normal(1, static_cast<TypeIO>(0.5L * coefficient[1]));
          CutFemWeight<TypeIO, TypeA>::operator()(s, normal, static_cast<TypeIO>(static_cast<long double>(coefficient[2]) + 0.5L * coefficient[1]), weightCF);
          return;
        }
        std::vector<double> positive;
        std::vector<double> delta;
        FELevelSetMoments::Integrate1D(coefficient[0], coefficient[1], coefficient[2], this->_qM, positive, delta, 0.0, 1.0, s == -1);
        GetWeightFromMoments((s == 0) ? positive : delta, weightCF);
      }
      else {
        throw std::invalid_argument("Expected reflected TRI7 coefficients or unit-interval LINE3 coefficients {a,b,c}.");
      }
    }

    void GetWeightFromCanonicalFE(const int& s, const std::vector<std::vector<std::vector<double>>>& coefficient, std::vector<TypeIO>& weightCF, double absoluteTolerance = 1.0e-12, double relativeTolerance = 1.0e-12) {
      if(this->_geomElemType != HEX || (s != 0 && s != -1)) throw std::invalid_argument("HEX27 weights require HEX and s=0 or s=-1.");
      std::vector<long double> checked;
      long double scale;
      FELevelSetMoments::Hex27Normalize(coefficient,checked,scale);
      bool planar = true;
      for(unsigned i = 0; i < 3; ++i) {
        for(unsigned j = 0; j < 3; ++j) {
          for(unsigned k = 0; k < 3; ++k) if(i+j+k > 1 && coefficient[i][j][k] != 0.) planar = false;
        }
      }
      if(planar && (coefficient[1][0][0] != 0. || coefficient[0][1][0] != 0. || coefficient[0][0][1] != 0.)) {
        std::vector<TypeIO> normal = {static_cast<TypeIO>(0.5L*coefficient[1][0][0]),static_cast<TypeIO>(0.5L*coefficient[0][1][0]),static_cast<TypeIO>(0.5L*coefficient[0][0][1])};
        CutFemWeight<TypeIO, TypeA>::operator()(s, normal,static_cast<TypeIO>(static_cast<long double>(coefficient[0][0][0])+0.5L*coefficient[1][0][0]+0.5L*coefficient[0][1][0]+0.5L*coefficient[0][0][1]),weightCF);
        return;
      }
      std::vector<std::vector<std::vector<double>>> positive, delta;
      FELevelSetMoments::IntegrateHex27UnitSelected(coefficient,this->_qM,positive,delta,absoluteTolerance,relativeTolerance,s == -1,24,_hex27IntegrationDirection,200000,nullptr,_hex27DirectionValidation);
      GetWeightFromMoments(s == 0 ? positive : delta,weightCF);
    }

    void GetWeightFromNodalValues(const int& s, const std::vector<double>& nodalValues, std::vector<TypeIO>& weightCF, double absoluteTolerance = 1.0e-12, double relativeTolerance = 1.0e-12, bool referenceSurface = false) {
      if(referenceSurface) {
        if(s != 0 && s != -1) throw std::invalid_argument("Reference surface nodal weights require s=0 or s=-1.");
        std::vector<TypeIO> positive;
        std::vector<TypeIO> surface;
        GetWeightsFromNodalValues(nodalValues,positive,surface,absoluteTolerance,relativeTolerance,true);
        weightCF = s == 0 ? positive : surface;
        return;
      }
      if(this->_geomElemType == TET) {
        if(s != 0 && s != -1) throw std::invalid_argument("TET15 nodal weights require s=0 or s=-1.");
        std::vector<TypeIO> positive;
        std::vector<TypeIO> delta;
        GetWeightsFromNodalValues(nodalValues,positive,delta,absoluteTolerance,relativeTolerance);
        weightCF = s == 0 ? positive : delta;
      }
      else if(this->_geomElemType == HEX) {
        std::vector<std::vector<std::vector<double>>> coefficient;
        FELevelSetMoments::Hex27UnitCoefficientsFromNodalValues(nodalValues,coefficient);
        GetWeightFromCanonicalFE(s,coefficient,weightCF,absoluteTolerance,relativeTolerance);
      }
      else if (this->_geomElemType == QUAD) {
        std::vector<std::vector<double>> coefficient;
        FELevelSetMoments::Quad9UnitCoefficientsFromNodalValues(nodalValues, coefficient);
        GetWeightFromCanonicalFE(s, coefficient, weightCF, absoluteTolerance, relativeTolerance);
      }
      else if (this->_geomElemType == TRI) {
        if(s != 0 && s != -1) throw std::invalid_argument("TRI7 nodal weights require s=0 or s=-1.");
        std::vector<double> coefficient;
        FELevelSetMoments::Tri7ReflectedCoefficientsFromNodalValues(nodalValues, coefficient);
        if(coefficient[3] == 0.0 && coefficient[4] == 0.0 && coefficient[5] == 0.0 && coefficient[6] == 0.0 && (coefficient[1] != 0.0 || coefficient[2] != 0.0)) {
          GetWeightFromCanonicalFE(s, coefficient, weightCF, absoluteTolerance, relativeTolerance);
          return;
        }
        std::vector<std::vector<double>> positive;
        std::vector<std::vector<double>> delta;
        FELevelSetMoments::IntegrateTri7NodalSelected(nodalValues, this->_qM, positive, delta, absoluteTolerance, relativeTolerance, s == -1, false, _tri7DirectionValidation, _tri7IntegrationDirection, nullptr, true);
        GetWeightFromMoments(s == 0 ? positive : delta, weightCF);
      }
      else if (this->_geomElemType == LINE) {
        double a = 0.0;
        double b = 0.0;
        double c = 0.0;
        FELevelSetMoments::Line3CoefficientsFromNodalValues(nodalValues, a, b, c);
        std::vector<double> coefficient = {a, b, c};
        GetWeightFromCanonicalFE(s, coefficient, weightCF, absoluteTolerance, relativeTolerance);
      }
      else {
        throw std::invalid_argument("FE nodal weights support LINE3, QUAD9, TRI7, HEX27, and TET15.");
      }
    }

    void SetHex27DirectionValidation(const unsigned mode) {
      if(mode > 2) throw std::invalid_argument("HEX27 direction validation must be 0, 1, or 2.");
      if(mode == _hex27DirectionValidation) return;
      _hex27DirectionValidation = mode;
      _nodalMomentCacheValid = false;
      _nodalWeightCacheValid = false;
    }

    void SetHex27IntegrationDirection(const unsigned order) {
      if(order > 6) throw std::invalid_argument("HEX27 integration direction must be between 0 and 6.");
      if(order == _hex27IntegrationDirection) return;
      _hex27IntegrationDirection = order;
      _nodalMomentCacheValid = false;
      _nodalWeightCacheValid = false;
    }

    void SetTri7DirectionValidation(const unsigned mode) {
      if(mode > 2) throw std::invalid_argument("TRI7 direction validation must be 0, 1, or 2.");
      if(mode == _tri7DirectionValidation) return;
      _tri7DirectionValidation = mode;
      _nodalMomentCacheValid = false;
      _nodalWeightCacheValid = false;
    }

    void SetTri7IntegrationDirection(const unsigned direction) {
      if(direction > 3) throw std::invalid_argument("TRI7 direction must be 0, 1, 2, or 3.");
      if(direction == _tri7IntegrationDirection) return;
      _tri7IntegrationDirection = direction;
      _nodalMomentCacheValid = false;
      _nodalWeightCacheValid = false;
    }

    void SetQuad9DirectionValidation(const unsigned mode) {
      if(mode > 2) throw std::invalid_argument("QUAD9 direction validation must be 0, 1, or 2.");
      if(mode == _quad9DirectionValidation) return;
      _quad9DirectionValidation = mode;
      _nodalMomentCacheValid = false;
      _nodalWeightCacheValid = false;
    }

    void SetQuad9IntegrationDirection(const unsigned order) {
      if(order > 2) throw std::invalid_argument("QUAD9 integration direction must be 0, 1, or 2.");
      if(order == _quad9IntegrationDirection) return;
      _quad9IntegrationDirection = order;
      _nodalMomentCacheValid = false;
      _nodalWeightCacheValid = false;
    }

    void SetTet15IntegrationDirection(const unsigned direction) {
      if(direction > 4) throw std::invalid_argument("TET15 integration direction must be between 0 and 4.");
      if(direction == _tet15IntegrationDirection) return;
      _tet15IntegrationDirection = direction;
      _nodalMomentCacheValid = false;
      _nodalWeightCacheValid = false;
    }

    void SetTet15DirectionValidation(const unsigned mode) {
      if(mode > 2) throw std::invalid_argument("TET15 direction validation must be 0, 1, or 2.");
      if(mode == _tet15DirectionValidation) return;
      _tet15DirectionValidation = mode;
      _nodalMomentCacheValid = false;
      _nodalWeightCacheValid = false;
    }

    void GetWeightsFromNodalValues(const std::vector<double>& nodalValues, std::vector<TypeIO>& positiveWeight, std::vector<TypeIO>& deltaWeight, double absoluteTolerance = 1.0e-12, double relativeTolerance = 1.0e-12, bool referenceSurface = false) {
      GetNodalWeightPair(nodalValues,positiveWeight,deltaWeight,absoluteTolerance,relativeTolerance,false,referenceSurface);
    }

    void GetWeightsFromNodalValuesByLineSlices(const std::vector<double>& nodalValues, std::vector<TypeIO>& positiveWeight, std::vector<TypeIO>& deltaWeight, double absoluteTolerance = 1.0e-12, double relativeTolerance = 1.0e-12, bool referenceSurface = false) {
      GetNodalWeightPair(nodalValues,positiveWeight,deltaWeight,absoluteTolerance,relativeTolerance,true,referenceSurface);
    }

    // For affine maps, include |J^-T grad G| in the moments and retain |det J| in the assembly weight.
    void GetPhysicalSurfaceWeightsFromNodalValues(const std::vector<double>& nodalValues, const std::vector<std::vector<double>>& inverseJacobian, std::vector<TypeIO>& positiveWeight, std::vector<TypeIO>& surfaceWeight, double absoluteTolerance = 1.0e-12, double relativeTolerance = 1.0e-12) {
      GetNodalWeightPair(nodalValues,positiveWeight,surfaceWeight,absoluteTolerance,relativeTolerance,false,true,&inverseJacobian);
    }

    void GetPhysicalSurfaceWeightsFromNodalValuesByLineSlices(const std::vector<double>& nodalValues, const std::vector<std::vector<double>>& inverseJacobian, std::vector<TypeIO>& positiveWeight, std::vector<TypeIO>& surfaceWeight, double absoluteTolerance = 1.0e-12, double relativeTolerance = 1.0e-12) {
      GetNodalWeightPair(nodalValues,positiveWeight,surfaceWeight,absoluteTolerance,relativeTolerance,true,true,&inverseJacobian);
    }

    // Geometry moments include |cof(J) grad G| and use the reference Gauss weight without |det J| during surface assembly.
    void GetPhysicalSurfaceWeightsFromNodalGeometry(const std::vector<double>& nodalValues, const std::vector<std::vector<double>>& coordinates, std::vector<TypeIO>& positiveWeight, std::vector<TypeIO>& surfaceWeight, double absoluteTolerance = 1.0e-12, double relativeTolerance = 1.0e-12) {
      GetNodalWeightPair(nodalValues,positiveWeight,surfaceWeight,absoluteTolerance,relativeTolerance,false,true,nullptr,&coordinates);
    }

    void GetPhysicalSurfaceWeightsFromNodalGeometryByLineSlices(const std::vector<double>& nodalValues, const std::vector<std::vector<double>>& coordinates, std::vector<TypeIO>& positiveWeight, std::vector<TypeIO>& surfaceWeight, double absoluteTolerance = 1.0e-12, double relativeTolerance = 1.0e-12) {
      GetNodalWeightPair(nodalValues,positiveWeight,surfaceWeight,absoluteTolerance,relativeTolerance,true,true,nullptr,&coordinates);
    }

    void UpdateQuadratureRule(const unsigned &qM) {
      const bool retainNodalMoments = _nodalMomentCacheValid;
      CutFemWeight<TypeIO, TypeA>::UpdateQuadratureRule(qM);
      _fPair.resize(this->_L, 2);
      _nodalMomentCacheValid = retainNodalMoments;
      _nodalWeightCacheValid = false;
    }

    unsigned GetCutFEMQuadratureOrder() {
      return CutFemWeight<TypeIO, TypeA>::GetCutFEMQuadratureOrder();
    }

    unsigned GetGaussQuadratureOrder() {
      return CutFemWeight<TypeIO, TypeA>::GetGaussQuadratureOrder();
    }

    unsigned GetDimension() {
      return CutFemWeight<TypeIO, TypeA>::GetDimension();
    }

    unsigned GetGaussQuadraturePointNumber() {
      return CutFemWeight<TypeIO, TypeA>::GetGaussQuadraturePointNumber();
    }

    const double* GetGaussWeightPointer() {
      return CutFemWeight<TypeIO, TypeA>::GetGaussWeightPointer();
    }

    const double* GetGaussCoordinatePointer(const unsigned &k) {
      return CutFemWeight<TypeIO, TypeA>::GetGaussCoordinatePointer(k);
    }

    const unsigned& GetCounter() {
      return CutFemWeight<TypeIO, TypeA>::GetCounter();
    }

    void ClearMap() {
      CutFemWeight<TypeIO, TypeA>::ClearMap();
      _nodalMomentCacheValid = false;
      _nodalWeightCacheValid = false;
    }

  private:
    void GetNodalWeightPair(const std::vector<double>& nodalValues, std::vector<TypeIO>& positiveWeight, std::vector<TypeIO>& deltaWeight, double absoluteTolerance, double relativeTolerance, bool useLineSlices, bool referenceSurface, const std::vector<std::vector<double>>* inverseJacobian = nullptr, const std::vector<std::vector<double>>* geometryCoordinates = nullptr) {
      if(&positiveWeight == &deltaWeight) {
        throw std::invalid_argument("Heaviside and delta weights require distinct output vectors.");
      }
      if(!(absoluteTolerance > 0.) || !(relativeTolerance >= 0.) || !std::isfinite(absoluteTolerance) || !std::isfinite(relativeTolerance)) {
        throw std::invalid_argument("Invalid nodal integration tolerance.");
      }
      if((inverseJacobian != nullptr || geometryCoordinates != nullptr) && !referenceSurface) throw std::invalid_argument("Surface geometry requires surface moments.");
      if(inverseJacobian != nullptr && geometryCoordinates != nullptr) throw std::invalid_argument("Specify a surface inverse Jacobian or geometry nodes.");
      if(inverseJacobian != nullptr) FELevelSetMoments::ValidateSurfaceInverseJacobian(*inverseJacobian,this->_dim);
      if(geometryCoordinates != nullptr) FELevelSetMoments::ValidateSurfaceGeometryCoordinates(*geometryCoordinates,this->_dim,this->_geomElemType == LINE ? 3U : (this->_geomElemType == TRI ? 7U : (this->_geomElemType == QUAD ? 9U : (this->_geomElemType == TET ? 15U : 27U))));
      useLineSlices = useLineSlices && (this->_geomElemType == LINE || this->_geomElemType == QUAD || this->_geomElemType == TRI || this->_geomElemType == HEX || this->_geomElemType == TET);
      const bool sameInverseJacobian = inverseJacobian == nullptr ? _nodalInverseJacobian.empty() : _nodalInverseJacobian == *inverseJacobian;
      const bool sameCoordinates = geometryCoordinates == nullptr ? _nodalGeometryCoordinates.empty() : _nodalGeometryCoordinates == *geometryCoordinates;
      const bool sameGeometry = sameInverseJacobian && sameCoordinates;
      bool sameField = sameGeometry && referenceSurface == _nodalMomentUsesReferenceSurface && useLineSlices == _nodalMomentUsesLineSlices && _nodalMomentCacheValid && nodalValues == _nodalValues && absoluteTolerance == _nodalAbsoluteTolerance && relativeTolerance == _nodalRelativeTolerance;
      if(sameField && _nodalWeightCacheValid) {
        positiveWeight = _nodalPositiveWeight;
        deltaWeight = _nodalDeltaWeight;
        return;
      }
      _nodalWeightCacheValid = false;
      if(sameField && this->_qM <= _nodalMomentDegree) {
        _fPair = _nodalMomentCache.topRows(this->_L);
      }
      else {
        IntegrateNodalMomentPair(nodalValues, absoluteTolerance, relativeTolerance, useLineSlices, referenceSurface, inverseJacobian, geometryCoordinates);
        _nodalMomentCache = _fPair;
        _nodalValues = nodalValues;
        _nodalMomentDegree = this->_qM;
        _nodalAbsoluteTolerance = absoluteTolerance;
        _nodalRelativeTolerance = relativeTolerance;
        _nodalMomentUsesLineSlices = useLineSlices;
        _nodalMomentUsesReferenceSurface = referenceSurface;
        if(inverseJacobian == nullptr) _nodalInverseJacobian.clear();
        else _nodalInverseJacobian = *inverseJacobian;
        if(geometryCoordinates == nullptr) _nodalGeometryCoordinates.clear();
        else _nodalGeometryCoordinates = *geometryCoordinates;
        _nodalMomentCacheValid = true;
      }
      ++this->_cntCall;
      EvaluateMomentWeightPair(positiveWeight, deltaWeight);
      _nodalPositiveWeight = positiveWeight;
      _nodalDeltaWeight = deltaWeight;
      _nodalWeightCacheValid = true;
    }

    void StoreNodalMomentPair(const std::vector<std::vector<std::vector<double>>>& positive, const std::vector<std::vector<std::vector<double>>>& delta) {
      unsigned count = 0;
      for(unsigned q = 0; q <= this->_qM; ++q) {
        for(int ii = q; ii >= 0; --ii) {
          for(int jj = q-static_cast<unsigned>(ii); jj >= 0; --jj) {
            unsigned i = static_cast<unsigned>(ii), j = static_cast<unsigned>(jj), k = q-i-j;
            _fPair(count,0) = static_cast<TypeA>(positive[i][j][k]);
            _fPair(count,1) = static_cast<TypeA>(delta[i][j][k]);
            ++count;
          }
        }
      }
    }

    void StoreNodalMomentPair(const std::vector<std::vector<double>>& positive, const std::vector<std::vector<double>>& delta) {
      unsigned count = 0;
      for(unsigned q = 0; q <= this->_qM; ++q) {
        for(unsigned j = 0; j <= q; ++j) {
          unsigned i = q - j;
          _fPair(count, 0) = static_cast<TypeA>(positive[i][j]);
          _fPair(count, 1) = static_cast<TypeA>(delta[i][j]);
          ++count;
        }
      }
    }

    void IntegratePlanarMomentPair(const std::vector<TypeIO>& normal, const TypeIO& offset) {
      std::vector<TypeA> a(this->_dim);
      for(unsigned d = 0; d < this->_dim; ++d) a[d] = static_cast<TypeA>(normal[d]);
      TypeA constant = static_cast<TypeA>(offset);
      if(this->_dim == 3) {
        unsigned count = 0;
        for(unsigned q = 0; q <= this->_qM; ++q) {
          for(int ii = q; ii >= 0; --ii) {
            for(int jj = q-static_cast<unsigned>(ii); jj >= 0; --jj) {
              unsigned i = static_cast<unsigned>(ii), j = static_cast<unsigned>(jj);
              std::vector<unsigned> powers = {i,j,q-i-j};
              _fPair(count,0) = (*this->_obj)(0,powers,a,constant);
              _fPair(count,1) = (*this->_obj)(-1,powers,a,constant);
              ++count;
            }
          }
        }
      }
      else if(this->_dim == 2) {
        unsigned count = 0;
        for(unsigned q = 0; q <= this->_qM; ++q) {
          for(unsigned j = 0; j <= q; ++j) {
            std::vector<unsigned> powers = {q - j, j};
            _fPair(count, 0) = (*this->_obj)(0, powers, a, constant);
            _fPair(count, 1) = (*this->_obj)(-1, powers, a, constant);
            ++count;
          }
        }
      }
      else {
        for(unsigned i = 0; i <= this->_qM; ++i) {
          std::vector<unsigned> powers(1, i);
          _fPair(i, 0) = (*this->_obj)(0, powers, a, constant);
          _fPair(i, 1) = (*this->_obj)(-1, powers, a, constant);
        }
      }
    }

    void ScalePlanarReferenceSurfaceMoments(const std::vector<TypeIO>& normal, const std::vector<std::vector<double>>* inverseJacobian = nullptr) {
      TypeA squaredNorm = 0;
      if(inverseJacobian == nullptr) {
        for(unsigned d = 0; d < this->_dim; ++d) squaredNorm += static_cast<TypeA>(normal[d])*static_cast<TypeA>(normal[d]);
      }
      else {
        for(unsigned d = 0; d < this->_dim; ++d) {
          TypeA component = 0;
          for(unsigned j = 0; j < this->_dim; ++j) component += static_cast<TypeA>((*inverseJacobian)[j][d])*static_cast<TypeA>(normal[j]);
          squaredNorm += component*component;
        }
      }
      TypeA norm = Eigen::numext::sqrt(squaredNorm);
      for(unsigned i = 0; i < this->_L; ++i) _fPair(i,1) *= norm;
    }

    void IntegrateLineReferenceSurfaceMomentPair(double a, double b, double c, bool zeroLeft, bool zeroRight) {
      long double scale = std::max(std::abs(static_cast<long double>(a)),std::max(std::abs(static_cast<long double>(b)),std::abs(static_cast<long double>(c))));
      if(scale == 0.0L) throw std::domain_error("The delta of an identically zero field is undefined.");
      long double quadratic = static_cast<long double>(a)/scale;
      long double linear = static_cast<long double>(b)/scale;
      long double constant = static_cast<long double>(c)/scale;
      std::vector<long double> positive;
      std::vector<long double> delta;
      FELevelSetMoments::Integrate1DInternal(quadratic,linear,constant,this->_qM,0.0L,1.0L,false,false,positive,delta);
      std::vector<long double> roots;
      if(zeroLeft && zeroRight) roots = {0.0L,1.0L};
      else if(zeroLeft) {
        roots.push_back(0.0L);
        if(quadratic != 0.0L) roots.push_back(-linear/quadratic);
      }
      else if(zeroRight) {
        roots.push_back(1.0L);
        if(quadratic != 0.0L) roots.push_back(constant/quadratic);
      }
      else roots = FELevelSetMoments::QuadraticRoots(quadratic,linear,constant);
      for(unsigned rootIndex = 0; rootIndex < roots.size(); ++rootIndex) {
        long double point = roots[rootIndex];
        if(point < 0.0L || point > 1.0L) continue;
        long double derivative = zeroRight && point == 1.0L ? quadratic-constant : 2.0L*quadratic*point+linear;
        if(derivative == 0.0L) throw std::domain_error("A repeated 1D zero has no regular surface moment.");
        long double weight = point == 0.0L || point == 1.0L ? 0.25L : 0.5L;
        long double power = 1.0L;
        for(unsigned i = 0; i <= this->_qM; ++i) {
          delta[i] += weight*power;
          power *= point;
        }
      }
      for(unsigned i = 0; i <= this->_qM; ++i) {
        _fPair(i,0) = static_cast<TypeA>(static_cast<double>(positive[i]));
        _fPair(i,1) = static_cast<TypeA>(delta[i]);
      }
    }

    void IntegrateNodalMomentPair(const std::vector<double>& nodalValues, double absoluteTolerance, double relativeTolerance, bool useLineSlices, bool referenceSurface, const std::vector<std::vector<double>>* inverseJacobian = nullptr, const std::vector<std::vector<double>>* geometryCoordinates = nullptr) {
      if(this->_geomElemType == TET) {
        std::vector<std::vector<std::vector<double>>> coefficient;
        FELevelSetMoments::Tet15ReferenceCoefficientsFromNodalValues(nodalValues,coefficient);
        bool planar = true;
        for(unsigned i = 0; i < 5; ++i) {
          for(unsigned j = 0; j < 5; ++j) {
            for(unsigned k = 0; k < 5; ++k) if(i+j+k > 1 && coefficient[i][j][k] != 0.0) planar = false;
          }
        }
        if(geometryCoordinates == nullptr && planar && (coefficient[1][0][0] != 0.0 || coefficient[0][1][0] != 0.0 || coefficient[0][0][1] != 0.0)) {
          std::vector<TypeIO> normal = {static_cast<TypeIO>(coefficient[1][0][0]),static_cast<TypeIO>(coefficient[0][1][0]),static_cast<TypeIO>(coefficient[0][0][1])};
          this->_obj->SetBaseType(1);
          IntegratePlanarMomentPair(normal,static_cast<TypeIO>(coefficient[0][0][0]));
          if(referenceSurface) ScalePlanarReferenceSurfaceMoments(normal,inverseJacobian);
          this->_obj->SetBaseType(this->_tetBaseType);
        }
        else {
          std::vector<std::vector<std::vector<double>>> positive;
          std::vector<std::vector<std::vector<double>>> delta;
          FELevelSetMoments::IntegrateTet15NodalByTriangleSlices(nodalValues,this->_qM,positive,delta,absoluteTolerance,relativeTolerance,true,24,_tet15IntegrationDirection,200000,nullptr,_tet15DirectionValidation,nullptr,referenceSurface,inverseJacobian,geometryCoordinates);
          StoreNodalMomentPair(positive,delta);
        }
      }
      else if(this->_geomElemType == HEX) {
        std::vector<std::vector<std::vector<double>>> coefficient;
        FELevelSetMoments::Hex27UnitCoefficientsFromNodalValues(nodalValues,coefficient);
        bool planar = true;
        for(unsigned i = 0; i < 3; ++i) {
          for(unsigned j = 0; j < 3; ++j) {
            for(unsigned k = 0; k < 3; ++k) if(i+j+k > 1 && coefficient[i][j][k] != 0.) planar = false;
          }
        }
        if(geometryCoordinates == nullptr && planar && (coefficient[1][0][0] != 0. || coefficient[0][1][0] != 0. || coefficient[0][0][1] != 0.)) {
          std::vector<TypeIO> normal = {static_cast<TypeIO>(0.5L*coefficient[1][0][0]),static_cast<TypeIO>(0.5L*coefficient[0][1][0]),static_cast<TypeIO>(0.5L*coefficient[0][0][1])};
          IntegratePlanarMomentPair(normal,static_cast<TypeIO>(static_cast<long double>(coefficient[0][0][0])+0.5L*coefficient[1][0][0]+0.5L*coefficient[0][1][0]+0.5L*coefficient[0][0][1]));
          if(referenceSurface) ScalePlanarReferenceSurfaceMoments(normal,inverseJacobian);
          return;
        }
        std::vector<std::vector<std::vector<double>>> positive, delta;
        FELevelSetMoments::IntegrateHex27NodalSelected(nodalValues,this->_qM,positive,delta,absoluteTolerance,relativeTolerance,true,_hex27DirectionValidation,_hex27IntegrationDirection,nullptr,!useLineSlices,nullptr,referenceSurface,inverseJacobian,geometryCoordinates);
        StoreNodalMomentPair(positive,delta);
      }
      else if(this->_geomElemType == QUAD) {
        std::vector<std::vector<double>> coefficient;
        FELevelSetMoments::Quad9UnitCoefficientsFromNodalValues(nodalValues, coefficient);
        bool planar = true;
        for(unsigned i = 0; i < 3; ++i) {
          for(unsigned j = 0; j < 3; ++j) {
            if(i + j > 1 && coefficient[i][j] != 0.) planar = false;
          }
        }
        if(geometryCoordinates == nullptr && planar && (coefficient[1][0] != 0. || coefficient[0][1] != 0.)) {
          std::vector<TypeIO> normal = {static_cast<TypeIO>(0.5L * coefficient[1][0]), static_cast<TypeIO>(0.5L * coefficient[0][1])};
          IntegratePlanarMomentPair(normal, static_cast<TypeIO>(static_cast<long double>(coefficient[0][0]) + 0.5L * coefficient[1][0] + 0.5L * coefficient[0][1]));
          if(referenceSurface) ScalePlanarReferenceSurfaceMoments(normal,inverseJacobian);
          return;
        }
        std::vector<std::vector<double>> positive;
        std::vector<std::vector<double>> delta;
        FELevelSetMoments::IntegrateQuad9NodalSelected(nodalValues, this->_qM, positive, delta, absoluteTolerance, relativeTolerance, true, false, _quad9DirectionValidation, _quad9IntegrationDirection, nullptr, !useLineSlices, referenceSurface, inverseJacobian, geometryCoordinates);
        StoreNodalMomentPair(positive, delta);
      }
      else if(this->_geomElemType == TRI) {
        std::vector<double> coefficient;
        FELevelSetMoments::Tri7ReflectedCoefficientsFromNodalValues(nodalValues, coefficient);
        if(geometryCoordinates == nullptr && coefficient[3] == 0. && coefficient[4] == 0. && coefficient[5] == 0. && coefficient[6] == 0. && (coefficient[1] != 0. || coefficient[2] != 0.)) {
          std::vector<TypeIO> normal = {static_cast<TypeIO>(-coefficient[1]), static_cast<TypeIO>(coefficient[2])};
          IntegratePlanarMomentPair(normal, static_cast<TypeIO>(static_cast<long double>(coefficient[0]) + coefficient[1]));
          if(referenceSurface) ScalePlanarReferenceSurfaceMoments(normal,inverseJacobian);
          return;
        }
        std::vector<std::vector<double>> positive;
        std::vector<std::vector<double>> delta;
        FELevelSetMoments::IntegrateTri7NodalSelected(nodalValues, this->_qM, positive, delta, absoluteTolerance, relativeTolerance, true, false, _tri7DirectionValidation, _tri7IntegrationDirection, nullptr, !useLineSlices, referenceSurface, inverseJacobian, geometryCoordinates);
        StoreNodalMomentPair(positive, delta);
      }
      else if(this->_geomElemType == LINE) {
        double a = 0.;
        double b = 0.;
        double c = 0.;
        FELevelSetMoments::Line3CoefficientsFromNodalValues(nodalValues, a, b, c);
        if(a == 0. && b != 0.) {
          std::vector<TypeIO> normal(1, static_cast<TypeIO>(0.5L * b));
          IntegratePlanarMomentPair(normal, static_cast<TypeIO>(static_cast<long double>(c) + 0.5L * b));
          if(referenceSurface) ScalePlanarReferenceSurfaceMoments(normal,inverseJacobian);
          return;
        }
        if(referenceSurface) {
          IntegrateLineReferenceSurfaceMomentPair(a,b,c,nodalValues[0] == 0.0,nodalValues[1] == 0.0);
          if(inverseJacobian != nullptr) {
            const TypeA factor = static_cast<TypeA>(std::abs((*inverseJacobian)[0][0]));
            for(unsigned i = 0; i < this->_L; ++i) _fPair(i,1) *= factor;
          }
          return;
        }
        std::vector<double> positive;
        std::vector<double> delta;
        if(useLineSlices) FELevelSetMoments::IntegrateLine3NodalByLineSlices(nodalValues,this->_qM,positive,delta,absoluteTolerance,relativeTolerance);
        else FELevelSetMoments::Integrate1D(a, b, c, this->_qM, positive, delta);
        for(unsigned i = 0; i <= this->_qM; ++i) {
          _fPair(i, 0) = static_cast<TypeA>(positive[i]);
          _fPair(i, 1) = static_cast<TypeA>(delta[i]);
        }
      }
      else {
        throw std::invalid_argument("Paired nodal weights support LINE3, QUAD9, TRI7, HEX27, and TET15.");
      }
    }


    void EvaluateMomentWeightPair(std::vector<TypeIO>& positiveWeight, std::vector<TypeIO>& deltaWeight) {
      positiveWeight.resize(this->_gn);
      deltaWeight.resize(this->_gn);
      for(unsigned ig = 0; ig < this->_gn; ++ig) {
        const TypeA* basisRow = this->_momentToWeight[0].data()+ig*this->_L;
        TypeA positive = 0.;
        TypeA delta = 0.;
        for(unsigned i = 0; i < this->_L; ++i) {
          positive += basisRow[i]*_fPair(i,0);
          delta += basisRow[i]*_fPair(i,1);
        }
        positiveWeight[ig] = static_cast<TypeIO>(positive);
        deltaWeight[ig] = static_cast<TypeIO>(delta);
      }
    }

    Eigen::Matrix<TypeA, Eigen::Dynamic, 2> _fPair;
    Eigen::Matrix<TypeA, Eigen::Dynamic, 2> _nodalMomentCache;
    unsigned _hex27DirectionValidation = 1;
    unsigned _hex27IntegrationDirection = 0;
    unsigned _tet15DirectionValidation = 1;
    unsigned _tet15IntegrationDirection = 0;
    unsigned _quad9DirectionValidation = 1;
    unsigned _quad9IntegrationDirection = 0;
    unsigned _tri7DirectionValidation = 1;
    unsigned _tri7IntegrationDirection = 0;
    bool _nodalMomentUsesLineSlices = false;
    bool _nodalMomentUsesReferenceSurface = false;
    bool _nodalMomentCacheValid = false;
    bool _nodalWeightCacheValid = false;
    unsigned _nodalMomentDegree = 0;
    double _nodalAbsoluteTolerance = 0.;
    double _nodalRelativeTolerance = 0.;
    std::vector<double> _nodalValues;
    std::vector<std::vector<double>> _nodalInverseJacobian;
    std::vector<std::vector<double>> _nodalGeometryCoordinates;
    std::vector<TypeIO> _nodalPositiveWeight;
    std::vector<TypeIO> _nodalDeltaWeight;
};

template <class TypeIO, class TypeA>
void CutFemWeightLS<TypeIO, TypeA>::operator()(const int &s, const std::vector<TypeA> &a, const TypeA &d, const bool &xSquared, std::vector<TypeIO> &weightCF) {
  if(this->_geomElemType != QUAD) throw std::invalid_argument("Quadratic CutFemWeightLS requires QUAD");
  if(a.size() != 3) throw std::invalid_argument("Quadratic CutFemWeightLS requires a={k,b,c}");
  if(s < -1) throw std::invalid_argument("Quadratic CutFemWeightLS requires s >= -1");

  this->_cntCall++;

  unsigned count = 0;
  for(unsigned q = 0; q <= this->_qM; q++) {
    for(unsigned j = 0; j <= q; j++) {
      unsigned i = q - j;
      this->_f[count] = SquareParabolaMoment<TypeA>(s, {i, j}, a, d, xSquared);
      count++;
    }
  }

  this->EvaluateMomentWeights(weightCF);
}

template <class TypeIO, class TypeA>
void CutFemWeightLS<TypeIO, TypeA>::operator()(const int &s, const std::vector<TypeA> &aX, const TypeA &dX, const std::vector<TypeA> &aY, const TypeA &dY, const bool &xSquared, std::vector<TypeIO> &weightCF) {
  if(xSquared) operator()(s, aX, dX, true, weightCF);
  else operator()(s, aY, dY, false, weightCF);
}

#endif
