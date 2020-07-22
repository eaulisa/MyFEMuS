/*=========================================================================

 Program: FEMUS
 Module: ElemType
 Authors: Eugenio Aulisa, Giorgio Bornia

 Copyright (c) FEMTTU
 All rights reserved.

 This software is distributed WITHOUT ANY WARRANTY; without even
 the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

//----------------------------------------------------------------------------
// includes :
//----------------------------------------------------------------------------

#include "ElemType.hpp"
#include "GaussPoints.hpp"
#include "FETypeEnum.hpp"
#include "Elem.hpp"
#include "NumericVector.hpp"

using std::cout;
using std::endl;


#define  PHIFACE_ONLY_FOR_LAGRANGIAN_FAMILIES  1



namespace femus {

//BEGIN Base Class specialization for template and inheritance.

  void elem_type::GetGaussQuantities(const vector < vector < adept::adouble > >& vt, const unsigned& ig,
                                     adept::adouble &weight,
                                     boost::optional < vector < adept::adouble > & > gradphi,
                                     boost::optional < vector < adept::adouble > & > nablaphi) const {

//     GetGaussQuantities_type(vt, ig, weight, gradphi, nablaphi);

    if(_dim1) {
      static_cast<const elem_type_1D&>(*this).GetGaussQuantities_type(vt, ig, weight, gradphi, nablaphi);
    }
    else if(_dim2) {
      static_cast<const elem_type_2D&>(*this).GetGaussQuantities_type(vt, ig, weight, gradphi, nablaphi);
    }
    else {
      static_cast<const elem_type_3D&>(*this).GetGaussQuantities_type(vt, ig, weight, gradphi, nablaphi);
    }
  }


  void elem_type::GetGaussQuantities(const vector < vector < double > >& vt, const unsigned& ig,
                                     double& weight,
                                     boost::optional < vector < double > & > gradphi,
                                     boost::optional < vector < double > & > nablaphi) const {

    if(_dim1) {
      static_cast<const elem_type_1D&>(*this).GetGaussQuantities_type(vt, ig, weight, gradphi, nablaphi);
    }
    else if(_dim2) {
      static_cast<const elem_type_2D&>(*this).GetGaussQuantities_type(vt, ig, weight, gradphi, nablaphi);
    }
    else {
      static_cast<const elem_type_3D&>(*this).GetGaussQuantities_type(vt, ig, weight, gradphi, nablaphi);
    }

  }

  //END Base Class specialization for template and inheritance.


  //BEGIN Derived Class specialization for template and inheritance.
  template <class type>
  void elem_type_1D::GetGaussQuantities_type(const vector < vector < type > >& vt, const unsigned& ig,
                                             type& weight,
                                             boost::optional < vector < type >& > gradphi,
                                             boost::optional < vector < type > & > nablaphi) const {

    gradphi->resize(_nc * 1);
    if(nablaphi) nablaphi->resize(_nc * 1);

    type Jac = 0.;
    type JacI;

    const double* dxi = _dphidxi[ig];

    for(int inode = 0; inode < _nc; inode++, dxi++) {
      Jac += (*dxi) * vt[0][inode];
    }

    weight = Jac * _gauss.GetGaussWeightsPointer()[ig];

    JacI = 1 / Jac;

    dxi = _dphidxi[ig];
    const double* dxi2 = _d2phidxi2[ig];

    for(int inode = 0; inode < _nc; inode++, dxi++, dxi2++) {
      (*gradphi)[inode] = (*dxi) * JacI;
      if(nablaphi)(*nablaphi)[inode] = (*dxi2) * JacI * JacI;
    }

  }

  ///////////////////////////////////////////////////////

  template <class type>
  void elem_type_2D::GetGaussQuantities_type(const vector < vector < type > >& vt, const unsigned& ig,
                                             type & weight,
                                             boost::optional < vector < type >& > gradphi,
                                             boost::optional < vector < type > & > nablaphi) const {

    type Jac[2][2] = {{0, 0}, {0, 0}};

    const double* dxi = _dphidxi[ig];
    const double* deta = _dphideta[ig];
    typename std::vector< type >::const_iterator vt0i = vt[0].begin();
    typename std::vector< type >::const_iterator vt1i = vt[1].begin();
    unsigned inode;

    // CAREFUL!!! Do not use vt[0].end() to close this loop since vt[0].size() >= _nc
    for(inode = 0; inode < _nc; inode++, dxi++, deta++, vt0i++, vt1i++) {
      Jac[0][0] += (*dxi) * (*vt0i);
      Jac[0][1] += (*dxi) * (*vt1i);
      Jac[1][0] += (*deta) * (*vt0i);
      Jac[1][1] += (*deta) * (*vt1i);
    }

    type det = (Jac[0][0] * Jac[1][1] - Jac[0][1] * Jac[1][0]);
    weight = det * _gauss.GetGaussWeightsPointer()[ig];

    if(gradphi) {

      type JacI[2][2];
      JacI[0][0] = Jac[1][1] / det;
      JacI[0][1] = -Jac[0][1] / det;
      JacI[1][0] = -Jac[1][0] / det;
      JacI[1][1] = Jac[0][0] / det;

      gradphi->resize(_nc * 2);
      dxi = _dphidxi[ig];
      deta = _dphideta[ig];

      for(int inode = 0; inode < _nc; inode++, dxi++, deta++) {
        (*gradphi)[2 * inode + 0] = (*dxi) * JacI[0][0] + (*deta) * JacI[0][1];
        (*gradphi)[2 * inode + 1] = (*dxi) * JacI[1][0] + (*deta) * JacI[1][1];
      }

      if(nablaphi) {
        nablaphi->resize(_nc * 3);
        const double* dxi2 = _d2phidxi2[ig];
        const double* deta2 = _d2phideta2[ig];
        const double* dxideta = _d2phidxideta[ig];
        for(int inode = 0; inode < _nc; inode++, dxi2++, deta2++, dxideta++) {
          if(nablaphi) {
            (*nablaphi)[3 * inode + 0] =
              ((*dxi2)   * JacI[0][0] + (*dxideta) * JacI[0][1]) * JacI[0][0] +
              ((*dxideta) * JacI[0][0] + (*deta2)  * JacI[0][1]) * JacI[0][1];
            (*nablaphi)[3 * inode + 1] =
              ((*dxi2)   * JacI[1][0] + (*dxideta) * JacI[1][1]) * JacI[1][0] +
              ((*dxideta) * JacI[1][0] + (*deta2)  * JacI[1][1]) * JacI[1][1];
            (*nablaphi)[3 * inode + 2] =
              ((*dxi2)   * JacI[0][0] + (*dxideta) * JacI[0][1]) * JacI[1][0] +
              ((*dxideta) * JacI[0][0] + (*deta2)  * JacI[0][1]) * JacI[1][1];
          }
        }
      }
    }
  }

  ///////////////////////////////////////////////////////

  template <class type>
  void elem_type_3D::GetGaussQuantities_type(const vector < vector < type > >& vt, const unsigned& ig,
                                             type& weight,
                                             boost::optional < vector < type >& > gradphi,
                                             boost::optional < vector < type > & > nablaphi) const {

    gradphi->resize(_nc * 3);
    if(nablaphi) nablaphi->resize(_nc * 6);


    type Jac[3][3] = {{0., 0., 0.}, {0., 0., 0.}, {0., 0., 0.}};
    type JacI[3][3];

    const double* dxi = _dphidxi[ig];
    const double* deta = _dphideta[ig];
    const double* dzeta = _dphidzeta[ig];

    for(int inode = 0; inode < _nc; inode++, dxi++, deta++, dzeta++) {
      Jac[0][0] += (*dxi) * vt[0][inode];
      Jac[0][1] += (*dxi) * vt[1][inode];
      Jac[0][2] += (*dxi) * vt[2][inode];
      Jac[1][0] += (*deta) * vt[0][inode];
      Jac[1][1] += (*deta) * vt[1][inode];
      Jac[1][2] += (*deta) * vt[2][inode];
      Jac[2][0] += (*dzeta) * vt[0][inode];
      Jac[2][1] += (*dzeta) * vt[1][inode];
      Jac[2][2] += (*dzeta) * vt[2][inode];
    }

    type det = (Jac[0][0] * (Jac[1][1] * Jac[2][2] - Jac[1][2] * Jac[2][1]) +
                Jac[0][1] * (Jac[1][2] * Jac[2][0] - Jac[1][0] * Jac[2][2]) +
                Jac[0][2] * (Jac[1][0] * Jac[2][1] - Jac[1][1] * Jac[2][0]));

    JacI[0][0] = (-Jac[1][2] * Jac[2][1] + Jac[1][1] * Jac[2][2]) / det;
    JacI[0][1] = (Jac[0][2] * Jac[2][1] - Jac[0][1] * Jac[2][2]) / det;
    JacI[0][2] = (-Jac[0][2] * Jac[1][1] + Jac[0][1] * Jac[1][2]) / det;
    JacI[1][0] = (Jac[1][2] * Jac[2][0] - Jac[1][0] * Jac[2][2]) / det;
    JacI[1][1] = (-Jac[0][2] * Jac[2][0] + Jac[0][0] * Jac[2][2]) / det;
    JacI[1][2] = (Jac[0][2] * Jac[1][0] - Jac[0][0] * Jac[1][2]) / det;
    JacI[2][0] = (-Jac[1][1] * Jac[2][0] + Jac[1][0] * Jac[2][1]) / det;
    JacI[2][1] = (Jac[0][1] * Jac[2][0] - Jac[0][0] * Jac[2][1]) / det;
    JacI[2][2] = (-Jac[0][1] * Jac[1][0] + Jac[0][0] * Jac[1][1]) / det;

    weight = det * _gauss.GetGaussWeightsPointer()[ig];

    dxi = _dphidxi[ig];
    deta = _dphideta[ig];
    dzeta = _dphidzeta[ig];

    const double* dxi2 = _d2phidxi2[ig];
    const double* deta2 = _d2phideta2[ig];
    const double* dzeta2 = _d2phidzeta2[ig];
    const double* dxideta = _d2phidxideta[ig];
    const double* detadzeta = _d2phidetadzeta[ig];
    const double* dzetadxi = _d2phidzetadxi[ig];

    for(int inode = 0; inode < _nc; inode++, dxi++, deta++, dzeta++, dxi2++, deta2++, dzeta2++, dxideta++, detadzeta++, dzetadxi++) {

      (*gradphi)[3 * inode + 0] = (*dxi) * JacI[0][0] + (*deta) * JacI[0][1] + (*dzeta) * JacI[0][2];
      (*gradphi)[3 * inode + 1] = (*dxi) * JacI[1][0] + (*deta) * JacI[1][1] + (*dzeta) * JacI[1][2];
      (*gradphi)[3 * inode + 2] = (*dxi) * JacI[2][0] + (*deta) * JacI[2][1] + (*dzeta) * JacI[2][2];

      if(nablaphi) {
        (*nablaphi)[6 * inode + 0] =
          ((*dxi2)    * JacI[0][0] + (*dxideta)  * JacI[0][1] + (*dzetadxi) * JacI[0][2]) * JacI[0][0] +
          ((*dxideta) * JacI[0][0] + (*deta2)    * JacI[0][1] + (*detadzeta) * JacI[0][2]) * JacI[0][1] +
          ((*dzetadxi) * JacI[0][0] + (*detadzeta) * JacI[0][1] + (*dzeta2)   * JacI[0][2]) * JacI[0][2];
        (*nablaphi)[6 * inode + 1] =
          ((*dxi2)    * JacI[1][0] + (*dxideta)  * JacI[1][1] + (*dzetadxi) * JacI[1][2]) * JacI[1][0] +
          ((*dxideta) * JacI[1][0] + (*deta2)    * JacI[1][1] + (*detadzeta) * JacI[1][2]) * JacI[1][1] +
          ((*dzetadxi) * JacI[1][0] + (*detadzeta) * JacI[1][1] + (*dzeta2)   * JacI[1][2]) * JacI[1][2];
        (*nablaphi)[6 * inode + 2] =
          ((*dxi2)    * JacI[2][0] + (*dxideta)  * JacI[2][1] + (*dzetadxi) * JacI[2][2]) * JacI[2][0] +
          ((*dxideta) * JacI[2][0] + (*deta2)    * JacI[2][1] + (*detadzeta) * JacI[2][2]) * JacI[2][1] +
          ((*dzetadxi) * JacI[2][0] + (*detadzeta) * JacI[2][1] + (*dzeta2)   * JacI[2][2]) * JacI[2][2];
        (*nablaphi)[6 * inode + 3] =
          ((*dxi2)    * JacI[0][0] + (*dxideta)  * JacI[0][1] + (*dzetadxi) * JacI[0][2]) * JacI[1][0] +
          ((*dxideta) * JacI[0][0] + (*deta2)    * JacI[0][1] + (*detadzeta) * JacI[0][2]) * JacI[1][1] +
          ((*dzetadxi) * JacI[0][0] + (*detadzeta) * JacI[0][1] + (*dzeta2)   * JacI[0][2]) * JacI[1][2];
        (*nablaphi)[6 * inode + 4] =
          ((*dxi2)    * JacI[1][0] + (*dxideta)  * JacI[1][1] + (*dzetadxi) * JacI[1][2]) * JacI[2][0] +
          ((*dxideta) * JacI[1][0] + (*deta2)    * JacI[1][1] + (*detadzeta) * JacI[1][2]) * JacI[2][1] +
          ((*dzetadxi) * JacI[1][0] + (*detadzeta) * JacI[1][1] + (*dzeta2)   * JacI[1][2]) * JacI[2][2];
        (*nablaphi)[6 * inode + 5] =
          ((*dxi2)    * JacI[2][0] + (*dxideta)  * JacI[2][1] + (*dzetadxi) * JacI[2][2]) * JacI[0][0] +
          ((*dxideta) * JacI[2][0] + (*deta2)    * JacI[2][1] + (*detadzeta) * JacI[2][2]) * JacI[0][1] +
          ((*dzetadxi) * JacI[2][0] + (*detadzeta) * JacI[2][1] + (*dzeta2)   * JacI[2][2]) * JacI[0][2];
      }
    }

  }
  //END Derived Class specialization for template and inheritance.

} //end namespace femus









