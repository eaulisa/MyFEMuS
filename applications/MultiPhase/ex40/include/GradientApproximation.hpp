#pragma once

#include <vector>
#include <cmath>
#include <stdexcept>

// -----------------------------------------------------------------------------
// Solve a 2x2 linear system
//
// [a00 a01] [x0] = [b0]
// [a10 a11] [x1]   [b1]
//
// Returns true if successful, false if nearly singular.
// -----------------------------------------------------------------------------
inline bool solveSmallSystem2x2(const double a00, const double a01,
                                const double a10, const double a11,
                                const double b0,  const double b1,
                                std::vector<double>& x)
{
  const double det = a00 * a11 - a01 * a10;
  const double scale = std::abs(a00) + std::abs(a01) + std::abs(a10) + std::abs(a11);
  const double tol = (scale > 0.0) ? 1.e-14 * scale * scale : 1.e-14;

  if(std::abs(det) <= tol) return false;

  if(x.size() != 2u) x.resize(2u);

  const double invDet = 1.0 / det;
  x[0] = ( a11 * b0 - a01 * b1) * invDet;
  x[1] = (-a10 * b0 + a00 * b1) * invDet;

  return true;
}

// -----------------------------------------------------------------------------
// Solve a 3x3 linear system
//
// A x = b
//
// Returns true if successful, false if nearly singular.
// -----------------------------------------------------------------------------
inline bool solveSmallSystem3x3(const double a00, const double a01, const double a02,
                                const double a10, const double a11, const double a12,
                                const double a20, const double a21, const double a22,
                                const double b0,  const double b1,  const double b2,
                                std::vector<double>& x)
{
  const double c00 =  a11 * a22 - a12 * a21;
  const double c01 = -a10 * a22 + a12 * a20;
  const double c02 =  a10 * a21 - a11 * a20;

  const double c10 = -a01 * a22 + a02 * a21;
  const double c11 =  a00 * a22 - a02 * a20;
  const double c12 = -a00 * a21 + a01 * a20;

  const double c20 =  a01 * a12 - a02 * a11;
  const double c21 = -a00 * a12 + a02 * a10;
  const double c22 =  a00 * a11 - a01 * a10;

  const double det = a00 * c00 + a01 * c01 + a02 * c02;

  const double scale =
      std::abs(a00) + std::abs(a01) + std::abs(a02) +
      std::abs(a10) + std::abs(a11) + std::abs(a12) +
      std::abs(a20) + std::abs(a21) + std::abs(a22);

  const double tol = (scale > 0.0) ? 1.e-14 * scale * scale * scale : 1.e-14;

  if(std::abs(det) <= tol) return false;

  if(x.size() != 3u) x.resize(3u);

  const double invDet = 1.0 / det;

  x[0] = (c00 * b0 + c10 * b1 + c20 * b2) * invDet;
  x[1] = (c01 * b0 + c11 * b1 + c21 * b2) * invDet;
  x[2] = (c02 * b0 + c12 * b1 + c22 * b2) * invDet;

  return true;
}

// -----------------------------------------------------------------------------
// Compute one constant gradient estimate on one element.
//
// INPUT:
//   xloc[d][a] : coordinate component d of local node a
//   uloc[a]    : nodal value at local node a
//
// OUTPUT:
//   grad[d]    : estimated constant gradient on the element
//
// This is an affine least-squares reconstruction.
// Works for arbitrary element types in 2D or 3D.
// No shape functions are used.
// -----------------------------------------------------------------------------
inline void computeElementGradientFromLocalData(
    const std::vector<std::vector<double>>& xloc,
    const std::vector<double>& uloc,
    std::vector<double>& grad)
{
  const std::size_t dim = xloc.size();
  if(dim != 2u && dim != 3u) {
    throw std::runtime_error("computeElementGradientFromLocalData: dim must be 2 or 3");
  }

  const std::size_t nloc = uloc.size();
  if(nloc == 0u) {
    throw std::runtime_error("computeElementGradientFromLocalData: empty uloc");
  }

  for(std::size_t d = 0; d < dim; ++d) {
    if(xloc[d].size() != nloc) {
      throw std::runtime_error("computeElementGradientFromLocalData: inconsistent xloc[d].size()");
    }
  }

  if(grad.size() != dim) grad.resize(dim);
  for(std::size_t d = 0; d < dim; ++d) grad[d] = 0.0;

  if(nloc < dim + 1u) return;

  // ---- centroid and average value
  double xc0 = 0.0, xc1 = 0.0, xc2 = 0.0;
  double ubar = 0.0;

  if(dim == 2u) {
    for(std::size_t a = 0; a < nloc; ++a) {
      xc0 += xloc[0][a];
      xc1 += xloc[1][a];
      ubar += uloc[a];
    }
    const double invN = 1.0 / static_cast<double>(nloc);
    xc0 *= invN;
    xc1 *= invN;
    ubar *= invN;

    // Build normal equations M grad = rhs
    double m00 = 0.0, m01 = 0.0, m11 = 0.0;
    double r0s = 0.0, r1s = 0.0;

    for(std::size_t a = 0; a < nloc; ++a) {
      const double dx = xloc[0][a] - xc0;
      const double dy = xloc[1][a] - xc1;
      const double du = uloc[a] - ubar;

      m00 += dx * dx;
      m01 += dx * dy;
      m11 += dy * dy;

      r0s += dx * du;
      r1s += dy * du;
    }

    if(!solveSmallSystem2x2(m00, m01,
                            m01, m11,
                            r0s, r1s,
                            grad))
    {
      grad[0] = 0.0;
      grad[1] = 0.0;
    }

    return;
  }

  // ---- 3D
  for(std::size_t a = 0; a < nloc; ++a) {
    xc0 += xloc[0][a];
    xc1 += xloc[1][a];
    xc2 += xloc[2][a];
    ubar += uloc[a];
  }
  const double invN = 1.0 / static_cast<double>(nloc);
  xc0 *= invN;
  xc1 *= invN;
  xc2 *= invN;
  ubar *= invN;

  double m00 = 0.0, m01 = 0.0, m02 = 0.0;
  double m11 = 0.0, m12 = 0.0, m22 = 0.0;
  double r0s = 0.0, r1s = 0.0, r2s = 0.0;

  for(std::size_t a = 0; a < nloc; ++a) {
    const double dx = xloc[0][a] - xc0;
    const double dy = xloc[1][a] - xc1;
    const double dz = xloc[2][a] - xc2;
    const double du = uloc[a] - ubar;

    m00 += dx * dx;
    m01 += dx * dy;
    m02 += dx * dz;
    m11 += dy * dy;
    m12 += dy * dz;
    m22 += dz * dz;

    r0s += dx * du;
    r1s += dy * du;
    r2s += dz * du;
  }

  if(!solveSmallSystem3x3(m00, m01, m02,
                          m01, m11, m12,
                          m02, m12, m22,
                          r0s, r1s, r2s,
                          grad))
  {
    grad[0] = 0.0;
    grad[1] = 0.0;
    grad[2] = 0.0;
  }
}
