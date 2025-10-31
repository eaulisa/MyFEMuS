 
#pragma once
#include <cmath>

class Mollifier {
public:
  Mollifier() = default;
  explicit Mollifier(double eps) { SetConstants(eps); }

  inline void SetConstants(double eps) noexcept {
    _eps = eps;

    // precompute inverse powers (faster than std::pow)
    const double inv  = 1.0 / eps;
    const double inv2 = inv * inv;
    const double inv3 = inv2 * inv;
    const double inv5 = inv3 * inv2;   // inv^5
    const double inv7 = inv5 * inv2;   // inv^7
    const double inv9 = inv7 * inv2;   // inv^9

    _a0 = 0.5;                 // 128/256
    _a1 = 1.23046875  * inv;   // 315/256 * eps^-1
    _a3 = -1.640625   * inv3;  // 420/256 * eps^-3
    _a5 = 1.4765625   * inv5;  // 378/256 * eps^-5
    _a7 = -0.703125   * inv7;  // 180/256 * eps^-7
    _a9 = 0.13671875  * inv9;  //  35/256 * eps^-9
  }

  // Smooth, piecewise unit step centered at 0 with width 2*eps
  inline double GetSmoothStepFunction(double dg1) const noexcept {
    if (dg1 < -_eps) return 0.0;
    if (dg1 >  _eps) return 1.0;

    const double dg2 = dg1 * dg1;
    return _a0 + dg1 * (_a1 + dg2 * (_a3 + dg2 * (_a5 + dg2 * (_a7 + dg2 * _a9))));
  }

  // Heaviside on [-1,1] via stretched smoothstep; exact 0/1 outside.
  inline double Sigmoid(double dg1) const noexcept {
    if (dg1 < -_eps) return -1.0;
    if (dg1 >  _eps) return 1.0;
    return -1. + 2. * GetSmoothStepFunction(dg1);
  }

  inline double eps() const noexcept { return _eps; }

private:
  double _eps{1.0};
  double _a0{0.5}, _a1{0.0}, _a3{0.0}, _a5{0.0}, _a7{0.0}, _a9{0.0};
};
