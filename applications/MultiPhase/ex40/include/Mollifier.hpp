
#pragma once
#include <cmath>

class Mollifier {
  public:
    Mollifier() = default;
    explicit Mollifier(double eps) {
      // Initialize all families to the same epsilon so they're ready to use
      SetConstantsC1(eps);
      SetConstantsC2(eps);
      SetConstantsC3(eps);
      SetConstantsC4(eps);
    }

    // ---------------- C^4 (nonic) : zero value + first 4 derivatives at ±eps ----------------
    // S4(u) = 126u^5 - 420u^6 + 540u^7 - 315u^8 + 70u^9 (u in [0,1])
    // Expanded on s in [-1,1], s = dg1/eps:
    // S4(s) = 1/2 + (315/256)s + (-105/64)s^3 + (189/128)s^5 + (-45/64)s^7 + (35/256)s^9
    inline void SetConstantsC4(double eps) noexcept {
      _eps = eps;

      const double inv  = 1.0 / eps;
      const double inv2 = inv * inv;
      const double inv3 = inv2 * inv;
      const double inv5 = inv3 * inv2;   // inv^5
      const double inv7 = inv5 * inv2;   // inv^7
      const double inv9 = inv7 * inv2;   // inv^9

      // Coefficients are the expanded S4(s) coefficients, scaled by inv^k for use in dg1^k
      _a0 = 0.5;
      _a1 = 1.23046875  * inv;   // 315/256 * inv
      _a3 = -1.640625   * inv3;  // -105/64  * inv^3
      _a5 = 1.4765625   * inv5;  // 189/128  * inv^5
      _a7 = -0.703125   * inv7;  // -45/64   * inv^7
      _a9 = 0.13671875  * inv9;  // 35/256   * inv^9
    }

    // Smooth unit step on [-eps, eps]; exact 0/1 outside; C^4 at the joins
    inline double GetSmoothStepFunctionC4(double dg1) const noexcept {
      if (dg1 < -_eps) return 0.0;
      if (dg1 >  _eps) return 1.0;
      const double dg2 = dg1 * dg1;
      return _a0 + dg1 * (_a1 + dg2 * (_a3 + dg2 * (_a5 + dg2 * (_a7 + dg2 * _a9))));
    }

    // ---------------- C^2 (quintic) : zero value + first 2 derivatives at ±eps ----------------
    // S2(u) = 6u^5 - 15u^4 + 10u^3
    // Expanded on s in [-1,1]: S2(s) = 1/2 + (15/16)s + (-5/8)s^3 + (3/16)s^5
    inline void SetConstantsC2(double eps) noexcept {
      _eps = eps;

      const double inv  = 1.0 / eps;
      const double inv2 = inv * inv;
      const double inv3 = inv2 * inv;
      const double inv5 = inv3 * inv2;

      _b0 = 0.5;
      _b1 = (15.0 / 16.0) * inv;     // s term → inv
      _b3 = (-5.0  / 8.0)  * inv3;   // s^3 term → inv^3
      _b5 = (3.0   / 16.0) * inv5;   // s^5 term → inv^5
    }

    inline double GetSmoothStepFunctionC2(double dg1) const noexcept {
      if (dg1 <= -_eps) return 0.0;
      if (dg1 >=  _eps) return 1.0;
      const double dg2 = dg1 * dg1;
      return _b0 + dg1 * (_b1 + dg2 * (_b3 + dg2 * _b5));
    }

    // ---------------- C^1 (cubic) : zero value + first derivative at ±eps ----------------
    // S1(u) = 3u^2 - 2u^3
    // Expanded on s in [-1,1]: S1(s) = 1/2 + (3/4)s + (-1/4)s^3
    inline void SetConstantsC1(double eps) noexcept {
      _eps = eps;

      const double inv  = 1.0 / eps;
      const double inv2 = inv * inv;
      const double inv3 = inv2 * inv;

      _c0 = 0.5;
      _c1 = (3.0 / 4.0) * inv;      // s term → inv
      _c3 = (-1.0 / 4.0) * inv3;    // s^3 term → inv^3
    }

    inline double GetSmoothStepFunctionC1(double dg1) const noexcept {
      if (dg1 <= -_eps) return 0.0;
      if (dg1 >=  _eps) return 1.0;
      const double dg2 = dg1 * dg1;
      // Horner: c0 + dg1*(c1 + dg2*c3)
      return _c0 + dg1 * (_c1 + dg2 * _c3);
    }

    // ---------------- C^3 (septic) : zero value + first 3 derivatives at ±eps ----------------
    // S3(u) = 35u^4 - 84u^5 + 70u^6 - 20u^7
    // Expanded on s in [-1,1]: S3(s) = 1/2 + (35/32)s + (-35/32)s^3 + (21/32)s^5 + (-5/32)s^7
    inline void SetConstantsC3(double eps) noexcept {
      _eps = eps;

      const double inv  = 1.0 / eps;
      const double inv2 = inv * inv;
      const double inv3 = inv2 * inv;
      const double inv5 = inv3 * inv2;   // inv^5
      const double inv7 = inv5 * inv2;   // inv^7

      _d0 = 0.5;
      _d1 = (35.0 / 32.0) * inv;     // s term → inv
      _d3 = (-35.0 / 32.0) * inv3;   // s^3 term → inv^3
      _d5 = (21.0 / 32.0) * inv5;    // s^5 term → inv^5
      _d7 = (-5.0  / 32.0) * inv7;   // s^7 term → inv^7
    }

    inline double GetSmoothStepFunctionC3(double dg1) const noexcept {
      if (dg1 <= -_eps) return 0.0;
      if (dg1 >=  _eps) return 1.0;
      const double dg2 = dg1 * dg1;
      // Horner: d0 + dg1*(d1 + dg2*(d3 + dg2*(d5 + dg2*d7)))
      return _d0 + dg1 * (_d1 + dg2 * (_d3 + dg2 * (_d5 + dg2 * _d7)));
    }

    // ---------------- Symmetric sigmoids in [-1,1] from the step variants ----------------
    inline double SigmoidC4(double dg1) const noexcept {
      if (dg1 < -_eps) return -1.0;
      if (dg1 >  _eps) return  1.0;
      return -1. + 2. * GetSmoothStepFunctionC4(dg1);
    }
    inline double SigmoidC3(double dg1) const noexcept {
      if (dg1 < -_eps) return -1.0;
      if (dg1 >  _eps) return  1.0;
      return -1. + 2. * GetSmoothStepFunctionC3(dg1);
    }
    inline double SigmoidC2(double dg1) const noexcept {
      if (dg1 < -_eps) return -1.0;
      if (dg1 >  _eps) return  1.0;
      return -1. + 2. * GetSmoothStepFunctionC2(dg1);
    }
    inline double SigmoidC1(double dg1) const noexcept {
      if (dg1 < -_eps) return -1.0;
      if (dg1 >  _eps) return  1.0;
      return -1. + 2. * GetSmoothStepFunctionC1(dg1);
    }

    inline double eps() const noexcept { return _eps; }

    inline double c1() const noexcept { return _c1; }

  private:
    // Shared half-width of the transition (exact 0/1 outside ±eps)
    double _eps{1.0};

    // C^4 coefficients (S4(s) expanded and scaled to dg1 powers)
    double _a0{0.5}, _a1{0.0}, _a3{0.0}, _a5{0.0}, _a7{0.0}, _a9{0.0};

    // C^2 coefficients
    double _b0{0.5}, _b1{0.0}, _b3{0.0}, _b5{0.0};

    // C^1 coefficients
    double _c0{0.5}, _c1{0.0}, _c3{0.0};

    // C^3 coefficients
    double _d0{0.5}, _d1{0.0}, _d3{0.0}, _d5{0.0}, _d7{0.0};
};
