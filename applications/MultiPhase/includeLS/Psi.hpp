
// #pragma once

// #include <vector>
// #include <stdexcept>
// #include <cmath>
// #include "Mollifier.hpp"

// // assumes Mollifier exists and has: Mollifier(double eps); double SigmoidC1(double d) const;

// // struct PsiBall {
// //   std::vector<double> _c; // center
// //   double              _r;
// //   double            _eps;
// //   Mollifier           _m;

// //   PsiBall(const std::vector<double>& center, const double radius, const Mollifier &m)
// //     : _c(center), _r(radius), _m(m)
// //   {
// //     if (_c.empty())   throw std::runtime_error("PsiBall: center is empty");
// //     if (!(_r > 0.0))  throw std::runtime_error("PsiBall: r must be > 0");
// //   }

// //   double operator()(const std::vector<double>& x) const {
// //     if (x.size() != _c.size()) {
// //       throw std::runtime_error("PsiBall::operator(): x.size() != center.size()");
// //     }

// //     double s2 = 0.0;
// //     for (std::size_t k = 0; k < _c.size(); ++k) {
// //       const double d = x[k] - _c[k];
// //       s2 += d * d;
// //     }

// //     const double d = _r - std::sqrt(s2);
// //     return _m.Sigmoid(d);
// //   }
// // };

// #pragma once

// #include <cmath>
// #include <stdexcept>
// #include <vector>

// #include "Mollifier.hpp"

// struct PsiBall {
//   std::vector<double> _c;
//   double              _r;
//   Mollifier           _m;

//   PsiBall(const std::vector<double>& center,
//           const double radius,
//           const Mollifier& m)
//     : _c(center), _r(radius), _m(m) {
//     if (_c.empty()) {
//       throw std::runtime_error("PsiBall: center is empty");
//     }

//     if (!(_r > 0.0)) {
//       throw std::runtime_error("PsiBall: r must be > 0");
//     }
//   }

//   double LevelSet(const std::vector<double>& x) const {
//     if (x.size() != _c.size()) {
//       throw std::runtime_error(
//         "PsiBall::LevelSet(): x.size() != center.size()");
//     }

//     double r2 = 0.0;

//     for (std::size_t k = 0; k < _c.size(); ++k) {
//       const double dx = x[k] - _c[k];
//       r2 += dx * dx;
//     }

//     return _r - std::sqrt(r2);
//   }

//   double SignedDistance(const std::vector<double>& x) const {
//     return LevelSet(x);
//   }

//   std::vector<double> Normal(const std::vector<double>& x) const {
//     if (x.size() != _c.size()) {
//       throw std::runtime_error(
//         "PsiBall::Normal(): x.size() != center.size()");
//     }

//     std::vector<double> n(_c.size(), 0.0);
//     double norm2 = 0.0;

//     for (std::size_t k = 0; k < _c.size(); ++k) {
//       n[k] = x[k] - _c[k];
//       norm2 += n[k] * n[k];
//     }

//     const double norm = std::sqrt(norm2);

//     if (norm < 1.e-14) {
//       throw std::runtime_error(
//         "PsiBall::Normal(): normal undefined at center");
//     }

//     for (std::size_t k = 0; k < _c.size(); ++k) {
//       n[k] /= norm;
//     }

//     return n;
//   }

//   double Curvature(const std::vector<double>& x) const {
//     if (x.size() != _c.size()) {
//       throw std::runtime_error(
//         "PsiBall::Curvature(): x.size() != center.size()");
//     }

//     return 1.0 / _r;
//   }

//   double Curvature() const {
//     return 1.0 / _r;
//   }

//   double operator()(const std::vector<double>& x) const {
//     return _m.Sigmoid(LevelSet(x));
//   }
// };

// struct PsiEllipse {
//   std::vector<double> _c;
//   double              _a;
//   double              _b;
//   Mollifier           _m;

//   PsiEllipse(const std::vector<double>& center,
//              const double a,
//              const double b,
//              const Mollifier& m)
//     : _c(center), _a(a), _b(b), _m(m) {
//     if (_c.size() < 2) {
//       throw std::runtime_error(
//         "PsiEllipse: center must have at least 2 components");
//     }

//     if (!(_a > 0.0)) {
//       throw std::runtime_error("PsiEllipse: a must be > 0");
//     }

//     if (!(_b > 0.0)) {
//       throw std::runtime_error("PsiEllipse: b must be > 0");
//     }
//   }

//   double LevelSet(const std::vector<double>& x) const {
//     if (x.size() != _c.size()) {
//       throw std::runtime_error(
//         "PsiEllipse::LevelSet(): x.size() != center.size()");
//     }

//     const double dx = x[0] - _c[0];
//     const double dy = x[1] - _c[1];

//     return 1.0
//            - (dx * dx) / (_a * _a)
//            - (dy * dy) / (_b * _b);
//   }

//   std::vector<double> Normal(const double theta) const {
//     const double ct = std::cos(theta);
//     const double st = std::sin(theta);

//     const double den =
//       std::sqrt(_b * _b * ct * ct + _a * _a * st * st);

//     return {_b * ct / den, _a * st / den};
//   }

//   std::vector<double> Normal(const std::vector<double>& x) const {
//     if (x.size() != _c.size()) {
//       throw std::runtime_error(
//         "PsiEllipse::Normal(): x.size() != center.size()");
//     }

//     const double dx = x[0] - _c[0];
//     const double dy = x[1] - _c[1];

//     const double theta = std::atan2(dy / _b, dx / _a);

//     return Normal(theta);
//   }

//   double Curvature(const double theta) const {
//     const double ct = std::cos(theta);
//     const double st = std::sin(theta);

//     const double den =
//       _a * _a * st * st +
//       _b * _b * ct * ct;

//     return _a * _b / std::pow(den, 1.5);
//   }

//   double Curvature(const std::vector<double>& x) const {
//     if (x.size() != _c.size()) {
//       throw std::runtime_error(
//         "PsiEllipse::Curvature(): x.size() != center.size()");
//     }

//     const double dx = x[0] - _c[0];
//     const double dy = x[1] - _c[1];

//     const double theta = std::atan2(dy / _b, dx / _a);

//     return Curvature(theta);
//   }

//   double operator()(const std::vector<double>& x) const {
//     return LevelSet(x);
//   }
// };

// struct PsiStar {
//   std::vector<double> _c;
//   unsigned            _n;
//   Mollifier           _m;

//   PsiStar(const std::vector<double>& center,
//           const unsigned n,
//           const Mollifier& m)
//     : _c(center), _n(n), _m(m) {
//     if (_c.size() < 2) {
//       throw std::runtime_error(
//         "PsiStar: center must have at least 2 components");
//     }

//     if (_n == 0) {
//       throw std::runtime_error("PsiStar: n must be > 0");
//     }
//   }

//   double Radius(const double theta) const {
//     return 0.25
//            + 0.1 * std::cos(static_cast<double>(_n) * theta + 2.0);
//   }

//   double LevelSet(const std::vector<double>& x) const {
//     if (x.size() != _c.size()) {
//       throw std::runtime_error(
//         "PsiStar::LevelSet(): x.size() != center.size()");
//     }

//     const double dx = x[0] - _c[0];
//     const double dy = x[1] - _c[1];

//     const double rho   = std::sqrt(dx * dx + dy * dy);
//     const double theta = std::atan2(dy, dx);

//     return Radius(theta) - rho;
//   }

//   std::vector<double> Normal(const double theta) const {
//     const double n = static_cast<double>(_n);
//     const double alpha = n * theta + 2.0;

//     const double r = 0.25 + 0.1 * std::cos(alpha);
//     const double rp = -0.1 * n * std::sin(alpha);

//     const double den = std::sqrt(r * r + rp * rp);

//     const double nx =
//       (r * std::cos(theta) + rp * std::sin(theta)) / den;

//     const double ny =
//       (r * std::sin(theta) - rp * std::cos(theta)) / den;

//     return {nx, ny};
//   }

//   std::vector<double> Normal(const std::vector<double>& x) const {
//     if (x.size() != _c.size()) {
//       throw std::runtime_error(
//         "PsiStar::Normal(): x.size() != center.size()");
//     }

//     const double dx = x[0] - _c[0];
//     const double dy = x[1] - _c[1];

//     const double theta = std::atan2(dy, dx);

//     return Normal(theta);
//   }

//   double Curvature(const double theta) const {
//     const double n = static_cast<double>(_n);
//     const double alpha = n * theta + 2.0;

//     const double r = 0.25 + 0.1 * std::cos(alpha);
//     const double rp = -0.1 * n * std::sin(alpha);
//     const double rpp = -0.1 * n * n * std::cos(alpha);

//     return
//       (r * r + 2.0 * rp * rp - r * rpp) /
//       std::pow(r * r + rp * rp, 1.5);
//   }

//   double Curvature(const std::vector<double>& x) const {
//     if (x.size() != _c.size()) {
//       throw std::runtime_error(
//         "PsiStar::Curvature(): x.size() != center.size()");
//     }

//     const double dx = x[0] - _c[0];
//     const double dy = x[1] - _c[1];

//     const double theta = std::atan2(dy, dx);

//     return Curvature(theta);
//   }

//   double operator()(const std::vector<double>& x) const {
//     return LevelSet(x);
//   }
// };

// struct PsiSaye {
//   Mollifier _m;

//   static constexpr double _a = 0.025;
//   static constexpr double _pi = 3.141592653589793238462643383279502884;
//   static constexpr double _k = 2.0 * _pi;

//   PsiSaye(const Mollifier& m)
//     : _m(m)
//   {}

//   double InterfaceHeight(const double x) const {
//     return _a * std::cos(_k * x);
//   }

//   double LevelSet(const std::vector<double>& x) const {
//     if (x.size() != 2) {
//       throw std::runtime_error(
//         "PsiSaye::LevelSet(): Saye test is two-dimensional");
//     }

//     return x[1] - InterfaceHeight(x[0]);
//   }

//   std::vector<double> Normal(const std::vector<double>& x) const {
//     if (x.size() != 2) {
//       throw std::runtime_error(
//         "PsiSaye::Normal(): Saye test is two-dimensional");
//     }

//     const double sx = _a * _k * std::sin(_k * x[0]);
//     const double norm = std::sqrt(1.0 + sx * sx);

//     return {sx / norm, 1.0 / norm};
//   }

//   double Curvature(const std::vector<double>& x) const {
//     if (x.size() != 2) {
//       throw std::runtime_error(
//         "PsiSaye::Curvature(): Saye test is two-dimensional");
//     }

//     const double sx = _a * _k * std::sin(_k * x[0]);

//     return _a * _k * _k * std::cos(_k * x[0])
//            / std::pow(1.0 + sx * sx, 1.5);
//   }

//   double operator()(const std::vector<double>& x) const {
//     return LevelSet(x);
//   }

//   double Mollified(const std::vector<double>& x) const {
//     return _m.Sigmoid(LevelSet(x));
//   }
// };

// struct PsiRayleighTaylor {
//   Mollifier _m;

//   static constexpr double _a  = 0.005;
//   static constexpr double _pi = 3.141592653589793238462643383279502884;
//   static constexpr double _k  = 2.0 * _pi;

//   PsiRayleighTaylor(const Mollifier& m)
//     : _m(m) {

//   }

//   double InterfaceHeight(const double x) const {
//     return  _a * std::cos(_k * x);
//   }

//   double LevelSet(const std::vector<double>& x) const {
//     if (x.size() != 2) {
//       throw std::runtime_error(
//         "PsiRayleighTaylor::LevelSet(): "
//         "Rayleigh-Taylor test is two-dimensional");
//     }

//     return x[1] - InterfaceHeight(x[0]);
//   }

//   std::vector<double> Normal(
//     const std::vector<double>& x) const {

//     if (x.size() != 2) {
//       throw std::runtime_error(
//         "PsiRayleighTaylor::Normal(): "
//         "Rayleigh-Taylor test is two-dimensional");
//     }

//     // Interface:
//     //
//     // y = -h/2 + a cos(k x)
//     //
//     // f'(x) = -a k sin(k x)
//     //
//     // grad(phi) = (-f'(x), 1)
//     //           = (a k sin(k x), 1)

//     const double sx =
//       _a * _k * std::sin(_k * x[0]);

//     const double norm =
//       std::sqrt(1.0 + sx * sx);

//     return {
//       sx / norm,
//       1.0 / norm
//     };
//   }

//   double Curvature(
//     const std::vector<double>& x) const {

//     if (x.size() != 2) {
//       throw std::runtime_error(
//         "PsiRayleighTaylor::Curvature(): "
//         "Rayleigh-Taylor test is two-dimensional");
//     }

//     const double sx =
//       _a * _k * std::sin(_k * x[0]);

//     return
//       _a * _k * _k * std::cos(_k * x[0])
//       /
//       std::pow(1.0 + sx * sx, 1.5);
//   }

//   double operator()(
//     const std::vector<double>& x) const {

//     return LevelSet(x);
//   }

//   double Mollified(
//     const std::vector<double>& x) const {

//     return _m.Sigmoid(LevelSet(x));
//   }
// };

#pragma once

#include <cmath>
#include <stdexcept>
#include <vector>

enum class PsiType {
  BALL,
  SAYE,
  RAYLEIGH_TAYLOR,
  DAM_BREAK,
  BUBBLES_FREE_SURFACE,
  TWO_BUBBLES
};

struct PsiConfig {

  PsiType type = PsiType::BALL;

  // Ball
  std::vector<double> center;
  double radius = 0.;

  // Sinusoidal interfaces
  double amplitude = 0.;
  double waveNumber = 0.;

  // Dam break
  double damWidth  = 0.;
  double damHeight = 0.;

  // Rising bubble(s) below a free surface
  std::vector<std::vector<double>> bubbleCenters;
  std::vector<double> bubbleRadii;
  double freeSurfaceHeight = 0.;
};

class Psi {

  public:

    // ============================================================
    // Constructor
    // ============================================================

    explicit Psi(const PsiConfig& config)
      : _cfg(config) {

      Validate();
    }

    // ============================================================
    // Level-set value
    // ============================================================

    double LevelSet(const std::vector<double>& x) const {

      switch(_cfg.type) {

        // ========================================================
        // Ball / circle
        //
        // psi(x) = R - |x - xc|
        //
        // psi > 0 inside
        // psi = 0 on interface
        // psi < 0 outside
        // ========================================================

        case PsiType::BALL: {

          if(x.size() != _cfg.center.size()) {
            throw std::runtime_error(
              "Psi::LevelSet(): "
              "point dimension does not match ball dimension");
          }

          double r2 = 0.;

          for(std::size_t d = 0;
              d < x.size();
              ++d) {

            const double dx =
              x[d] - _cfg.center[d];

            r2 += dx * dx;
          }

          return
            _cfg.radius
            - std::sqrt(r2);
        }

        // ========================================================
        // Saye interface
        //
        // y = a cos(k x)
        //
        // psi(x,y) = y - a cos(k x)
        //
        // In 3D:
        //
        // psi(x,y,z) = y - a cos(k x)
        //
        // i.e. extrusion along z.
        // ========================================================

        case PsiType::SAYE: {

          CheckSinusoidalPoint(x);

          return
            x[1]
            - _cfg.amplitude
            * std::cos(
              _cfg.waveNumber * x[0]);
        }

        // ========================================================
        // Rayleigh-Taylor interface
        //
        // y = a cos(k x)
        //
        // psi(x,y) = y - a cos(k x)
        // ========================================================

        case PsiType::RAYLEIGH_TAYLOR: {

          CheckSinusoidalPoint(x);

          const unsigned dim = x.size();

          double psi = 0;
          if (dim == 2) psi = x[1] + _cfg.amplitude * std::cos(_cfg.waveNumber * x[0]);
          else psi = x[2] + _cfg.amplitude * (std::cos(_cfg.waveNumber * x[0]) + std::cos(_cfg.waveNumber * x[1]));

          return psi;
        }

        // ========================================================
        // Dam break
        //
        // Initial liquid column:
        //
        // 0 <= x <= damWidth
        // 0 <= y <= damHeight
        //
        // The left and bottom sides coincide with domain walls,
        // therefore only the two fluid-fluid interfaces are:
        //
        // x = damWidth
        // y = damHeight
        //
        // psi > 0 inside the liquid.
        //
        // This is not an exact Euclidean signed distance near the
        // upper-right corner, but it gives the exact zero level set.
        // ========================================================

        case PsiType::DAM_BREAK: {

          CheckDamBreakPoint(x);

          return
            std::min(
              _cfg.damWidth  - x[0],
              _cfg.damHeight - x[1]);
        }

        case PsiType::BUBBLES_FREE_SURFACE: {

          // Positive phase:
          //
          //   1) gas above the free surface
          //   2) gas inside any bubble
          //
          // The union of level-set regions is obtained with max().

          double psi =
            x[1] - _cfg.freeSurfaceHeight;

          for(std::size_t k = 0;
              k < _cfg.bubbleCenters.size();
              ++k) {

            double r2 = 0.;

            for(std::size_t d = 0;
                d < x.size();
                ++d) {

              const double dx =
                x[d] - _cfg.bubbleCenters[k][d];

              r2 += dx * dx;
            }

            const double psiBubble =
              _cfg.bubbleRadii[k]
              - std::sqrt(r2);

            psi = std::max(
                    psi,
                    psiBubble);
          }

          return psi;
        }



        case PsiType::TWO_BUBBLES: {

          double psi = - 1.e12;;

          for(std::size_t k = 0;
              k < _cfg.bubbleCenters.size();
          ++k) {

            double r2 = 0.;

            for(std::size_t d = 0;
                d < x.size();
            ++d) {

              const double dx =
              x[d] - _cfg.bubbleCenters[k][d];

              r2 += dx * dx;
            }

            const double psiBubble =
            _cfg.bubbleRadii[k]
            - std::sqrt(r2);

            psi = std::max(
              psi,
              psiBubble);
          }

          return psi;
        }
      }


      throw std::runtime_error(
        "Psi::LevelSet(): invalid PsiType");
    }

    // ============================================================
    // Function-call operator
    // ============================================================

    double operator()(
      const std::vector<double>& x) const {

      return LevelSet(x);
    }

    // ============================================================
    // Type
    // ============================================================

    PsiType Type() const {

      return _cfg.type;
    }

  private:

    PsiConfig _cfg;

    // ============================================================
    // Validation
    // ============================================================

    void Validate() const {

      switch(_cfg.type) {

        case PsiType::BALL: {

          if(_cfg.center.empty()) {
            throw std::runtime_error(
              "Psi: BALL center is empty");
          }

          if(_cfg.center.size() != 2 &&
              _cfg.center.size() != 3) {
            throw std::runtime_error(
              "Psi: BALL must be either 2D or 3D");
          }

          if(!(_cfg.radius > 0.)) {
            throw std::runtime_error(
              "Psi: BALL radius must be > 0");
          }

          break;
        }

        case PsiType::SAYE: {

          if(!(_cfg.waveNumber > 0.)) {
            throw std::runtime_error(
              "Psi: SAYE waveNumber must be > 0");
          }

          break;
        }

        case PsiType::RAYLEIGH_TAYLOR: {

          if(!(_cfg.waveNumber > 0.)) {
            throw std::runtime_error(
              "Psi: RAYLEIGH_TAYLOR "
              "waveNumber must be > 0");
          }

          break;
        }

        case PsiType::DAM_BREAK: {

          if(!(_cfg.damWidth > 0.)) {
            throw std::runtime_error(
              "Psi: DAM_BREAK damWidth must be > 0");
          }

          if(!(_cfg.damHeight > 0.)) {
            throw std::runtime_error(
              "Psi: DAM_BREAK damHeight must be > 0");
          }

          break;
        }
      }
    }

    // ============================================================
    // Point validation for sinusoidal interfaces
    // ============================================================

    void CheckSinusoidalPoint(
      const std::vector<double>& x) const {

      if(x.size() != 2 &&
          x.size() != 3) {

        throw std::runtime_error(
          "Psi::LevelSet(): "
          "sinusoidal interface must be evaluated "
          "in 2D or 3D");
      }
    }

    void CheckDamBreakPoint(
      const std::vector<double>& x) const {

      if(x.size() != 2) {

        throw std::runtime_error(
          "Psi::LevelSet(): "
          "DAM_BREAK is currently implemented only in 2D");
      }
    }
};
