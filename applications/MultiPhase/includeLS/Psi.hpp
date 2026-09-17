 
#pragma once

#include <vector>
#include <stdexcept>
#include <cmath>
#include "Mollifier.hpp"

// assumes Mollifier exists and has: Mollifier(double eps); double SigmoidC1(double d) const;

// struct PsiBall {
//   std::vector<double> _c; // center
//   double              _r;
//   double            _eps;
//   Mollifier           _m;


//   PsiBall(const std::vector<double>& center, const double radius, const Mollifier &m)
//     : _c(center), _r(radius), _m(m)
//   {
//     if (_c.empty())   throw std::runtime_error("PsiBall: center is empty");
//     if (!(_r > 0.0))  throw std::runtime_error("PsiBall: r must be > 0");
//   }

//   double operator()(const std::vector<double>& x) const {
//     if (x.size() != _c.size()) {
//       throw std::runtime_error("PsiBall::operator(): x.size() != center.size()");
//     }

//     double s2 = 0.0;
//     for (std::size_t k = 0; k < _c.size(); ++k) {
//       const double d = x[k] - _c[k];
//       s2 += d * d;
//     }

//     const double d = _r - std::sqrt(s2);
//     return _m.Sigmoid(d);
//   }
// };

#include <vector>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <algorithm>

struct PsiBall {
  std::vector<double> _c;
  double              _r;
  Mollifier           _m;

  PsiBall(const std::vector<double>& center, const double radius, const Mollifier& m)
    : _c(center), _r(radius), _m(m)
  {
    if (_c.empty()) throw std::runtime_error("PsiBall: center is empty");
    if (!(_r > 0.0)) throw std::runtime_error("PsiBall: r must be > 0");
  }

  double SignedDistance(const std::vector<double>& x) const {
    if (x.size() != _c.size()) {
      throw std::runtime_error("PsiBall::SignedDistance(): x.size() != center.size()");
    }

    double s2 = 0.0;

    for(std::size_t k = 0; k < _c.size(); ++k) {
      const double d = x[k] - _c[k];
      s2 += d * d;
    }

    return _r - std::sqrt(s2);
  }

  std::vector<double> Normal(const std::vector<double>& x) const {
    if(x.size() != _c.size()) throw std::runtime_error("PsiBall::Normal(): x.size() != center.size()");

    std::vector<double> n(_c.size(), 0.);

    double norm2 = 0.;

    for(std::size_t k = 0; k < _c.size(); ++k) {
      n[k] = x[k] - _c[k];
      norm2 += n[k] * n[k];
    }

    const double norm = std::sqrt(norm2);

    if(norm < 1.e-14) throw std::runtime_error("PsiBall::Normal(): normal undefined at center");

    for(std::size_t k = 0; k < _c.size(); ++k) n[k] /= norm;

    return n;
  }

  double Curvature(const std::vector<double>& x) const {
    return 1.0 / _r;
  }

  double Curvature() const {
    return 1.0 / _r;
}

  double operator()(const std::vector<double>& x) const {
    return _m.Sigmoid(SignedDistance(x));
  }
};


struct PsiEllipse {
  std::vector<double> _c;
  double              _a;
  double              _b;
  Mollifier           _m;

  PsiEllipse(const std::vector<double>& center, const double a, const double b, const Mollifier& m)
    : _c(center), _a(a), _b(b), _m(m)
  {
    if (_c.size() < 2) throw std::runtime_error("PsiEllipse: center must have at least 2 components");
    if (!(_a > 0.0)) throw std::runtime_error("PsiEllipse: a must be > 0");
    if (!(_b > 0.0)) throw std::runtime_error("PsiEllipse: b must be > 0");
  }

  double SignedDistance(const std::vector<double>& x) const {
    if (x.size() != _c.size()) {
      throw std::runtime_error("PsiEllipse::SignedDistance(): x.size() != center.size()");
    }

    const double px = std::abs(x[0] - _c[0]);
    const double py = std::abs(x[1] - _c[1]);

    double theta = std::atan2(_a * py, _b * px);

    for(unsigned iter = 0; iter < 30; ++iter) {
      const double ct = std::cos(theta);
      const double st = std::sin(theta);

      const double ex = _a * ct;
      const double ey = _b * st;

      const double dx = ex - px;
      const double dy = ey - py;

      const double dex = -_a * st;
      const double dey =  _b * ct;

      const double ddex = -_a * ct;
      const double ddey = -_b * st;

      const double f = dx * dex + dy * dey;

      const double df =
          dex * dex +
          dx * ddex +
          dey * dey +
          dy * ddey;

      if(std::abs(df) < 1.e-14) break;

      const double dtheta = f / df;

      theta -= dtheta;

      if(std::abs(dtheta) < 1.e-14) break;
    }

    const double ex = _a * std::cos(theta);
    const double ey = _b * std::sin(theta);

    const double dx = ex - px;
    const double dy = ey - py;

    const double distance = std::sqrt(dx * dx + dy * dy);

    const double q =
        px * px / (_a * _a) +
        py * py / (_b * _b);

    return (q <= 1.0) ? distance : -distance;
  }

  std::vector<double> Normal(const double theta) const {
    const double ct = std::cos(theta);
    const double st = std::sin(theta);
    const double den = std::sqrt(_b * _b * ct * ct + _a * _a * st * st);
    return {_b * ct / den, _a * st / den};
  }

  std::vector<double> Normal(const std::vector<double>& x) const {
    const double dx = x[0] - _c[0];
    const double dy = x[1] - _c[1];
    const double theta = std::atan2(dy / _b, dx / _a);
    return Normal(theta);
  }

  double Curvature(const double theta) const {
    const double ct = std::cos(theta);
    const double st = std::sin(theta);
    const double den = _a * _a * st * st + _b * _b * ct * ct;
    return _a * _b / std::pow(den, 1.5);
  }

  double Curvature(const std::vector<double>& x) const {
    const double dx = x[0] - _c[0];
    const double dy = x[1] - _c[1];
    const double theta = std::atan2(dy / _b, dx / _a);
    return Curvature(theta);
  }

  double operator()(const std::vector<double>& x) const {
    return _m.Sigmoid(SignedDistance(x));
  }
};


struct PsiStar {
  std::vector<double> _c;
  unsigned            _n;
  Mollifier           _m;

  PsiStar(const std::vector<double>& center, const unsigned n, const Mollifier& m)
    : _c(center), _n(n), _m(m)
  {
    if (_c.size() < 2) throw std::runtime_error("PsiStar: center must have at least 2 components");
    if (_n == 0) throw std::runtime_error("PsiStar: n must be > 0");
  }

  double Radius(const double theta) const {
    return 0.25 + 0.1 * std::cos(static_cast<double>(_n) * theta + 2.0);
  }

  double SignedDistance(const std::vector<double>& x) const {
    if (x.size() != _c.size()) {
      throw std::runtime_error("PsiStar::SignedDistance(): x.size() != center.size()");
    }

    const double px = x[0] - _c[0];
    const double py = x[1] - _c[1];

    const double rho = std::sqrt(px * px + py * py);
    const double thetaPoint = std::atan2(py, px);

    const unsigned nSamples = std::max(100u, 20u * _n);

    double theta = 0.0;
    double minDist2 = std::numeric_limits<double>::max();

    for(unsigned i = 0; i < nSamples; ++i) {
      const double t =
          -M_PI +
          2.0 * M_PI *
          static_cast<double>(i) /
          static_cast<double>(nSamples);

      const double r = Radius(t);

      const double qx = r * std::cos(t);
      const double qy = r * std::sin(t);

      const double dx = qx - px;
      const double dy = qy - py;

      const double dist2 = dx * dx + dy * dy;

      if(dist2 < minDist2) {
        minDist2 = dist2;
        theta = t;
      }
    }

    const double n = static_cast<double>(_n);

    for(unsigned iter = 0; iter < 30; ++iter) {
      const double alpha = n * theta + 2.0;

      const double r =
          0.25 + 0.1 * std::cos(alpha);

      const double rp =
          -0.1 * n * std::sin(alpha);

      const double rpp =
          -0.1 * n * n * std::cos(alpha);

      const double ct = std::cos(theta);
      const double st = std::sin(theta);

      const double qx = r * ct;
      const double qy = r * st;

      const double qpx = rp * ct - r * st;
      const double qpy = rp * st + r * ct;

      const double qppx =
          (rpp - r) * ct -
          2.0 * rp * st;

      const double qppy =
          (rpp - r) * st +
          2.0 * rp * ct;

      const double dx = qx - px;
      const double dy = qy - py;

      const double f =
          dx * qpx +
          dy * qpy;

      const double df =
          qpx * qpx +
          qpy * qpy +
          dx * qppx +
          dy * qppy;

      if(std::abs(df) < 1.e-14) break;

      const double dtheta = f / df;

      theta -= dtheta;

      if(std::abs(dtheta) < 1.e-14) break;
    }

    const double r = Radius(theta);

    const double qx = r * std::cos(theta);
    const double qy = r * std::sin(theta);

    const double dx = qx - px;
    const double dy = qy - py;

    const double distance = std::sqrt(dx * dx + dy * dy);

    const double rPoint = Radius(thetaPoint);

    return (rho <= rPoint) ? distance : -distance;
  }

  std::vector<double> Normal(const double theta) const {
    const double n = static_cast<double>(_n);
    const double alpha = n * theta + 2.0;
    const double r = 0.25 + 0.1 * std::cos(alpha);
    const double rp = -0.1 * n * std::sin(alpha);
    const double den = std::sqrt(r * r + rp * rp);
    const double nx = (r * std::cos(theta) + rp * std::sin(theta)) / den;
    const double ny = (r * std::sin(theta) - rp * std::cos(theta)) / den;
    return {nx, ny};
  }

  std::vector<double> Normal(const std::vector<double>& x) const {
    const double dx = x[0] - _c[0];
    const double dy = x[1] - _c[1];
    const double theta = std::atan2(dy, dx);
    return Normal(theta);
  }

  double Curvature(const double theta) const {
    const double n = static_cast<double>(_n);
    const double alpha = n * theta + 2.0;
    const double r = 0.25 + 0.1 * std::cos(alpha);
    const double rp = -0.1 * n * std::sin(alpha);
    const double rpp = -0.1 * n * n * std::cos(alpha);
    return (r * r + 2.0 * rp * rp - r * rpp) / std::pow(r * r + rp * rp, 1.5);
  }

  double Curvature(const std::vector<double>& x) const {
    const double dx = x[0] - _c[0];
    const double dy = x[1] - _c[1];
    const double theta = std::atan2(dy, dx);
    return Curvature(theta);
  }

  double operator()(const std::vector<double>& x) const {
    return _m.Sigmoid(SignedDistance(x));
  }
};