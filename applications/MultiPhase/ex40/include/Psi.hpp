 
#pragma once

#include <vector>
#include <stdexcept>
#include <cmath>
#include "Mollifier.hpp"

// assumes Mollifier exists and has: Mollifier(double eps); double SigmoidC1(double d) const;

struct PsiBall {
  std::vector<double> _c; // center
  double              _r;
  double            _eps;
  Mollifier           _m;


  PsiBall(const std::vector<double>& center, double radius, double eps)
    : _c(center), _r(radius),
    _eps(eps)
  {
    _m = Mollifier(_eps);
    if (_c.empty())   throw std::runtime_error("PsiBall: center is empty");
    if (!(_r > 0.0))  throw std::runtime_error("PsiBall: r must be > 0");
    if (!(_eps > 0.0)) throw std::runtime_error("PsiBall: eps must be > 0");
  }

  double operator()(const std::vector<double>& x) const {
    if (x.size() != _c.size()) {
      throw std::runtime_error("PsiBall::operator(): x.size() != center.size()");
    }

    double s2 = 0.0;
    for (std::size_t k = 0; k < _c.size(); ++k) {
      const double d = x[k] - _c[k];
      s2 += d * d;
    }

    const double d = _r - std::sqrt(s2);
    return _m.SigmoidC1(d);
  }
};
