#pragma once

#include <cassert>
#include <cmath>
#include <functional>
#include <vector>

#include "Shape.hpp"
#include "Mollifier.hpp"

enum class InflowVelKind : unsigned {
  Uniform = 0,
  Poiseuille = 1
};

const InflowVelKind inflowVelocityType = InflowVelKind::Poiseuille;

inline std::vector<double>
InflowVelocity(std::vector<double>& xp,
               const double time,
               const double period) noexcept {

  double u = 0.0;
  double v = 0.0;

  switch(inflowVelocityType) {

    case InflowVelKind::Uniform: {
      u = 0.0;
      v = -0.3;
      break;
    }

    case InflowVelKind::Poiseuille: {
      u = 0.0;
      v = -0.3 * 4.0 * (0.5 - xp[0]) * (0.5 + xp[0]);
      break;
    }
  }

  return {u, v};
}

class Boundary {
  public:
    Boundary() = default;
    virtual ~Boundary() = default;

    virtual double getValue(std::vector<double>& x, const double time, const double period) const {
      return -1.;
    };
};

class Inflow : public Boundary {
  public:
    using VelocityFunction =
      std::function<std::vector<double>(std::vector<double>&,
                                        double,
                                        double)>;

    Inflow(VelocityFunction velocity, const std::vector<Shape*> shape, std::vector<double>& timeOffset, const Mollifier& m)
      : _velocity(std::move(velocity)),
        _shape(shape),
        _timeOffset(timeOffset),
        _m(m) {
    }

    double getValue(std::vector<double>& x, const double time, const double period) const override {

      auto vel = _velocity(x, time, period);

      double min_distance = 1.e10;
      double value = -1.;

      for (unsigned s = 0; s < _shape.size(); s++) {
        for (unsigned t  = 0 ; t <= floor(time / _timeOffset[s]); t++) {
          std::vector<double> x0 (x.size(), 0.);
          for (int d = 0; d < x.size(); d++)
            x0[d] = x[d] - vel[d] * (time - t * _timeOffset[s]);

          double distance = _shape[s]->GetDistance(x0);

          if (fabs(distance) < fabs(min_distance)) {
            min_distance = distance;
            value = _m.Sigmoid(distance);
          }
        }
      }
      // upgrade backward integratio
      return value;
    }

  private:
    VelocityFunction _velocity;
    const std::vector<Shape*> _shape;
    const std::vector<double>& _timeOffset;
    const Mollifier& _m;
};
