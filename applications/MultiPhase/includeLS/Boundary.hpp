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

    virtual double getValue(std::vector<double>& x, const double time, const double period, const double dt) const {
      return -1.;
    };

    virtual void updateMarkers(std::vector<std::vector<std::vector<double>>>& X0, std::vector<std::vector<double>>& X, const double time, const double period, const double dt) const {
      return;
    }
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

  inline double getValue(std::vector<double>& x,
                  const double time,
                  const double period,
                  const double dt) const override {

    double min_distance = 1.e10;
    double value = -1.;

    for (unsigned s = 0; s < _shape.size(); ++s) {

      const unsigned nEmissions =
        static_cast<unsigned>(std::floor(time / _timeOffset[s]));

      for (unsigned t = 0; t <= nEmissions; ++t) {

        // Time at which this copy of the shape entered the domain
        const double injectionTime = static_cast<double>(t) * _timeOffset[s];

        // Backward integration interval: current time -> injection time
        const double integrationTime = time - injectionTime;

        // Always start from the ORIGINAL spatial point
        std::vector<double> xBack = x;

        const unsigned nSteps = std::max(1u, static_cast<unsigned>(std::ceil(std::abs(integrationTime) / std::abs(dt))));

        rk4(xBack, time, -integrationTime, period, nSteps);

        const double distance =
          _shape[s]->GetDistance(xBack);

        if (std::abs(distance) < std::abs(min_distance)) {
          min_distance = distance;
          value = _m.Sigmoid(distance);
        }
      }
    }

    return value;
  }

  void updateMarkers(std::vector<std::vector<std::vector<double>>>& X0, std::vector<std::vector<double>>& X, const double time, const double period, const double dt) const override {

    
    for (unsigned s = 0; s < _shape.size(); ++s) {

      for (unsigned m = 0; m < X0[s][0].size(); m++) {
        std::vector<double> x (X0[s].size());
        for (unsigned d = 0; d < X0[s].size(); d++)
          x[d] = X0[s][d][m];

        const unsigned nEmissions =
          static_cast<unsigned>(std::floor(time / _timeOffset[s]));

        for (unsigned t = (nEmissions >= 1) ? nEmissions - 1 : 0; t <= nEmissions; ++t) {

          // Time at which this copy of the shape entered the domain
          const double injectionTime = static_cast<double>(t) * _timeOffset[s];

          // Backward integration interval: current time -> injection time
          const double integrationTime = time - injectionTime;

          // Always start from the ORIGINAL spatial point
          std::vector<double> xForward = x;

          const unsigned nSteps = std::max(1u, static_cast<unsigned>(std::ceil(std::abs(integrationTime) / std::abs(dt))));

          rk4(xForward, injectionTime, integrationTime, period, nSteps);

          for (unsigned d = 0; d < X.size(); d++)
            X[d].push_back(xForward[d]);
        }
      }
    }
  }



inline void rk4(std::vector<double>& Xp,
                const double t0,
                const double dt,
                const double period,
                const unsigned nSteps) const {

  const std::size_t dim = Xp.size();

  if (dim != 2u && dim != 3u) {
    throw std::runtime_error("RungeKutta::rk4: Xp.size() must be 2 or 3");
  }

  if (!_velocity) {
    throw std::runtime_error("RungeKutta::rk4: velocity function is not set");
  }

  if (nSteps == 0u) {
    throw std::runtime_error("RungeKutta::rk4: nSteps must be greater than zero");
  }

  const double h      = dt / static_cast<double>(nSteps);
  const double halfh  = 0.5 * h;
  const double sixthh = h / 6.0;

  double t = t0;

  std::vector<double> X2(dim);
  std::vector<double> X3(dim);
  std::vector<double> X4(dim);

  for (unsigned step = 0; step < nSteps; ++step) {

    auto vel1 = _velocity(Xp, t, period);

    for (std::size_t d = 0; d < dim; ++d) {
      X2[d] = Xp[d] + halfh * vel1[d];
    }

    auto vel2 = _velocity(X2, t + halfh, period);

    for (std::size_t d = 0; d < dim; ++d) {
      X3[d] = Xp[d] + halfh * vel2[d];
    }

    auto vel3 = _velocity(X3, t + halfh, period);

    for (std::size_t d = 0; d < dim; ++d) {
      X4[d] = Xp[d] + h * vel3[d];
    }

    auto vel4 = _velocity(X4, t + h, period);

    if (vel1.size() != dim ||
        vel2.size() != dim ||
        vel3.size() != dim ||
        vel4.size() != dim) {
      throw std::runtime_error(
        "RungeKutta::rk4: velocity dimension does not match Xp dimension"
      );
    }

    for (std::size_t d = 0; d < dim; ++d) {
      Xp[d] += sixthh *
               (vel1[d]
              + 2.0 * vel2[d]
              + 2.0 * vel3[d]
              + vel4[d]);
    }

    t += h;
  }
}


  private:
    VelocityFunction _velocity;
    const std::vector<Shape*> _shape;
    const std::vector<double>& _timeOffset;
    const Mollifier& _m;
};
