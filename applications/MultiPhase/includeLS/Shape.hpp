#pragma once

#include "Mollifier.hpp"

#include <cassert>
#include <cmath>
#include <vector>

class Shape {
  public:
    Shape() = default;
    virtual ~Shape() = default;

    virtual double GetIndicator(const std::vector<double>& x,
                                const Mollifier& m) const = 0;
    virtual double GetDistance(const std::vector<double>& x) const = 0;

    virtual void GetMarkers(std::vector<std::vector<double>>& X, const double marker_density) const = 0;
};

class Circle : public Shape {
  public:
    Circle(const std::vector<double>& xc, const double r)
      : _xc(xc), _r(r) {
      assert(_xc.size() == _dim);
      assert(_r > 0.);
    }

    double GetIndicator(const std::vector<double>& x,
                        const Mollifier& m) const override {

      assert(x.size() == _dim);

      double dist = 0.;

      for(unsigned d = 0; d < _dim; d++) {
        dist += (x[d] - _xc[d]) * (x[d] - _xc[d]);
      }

      dist = std::sqrt(dist);

      return m.Sigmoid(_r - dist);
    }

    double GetDistance(const std::vector<double>& x) const override {

      assert(x.size() == _dim);

      double dist = 0.;

      for(unsigned d = 0; d < _dim; d++) {
        dist += (x[d] - _xc[d]) * (x[d] - _xc[d]);
      }

      dist = std::sqrt(dist);

      return _r - dist;
    }

    void GetMarkers(std::vector<std::vector<double>>& X, const double marker_density) const override
    {
      if(marker_density <= 0.) {
        return;
      }

      const double pi = std::acos(-1.);

      const double angular_density = marker_density * _r;

      unsigned nmarkers = 0;

      if(_dim == 2) {
        nmarkers = std::max(
            1u,
            static_cast<unsigned>(
                std::ceil(2. * pi * angular_density)
            )
        );
      }
      else if(_dim == 3) {
        nmarkers = std::max(
            1u,
            static_cast<unsigned>(
                std::ceil(4. * pi *
                          angular_density * angular_density)
            )
        );
      }

      for (int d = 0; d < _dim; ++d)
        X[d].insert(X[d].end(), nmarkers, 0.0);

      if(_dim == 2) {

        const double dtheta =
            2. * pi / static_cast<double>(nmarkers);

        for(unsigned i = 0; i < nmarkers; ++i) {

          const double theta =
              static_cast<double>(i) * dtheta;

          X[0][i] = _xc[0] + _r * std::cos(theta);
          X[1][i] = _xc[1] + _r * std::sin(theta);
        }
      }
      else if(_dim == 3) {

        const double golden_angle =
            pi * (3. - std::sqrt(5.));

        for(unsigned i = 0; i < nmarkers; ++i) {

          const double z =
              1. - 2. *
              (static_cast<double>(i) + 0.5) /
              static_cast<double>(nmarkers);

          const double rho =
              std::sqrt(std::max(0., 1. - z * z));

          const double theta =
              golden_angle * static_cast<double>(i);

          X[0][i] =
              _xc[0] + _r * rho * std::cos(theta);

          X[1][i] =
              _xc[1] + _r * rho * std::sin(theta);

          X[2][i] =
              _xc[2] + _r * z;
        }
      }
    }

    void SetCenter(const std::vector<double>& xc) {
      assert(xc.size() == _dim);
      _xc = xc;
    }

    void SetRadius(const double r) {
      assert(r > 0.);
      _r = r;
    }

    const std::vector<double>& GetCenter() const {
      return _xc;
    }

    double GetRadius() const {
      return _r;
    }

  private:
    std::vector<double> _xc;
    double _r;

    static constexpr unsigned _dim = 2;
};
