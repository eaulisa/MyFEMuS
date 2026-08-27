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
