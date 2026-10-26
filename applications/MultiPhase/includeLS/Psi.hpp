#pragma once

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <vector>

#include "BurstingBubble.hpp"

enum class PsiType {
  BALL,
  SAYE,
  RAYLEIGH_TAYLOR,
  DAM_BREAK,
  BUBBLES_FREE_SURFACE,
  TWO_BUBBLES,
  BURSTING_BUBBLE
};

struct PsiConfig {
  PsiType type = PsiType::BALL;
  std::vector<double> center;
  double radius = 0.;
  double amplitude = 0.;
  double waveNumber = 0.;
  double damWidth = 0.;
  double damHeight = 0.;
  std::vector<std::vector<double>> bubbleCenters;
  std::vector<double> bubbleRadii;
  double freeSurfaceHeight = 0.;
  double bondNumber = 0.;
  unsigned burstingPointLess = 200u;
  double burstingOuterRadius = 10.;
};

class Psi {

  public:

    explicit Psi(const PsiConfig& config) : _cfg(config) {
      if(_cfg.type == PsiType::BURSTING_BUBBLE) {
        BurstingBubbleConfig bb;
        bb.bondNumber = _cfg.bondNumber;
        bb.pointLess = _cfg.burstingPointLess;
        bb.outerRadius = _cfg.burstingOuterRadius;
        bb.freeSurfaceHeight = _cfg.freeSurfaceHeight;
        _burstingBubble = ComputeBurstingBubbleGeometry(bb);
      }
    }

    double LevelSet(const std::vector<double>& x) const {
      switch(_cfg.type) {

        case PsiType::BALL: {
          if(x.size() != _cfg.center.size()) {
            throw std::runtime_error("Psi::LevelSet(): point dimension does not match ball dimension");
          }
          double r2 = 0.;
          for(std::size_t d = 0; d < x.size(); ++d) {
            const double dx = x[d] - _cfg.center[d];
            r2 += dx * dx;
          }
          return _cfg.radius - std::sqrt(r2);
        }

        case PsiType::SAYE: {
          return x[1] - _cfg.amplitude * std::cos(_cfg.waveNumber * x[0]);
        }

        case PsiType::RAYLEIGH_TAYLOR: {
          const unsigned dim = x.size();
          double psi = 0.;
          if(dim == 2u) psi = x[1] + _cfg.amplitude * std::cos(_cfg.waveNumber * x[0]);
          else if(dim == 3u) psi = x[2] + _cfg.amplitude * (std::cos(_cfg.waveNumber * x[0]) + std::cos(_cfg.waveNumber * x[1]));
          else throw std::runtime_error("Psi::LevelSet(): RAYLEIGH_TAYLOR requires dimension 2 or 3");
          return psi;
        }

        case PsiType::DAM_BREAK: {
          return std::min(_cfg.damWidth - x[0], _cfg.damHeight - x[1]);
        }

        case PsiType::BUBBLES_FREE_SURFACE: {
          double psi = x[1] - _cfg.freeSurfaceHeight;
          if(_cfg.bubbleCenters.size() != _cfg.bubbleRadii.size()) {
            throw std::runtime_error("Psi::LevelSet(): inconsistent bubbleCenters and bubbleRadii sizes");
          }
          for(std::size_t k = 0; k < _cfg.bubbleCenters.size(); ++k) {
            if(_cfg.bubbleCenters[k].size() != x.size()) {
              throw std::runtime_error("Psi::LevelSet(): bubble center dimension mismatch");
            }
            double r2 = 0.;
            for(std::size_t d = 0; d < x.size(); ++d) {
              const double dx = x[d] - _cfg.bubbleCenters[k][d];
              r2 += dx * dx;
            }
            const double psiBubble = _cfg.bubbleRadii[k] - std::sqrt(r2);
            psi = std::max(psi, psiBubble);
          }
          return psi;
        }

        case PsiType::TWO_BUBBLES: {
          if(_cfg.bubbleCenters.size() != _cfg.bubbleRadii.size()) {
            throw std::runtime_error("Psi::LevelSet(): inconsistent bubbleCenters and bubbleRadii sizes");
          }
          double psi = -1.e12;
          for(std::size_t k = 0; k < _cfg.bubbleCenters.size(); ++k) {
            if(_cfg.bubbleCenters[k].size() != x.size()) {
              throw std::runtime_error("Psi::LevelSet(): bubble center dimension mismatch");
            }
            double r2 = 0.;
            for(std::size_t d = 0; d < x.size(); ++d) {
              const double dx = x[d] - _cfg.bubbleCenters[k][d];
              r2 += dx * dx;
            }
            const double psiBubble = _cfg.bubbleRadii[k] - std::sqrt(r2);
            psi = std::max(psi, psiBubble);
          }
          return psi;
        }

        case PsiType::BURSTING_BUBBLE: {
          return BurstingBubbleLevelSet(x, _burstingBubble);
        }
      }

      throw std::runtime_error("Psi::LevelSet(): invalid PsiType");
    }

    double operator()(const std::vector<double>& x) const {
      return LevelSet(x);
    }

    PsiType Type() const {
      return _cfg.type;
    }

  private:

    PsiConfig _cfg;
    BurstingBubbleGeometry _burstingBubble;
};