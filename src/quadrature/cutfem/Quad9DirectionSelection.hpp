#ifndef FEMUS_QUAD9_DIRECTION_SELECTION_HPP
#define FEMUS_QUAD9_DIRECTION_SELECTION_HPP

#include "AdaptiveLineSliceMoments.hpp"
#include "FELevelSetMoments.hpp"

// Unit-square QUAD9 moments: order 0 selects automatically, 1 integrates y first, and 2 integrates x first.
namespace FELevelSetMoments {

inline void Quad9SelectedQuadraticRange(long double a, long double b, long double c, long double& minimum, long double& maximum) {
  const long double endpoint = (a+b)+c;
  minimum = std::min(minimum,std::min(c,endpoint));
  maximum = std::max(maximum,std::max(c,endpoint));
  if(a != 0.0L && ((b < 0.0L && a > 0.0L) || (b > 0.0L && a < 0.0L)) && std::abs(b) < 2.0L*std::abs(a)) {
    const long double point = -b/(2.0L*a);
    const long double value = std::fma(std::fma(a,point,b),point,c);
    minimum = std::min(minimum,value);
    maximum = std::max(maximum,value);
  }
}

inline void Quad9SelectedDerivativeBounds(const std::vector<std::vector<long double>>& coefficient, unsigned axis, long double& minimum, long double& maximum) {
  if(coefficient.size() != 3 || axis > 1) throw std::invalid_argument("Invalid QUAD9 derivative input.");
  long double norm = 0.0L;
  for(unsigned i = 0; i < 3; ++i) {
    if(coefficient[i].size() != 3) throw std::invalid_argument("QUAD9 coefficients must be 3 by 3.");
    for(unsigned j = 0; j < 3; ++j) {
      if(!std::isfinite(coefficient[i][j])) throw std::invalid_argument("A QUAD9 coefficient is not finite.");
      norm += (axis == 0 ? i : j)*std::abs(coefficient[i][j]);
    }
  }
  minimum = std::numeric_limits<long double>::infinity();
  maximum = -minimum;
  for(unsigned side = 0; side < 2; ++side) {
    long double q[3];
    for(unsigned j = 0; j < 3; ++j) q[j] = axis == 0 ? coefficient[1][j]+2.0L*side*coefficient[2][j] : coefficient[j][1]+2.0L*side*coefficient[j][2];
    Quad9SelectedQuadraticRange(q[2],q[1],q[0],minimum,maximum);
  }
  const long double padding = 512.0L*std::numeric_limits<long double>::epsilon()*norm;
  minimum -= padding;
  maximum += padding;
}

inline unsigned Quad9SelectedOrder(const std::vector<std::vector<long double>>& coefficient, long double lower[2], long double upper[2]) {
  for(unsigned axis = 0; axis < 2; ++axis) {
    long double minimum;
    long double maximum;
    Quad9SelectedDerivativeBounds(coefficient,axis,minimum,maximum);
    upper[axis] = std::max(std::abs(minimum),std::abs(maximum));
    lower[axis] = minimum > 0.0L ? minimum : (maximum < 0.0L ? -maximum : 0.0L);
    const long double conditioning = 32.0L*std::sqrt(std::numeric_limits<long double>::epsilon())*upper[axis];
    if(lower[axis] <= conditioning) lower[axis] = 0.0L;
  }
  if(lower[0] != lower[1]) return lower[0] > lower[1] ? 2U : 1U;
  return upper[0] > upper[1] ? 2U : 1U;
}

inline void Quad9SelectedOneDirection(const std::vector<std::vector<long double>>& coefficient, long double scale, unsigned degree, unsigned order, long double derivativeBound, bool computeDelta, bool fullBoundaryDelta, long double absoluteTolerance, long double relativeTolerance, unsigned maximumDepth, unsigned maximumEvaluations, const std::vector<unsigned>* zeroBoundary, std::vector<long double>& values, std::vector<long double>& errors, unsigned& evaluations, unsigned& accepted, bool useTaylor, const std::vector<std::vector<std::vector<long double>>>* surfaceGradient = nullptr) {
  bool zeroEdge = false;
  if(zeroBoundary != nullptr) {
    for(unsigned side = 0; side < 4; ++side) zeroEdge = zeroEdge || (*zeroBoundary)[side] != 0;
  }
  bool useLine = !useTaylor || zeroEdge;
  if(!useLine) {
    const unsigned totalCount = (degree+1U)*(degree+2U);
    values.assign(totalCount,0.0L);
    errors.assign(totalCount,0.0L);
    bool completed = false;
    if(derivativeBound > 0.0L && !fullBoundaryDelta) {
      try { completed = IntegrateCanonicalGraph(coefficient,degree,order,false,computeDelta,derivativeBound,absoluteTolerance,absoluteTolerance*scale,relativeTolerance,maximumDepth,values,errors,accepted,surfaceGradient); }
      catch(const std::runtime_error&) { completed = false; }
    }
    if(!completed) {
      accepted = 0;
      const std::vector<long double> domain = {0.0L,1.0L,0.0L,1.0L};
      IntegrateNormalizedCell(coefficient,degree,0,domain,true,computeDelta,fullBoundaryDelta,absoluteTolerance,absoluteTolerance*scale,relativeTolerance,32,maximumDepth,order,values,errors,accepted,surfaceGradient);
    }
    const unsigned count = (degree+1U)*(degree+2U)/2U;
    if(computeDelta && values[count] == 0.0L) useLine = true;
  }
  if(useLine) {
    accepted = 0;
    LineSliceDirection(coefficient,degree,false,order,computeDelta,fullBoundaryDelta,absoluteTolerance,absoluteTolerance*scale,relativeTolerance,maximumDepth,maximumEvaluations,values,errors,evaluations,accepted,derivativeBound > 0.0L ? order : 0U,zeroBoundary,surfaceGradient);
  }
  const unsigned count = (degree+1U)*(degree+2U)/2U;
  if(useLine && order == 2) {
    for(unsigned i = 0; i <= degree; ++i) {
      for(unsigned j = i+1U; j <= degree-i; ++j) {
        const unsigned first = LineSliceUtilities::MomentIndex(degree,i,j);
        const unsigned second = LineSliceUtilities::MomentIndex(degree,j,i);
        std::swap(values[first],values[second]);
        std::swap(values[count+first],values[count+second]);
        std::swap(errors[first],errors[second]);
        std::swap(errors[count+first],errors[count+second]);
      }
    }
  }
  for(unsigned k = 0; k < values.size(); ++k) {
    const long double tolerance = absoluteTolerance*(k < count ? 1.0L : scale)+relativeTolerance*std::abs(values[k]);
    if(!std::isfinite(values[k]) || !std::isfinite(errors[k]) || errors[k] < 0.0L || errors[k] > tolerance) throw std::runtime_error("QUAD9 accumulated directional error exceeds tolerance.");
  }
}

inline void IntegrateQuad9UnitSelected(const std::vector<std::vector<double>>& coefficient, unsigned degree, std::vector<std::vector<double>>& positive, std::vector<std::vector<double>>& delta, double absoluteTolerance = 1.0e-12, double relativeTolerance = 1.0e-12, bool computeDelta = true, bool fullBoundaryDelta = false, unsigned maximumDepth = 28, unsigned integrationOrder = 0, unsigned maximumEvaluations = 200000, std::vector<double>* estimatedError = nullptr, unsigned validationMode = 1, const std::vector<unsigned>* zeroBoundary = nullptr, bool useTaylor = false, bool referenceSurface = false, const std::vector<std::vector<std::vector<long double>>>* surfaceGradient = nullptr) {
  if(&positive == &delta || degree > 64 || !(absoluteTolerance > 0.0) || !(relativeTolerance >= 0.0) || !std::isfinite(absoluteTolerance) || !std::isfinite(relativeTolerance) || maximumDepth > 60 || integrationOrder > 2 || validationMode > 2 || maximumEvaluations < 15) throw std::invalid_argument("Invalid QUAD9 selected-direction input.");
  if(zeroBoundary != nullptr && zeroBoundary->size() != 4) throw std::invalid_argument("QUAD9 requires four boundary flags.");
  std::vector<std::vector<long double>> normalized;
  long double scale;
  LineSliceUtilities::NormalizeCanonicalPolynomial(coefficient,normalized,scale);
  std::vector<std::vector<std::vector<long double>>> referenceGradient;
  if(referenceSurface && surfaceGradient == nullptr) {
    SurfaceGradientFromPolynomial(coefficient,0.5L,0.5L,referenceGradient);
    surfaceGradient = &referenceGradient;
  }
  long double norm = 0.0L;
  for(unsigned i = 0; i < 3; ++i) {
    for(unsigned j = 0; j < 3; ++j) norm += std::abs(normalized[i][j]);
  }
  if(norm == 0.0L && computeDelta) throw std::domain_error("The Dirac delta of an identically zero QUAD9 field is undefined.");
  const unsigned count = (degree+1U)*(degree+2U)/2U;
  std::vector<long double> values(2U*count,0.0L);
  std::vector<long double> errors(2U*count,0.0L);
  int phase = 0;
  LineSliceUtilities::CanonicalGraphOrder(normalized,0,phase);
  if(computeDelta && zeroBoundary != nullptr) {
    for(unsigned side = 0; side < 4; ++side) {
      if((*zeroBoundary)[side] != 0) phase = 0;
    }
  }
  unsigned evaluations = 0;
  unsigned accepted = 0;
  unsigned selected = 0;
  unsigned compared = 0;
  unsigned attempts = 0;
  if(phase != 0 || norm == 0.0L) {
    if(phase > 0 || norm == 0.0L) {
      for(unsigned i = 0; i <= degree; ++i) {
        for(unsigned j = 0; j <= degree-i; ++j) values[LineSliceUtilities::MomentIndex(degree,i,j)] = (norm == 0.0L ? 0.5L : 1.0L)/((i+1.0L)*(j+1.0L));
      }
    }
  }
  else {
    long double lower[2];
    long double upper[2];
    const unsigned automatic = Quad9SelectedOrder(normalized,lower,upper);
    const unsigned preferred = integrationOrder == 0 ? automatic : integrationOrder;
    const bool checkBoth = validationMode == 2 || (validationMode == 1 && computeDelta && lower[2U-preferred] == 0.0L);
    bool completed = false;
    std::string failure = "QUAD9 integration directions did not agree.";
    for(unsigned retry = 0; retry < (checkBoth ? 2U : 1U) && !completed; ++retry) {
      const long double factor = retry == 0 ? 1.0L : 0.125L;
      std::vector<std::vector<long double>> trial(2);
      std::vector<std::vector<long double>> trialError(2);
      unsigned trialAccepted[2] = {0,0};
      bool success[2] = {false,false};
      const unsigned number = checkBoth || integrationOrder == 0 ? 2U : 1U;
      for(unsigned candidate = 0; candidate < number; ++candidate) {
        const unsigned order = candidate == 0 ? preferred : 3U-preferred;
        ++attempts;
        try {
          Quad9SelectedOneDirection(normalized,scale,degree,order,lower[2U-order],computeDelta,fullBoundaryDelta,absoluteTolerance*factor,relativeTolerance*factor,maximumDepth,maximumEvaluations,zeroBoundary,trial[candidate],trialError[candidate],evaluations,trialAccepted[candidate],useTaylor,surfaceGradient);
          success[candidate] = true;
          if(!checkBoth) {
            values.swap(trial[candidate]);
            errors.swap(trialError[candidate]);
            accepted = trialAccepted[candidate];
            selected = order;
            completed = true;
            break;
          }
        }
        catch(const std::runtime_error& exception) {
          failure = exception.what();
        }
      }
      if(checkBoth && success[0] && success[1]) {
        bool agree = true;
        for(unsigned k = 0; k < values.size(); ++k) {
          const long double difference = std::abs(trial[0][k]-trial[1][k]);
          const long double error = std::max(difference,std::max(trialError[0][k],trialError[1][k]));
          const long double tolerance = absoluteTolerance*(k < count ? 1.0L : scale)+relativeTolerance*std::abs(trial[0][k]);
          agree = agree && error <= tolerance;
          errors[k] = error;
        }
        if(agree) {
          values.swap(trial[0]);
          accepted = trialAccepted[0]+trialAccepted[1];
          selected = preferred;
          compared = 1;
          completed = true;
        }
        else failure = "QUAD9 x-first and y-first moments disagree beyond tolerance.";
      }
    }
    if(!completed) throw std::runtime_error(failure);
  }
  std::vector<std::vector<double>> resultH(degree+1U);
  std::vector<std::vector<double>> resultD(degree+1U);
  long double errorH = 0.0L;
  long double errorD = 0.0L;
  for(unsigned i = 0; i <= degree; ++i) {
    resultH[i].resize(degree-i+1U);
    resultD[i].resize(degree-i+1U);
    for(unsigned j = 0; j <= degree-i; ++j) {
      const unsigned k = LineSliceUtilities::MomentIndex(degree,i,j);
      resultH[i][j] = static_cast<double>(values[k]);
      resultD[i][j] = computeDelta ? static_cast<double>(values[count+k]/scale) : 0.0;
      if(!std::isfinite(resultH[i][j]) || !std::isfinite(resultD[i][j])) throw std::overflow_error("A QUAD9 selected-direction moment exceeds double range.");
      errorH = std::max(errorH,errors[k]);
      if(computeDelta) errorD = std::max(errorD,errors[count+k]/scale);
    }
  }
  positive.swap(resultH);
  delta.swap(resultD);
  if(estimatedError != nullptr) *estimatedError = {static_cast<double>(errorH),static_cast<double>(errorD),static_cast<double>(accepted),useTaylor ? -1.0 : static_cast<double>(evaluations),static_cast<double>(selected),static_cast<double>(compared),static_cast<double>(attempts)};
}

inline void IntegrateQuad9NodalSelected(const std::vector<double>& nodalValues, unsigned degree, std::vector<std::vector<double>>& positive, std::vector<std::vector<double>>& delta, double absoluteTolerance = 1.0e-12, double relativeTolerance = 1.0e-12, bool computeDelta = true, bool fullBoundaryDelta = false, unsigned validationMode = 1, unsigned integrationOrder = 0, std::vector<double>* estimatedError = nullptr, bool useTaylor = false, bool referenceSurface = false) {
  std::vector<std::vector<double>> coefficient;
  LineSliceUtilities::Quad9UnitCoefficientsFromNodalValues(nodalValues,coefficient);
  std::vector<unsigned> boundary = {static_cast<unsigned>(nodalValues[0] == 0.0 && nodalValues[3] == 0.0 && nodalValues[7] == 0.0),static_cast<unsigned>(nodalValues[1] == 0.0 && nodalValues[2] == 0.0 && nodalValues[5] == 0.0),static_cast<unsigned>(nodalValues[0] == 0.0 && nodalValues[1] == 0.0 && nodalValues[4] == 0.0),static_cast<unsigned>(nodalValues[3] == 0.0 && nodalValues[2] == 0.0 && nodalValues[6] == 0.0)};
  IntegrateQuad9UnitSelected(coefficient,degree,positive,delta,absoluteTolerance,relativeTolerance,computeDelta,fullBoundaryDelta,useTaylor ? 24U : 28U,integrationOrder,200000,estimatedError,validationMode,&boundary,useTaylor,referenceSurface);
}

inline void IntegrateQuad9CanonicalSelected(const std::vector<std::vector<double>>& coefficient, unsigned degree, std::vector<std::vector<double>>& positive, std::vector<std::vector<double>>& delta, double absoluteTolerance = 1.0e-12, double relativeTolerance = 1.0e-12, double xLeft = -1.0, double xRight = 1.0, double yBottom = -1.0, double yTop = 1.0, bool computeDelta = true, bool fullBoundaryDelta = false, unsigned validationMode = 1, unsigned integrationOrder = 0, std::vector<double>* estimatedError = nullptr, bool useTaylor = false, bool referenceSurface = false) {
  std::vector<std::vector<double>> transformed;
  LineSliceUtilities::Quad9CoefficientsToUnitSquare(coefficient,transformed,xLeft,xRight,yBottom,yTop);
  std::vector<std::vector<std::vector<long double>>> referenceGradient;
  if(referenceSurface) SurfaceGradientFromPolynomial(transformed,1.0L/(xRight-xLeft),1.0L/(yTop-yBottom),referenceGradient);
  IntegrateQuad9UnitSelected(transformed,degree,positive,delta,absoluteTolerance,relativeTolerance,computeDelta,fullBoundaryDelta,useTaylor ? 24U : 28U,integrationOrder,200000,estimatedError,validationMode,nullptr,useTaylor,false,referenceSurface ? &referenceGradient : nullptr);
}

}

#endif

