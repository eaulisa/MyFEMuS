#ifndef FEMUS_TRI7_DIRECTION_SELECTION_HPP
#define FEMUS_TRI7_DIRECTION_SELECTION_HPP

#include "LineIntegration.hpp"
#include "TriangleIntegration.hpp"
#include "TangencyStripIntegration.hpp"

// TRI7 uses ordered triangular charts internally and returns moments on the standard FEMuS reference triangle.
namespace FELevelSetMoments {

inline void Tri7SelectedDerivativeBounds(const std::vector<std::vector<long double>>& coefficient, unsigned direction, long double& lower, long double& upper) {
  if(direction < 1 || direction > 3 || coefficient.size() != 3) throw std::invalid_argument("Invalid TRI7 derivative direction.");
  long double norm = 0.0L;
  for(unsigned i = 0; i < 3; ++i) {
    if(coefficient[i].size() != 3) throw std::invalid_argument("TRI7 requires a 3 by 3 coefficient tensor.");
    for(unsigned j = 0; j < 3; ++j) {
      if(!std::isfinite(coefficient[i][j])) throw std::invalid_argument("A TRI7 coefficient is not finite.");
      norm += std::abs(coefficient[i][j]);
    }
  }
  const long double c1 = coefficient[1][0], c2 = coefficient[0][1], c3 = coefficient[2][0], c4 = coefficient[1][1], c5 = coefficient[0][2], c6 = coefficient[2][1];
  const long double dx[6] = {c1,c1+c3,c1+c3+0.5L*c4,c1+2.0L*c3,c1+2.0L*c3+0.5L*c4+c6,c1+2.0L*c3+c4+c6};
  const long double dy[6] = {c2,c2+0.5L*c4,c2+0.5L*c4+c5,c2+c4+c6,c2+c4+c5,c2+c4+2.0L*c5-c6};
  long double minimum = std::numeric_limits<long double>::infinity();
  long double maximum = -minimum;
  for(unsigned k = 0; k < 6; ++k) {
    const long double value = direction == 1 ? dy[k] : (direction == 2 ? dx[k] : dx[k]+dy[k]);
    minimum = std::min(minimum,value);
    maximum = std::max(maximum,value);
  }
  const long double padding = 2048.0L*std::numeric_limits<long double>::epsilon()*norm;
  minimum -= padding;
  maximum += padding;
  lower = minimum > 0.0L ? minimum : (maximum < 0.0L ? -maximum : 0.0L);
  upper = std::max(std::abs(minimum),std::abs(maximum));
  if(lower <= 32.0L*std::sqrt(std::numeric_limits<long double>::epsilon())*upper) lower = 0.0L;
}

inline int Tri7SelectedPhase(const std::vector<std::vector<long double>>& c, bool computeDelta) {
  const long double power[4][4] = {{c[0][0],c[1][0]+c[0][1],c[2][0]+c[1][1]+c[0][2],0.0L},{c[1][0],2.0L*c[2][0]+c[1][1],c[2][1],0.0L},{c[2][0],c[2][1],0.0L,0.0L},{0.0L,0.0L,0.0L,0.0L}};
  long double minimum = std::numeric_limits<long double>::infinity();
  long double maximum = -minimum;
  long double norm = 0.0L;
  for(unsigned a = 0; a <= 3; ++a) {
    for(unsigned b = 0; b <= 3-a; ++b) norm += std::abs(power[a][b]);
  }
  for(unsigned i = 0; i <= 3; ++i) {
    for(unsigned j = 0; j <= 3-i; ++j) {
      long double value = 0.0L;
      for(unsigned a = 0; a <= i; ++a) {
        for(unsigned b = 0; b <= j; ++b) value += power[a][b]*LineSliceUtilities::Binomial(i,a)*LineSliceUtilities::Binomial(j,b)/(LineSliceUtilities::Binomial(3,a+b)*LineSliceUtilities::Binomial(a+b,a));
      }
      minimum = std::min(minimum,value);
      maximum = std::max(maximum,value);
    }
  }
  const long double padding = 1024.0L*std::numeric_limits<long double>::epsilon()*norm;
  if(minimum > padding || (!computeDelta && minimum >= 0.0L && maximum > 0.0L)) return 1;
  if(maximum < -padding || (!computeDelta && maximum <= 0.0L && minimum < 0.0L)) return -1;
  return 0;
}

inline void Tri7SelectedDiagonalPolynomial(const std::vector<std::vector<long double>>& input, std::vector<std::vector<long double>>& output) {
  const long double c0 = input[0][0], c1 = input[1][0], c2 = input[0][1], c3 = input[2][0], c4 = input[1][1], c5 = input[0][2], c6 = input[2][1];
  output.assign(3,std::vector<long double>(3,0.0L));
  output[0][0] = c0+c1+c3;
  output[1][0] = c2+c4+c6;
  output[0][1] = -c1-c2-2.0L*c3-c4-c6;
  output[2][0] = c5-c6;
  output[1][1] = -c4-2.0L*c5;
  output[0][2] = c3+c4+c5+c6;
  output[2][1] = c6;
  output[1][2] = -c6;
}

inline void Tri7SelectedRestoreDiagonal(unsigned degree, const std::vector<long double>& raw, const std::vector<long double>& rawError, std::vector<long double>& values, std::vector<long double>& errors) {
  const unsigned count = (degree+1U)*(degree+2U)/2U;
  values.assign(2U*count,0.0L);
  errors.assign(2U*count,0.0L);
  for(unsigned i = 0; i <= degree; ++i) {
    for(unsigned j = 0; j <= degree-i; ++j) {
      const unsigned target = LineSliceUtilities::MomentIndex(degree,i,j);
      long double correction[2] = {0.0L,0.0L};
      long double absoluteSum[2] = {0.0L,0.0L};
      for(unsigned a = 0; a <= i; ++a) {
        for(unsigned b = 0; b <= j; ++b) {
          const unsigned source = LineSliceUtilities::MomentIndex(degree,j-b,a+b);
          long double multiplier = LineSliceUtilities::Binomial(i,a)*LineSliceUtilities::Binomial(j,b);
          if((a+b)%2U != 0) multiplier = -multiplier;
          for(unsigned pair = 0; pair < 2; ++pair) {
            const unsigned destination = pair*count+target;
            const long double term = multiplier*raw[pair*count+source];
            const long double adjusted = term-correction[pair];
            const long double next = values[destination]+adjusted;
            correction[pair] = (next-values[destination])-adjusted;
            values[destination] = next;
            errors[destination] += std::abs(multiplier)*rawError[pair*count+source];
            absoluteSum[pair] += std::abs(term);
          }
        }
      }
      for(unsigned pair = 0; pair < 2; ++pair) errors[pair*count+target] += (8.0L+2.0L*(i+1U)*(j+1U))*std::numeric_limits<long double>::epsilon()*absoluteSum[pair];
    }
  }
}

inline void Tri7SelectedOneDirection(const std::vector<std::vector<long double>>& normalized, long double scale, unsigned degree, unsigned direction, long double derivativeBound, bool computeDelta, bool fullBoundaryDelta, long double absoluteTolerance, long double relativeTolerance, unsigned maximumDepth, unsigned maximumEvaluations, const std::vector<unsigned>* zeroBoundary, bool useTaylor, std::vector<long double>& values, std::vector<long double>& errors, unsigned& evaluations, unsigned& accepted, const std::vector<std::vector<std::vector<long double>>>* surfaceGradient = nullptr, bool fastOnly = false) {
  const unsigned count = (degree+1U)*(degree+2U)/2U;
  const unsigned order = direction == 2 ? 2U : 1U;
  const long double amplification = direction == 3 ? LineSliceUtilities::Power(2.0L,degree) : 1.0L;
  std::vector<std::vector<long double>> coefficient;
  if(direction == 3) Tri7SelectedDiagonalPolynomial(normalized,coefficient);
  else coefficient = normalized;
  std::vector<std::vector<std::vector<long double>>> diagonalGradient;
  if(direction == 3 && surfaceGradient != nullptr) {
    SurfaceGradientDiagonal(*surfaceGradient,diagonalGradient);
    surfaceGradient = &diagonalGradient;
  }
  std::vector<unsigned> boundary;
  bool zeroEdge = false;
  if(zeroBoundary != nullptr) {
    boundary = *zeroBoundary;
    if(direction == 3) boundary = {(*zeroBoundary)[1],(*zeroBoundary)[2],(*zeroBoundary)[0]};
    for(unsigned k = 0; k < 3; ++k) zeroEdge = zeroEdge || boundary[k] != 0;
  }
  std::vector<long double> raw(2U*count,0.0L), rawError(2U*count,0.0L);
  const long double absoluteH = absoluteTolerance/amplification;
  const long double absoluteD = (absoluteTolerance*scale)/amplification;
  const long double relative = relativeTolerance/amplification;
  bool line = !useTaylor || zeroEdge;
  if(!line) {
    bool completed = false;
    if(derivativeBound > 0.0L && !fullBoundaryDelta) {
      try { completed = LineIntegration::IntegrateFastGraph(coefficient,degree,order,true,computeDelta,derivativeBound,absoluteH,absoluteD,relative,maximumDepth,raw,rawError,accepted,surfaceGradient); }
      catch(const std::runtime_error&) { completed = false; }
    }
    if(!completed && fastOnly) throw std::runtime_error("TRI graph integration requires a parent-coordinate chart.");
    if(!completed) {
      accepted = 0;
      const std::vector<long double> domain = {0.0L,1.0L,0.0L,1.0L};
      try { LineIntegration::IntegrateTaylor(coefficient,degree,2,domain,true,computeDelta,fullBoundaryDelta,absoluteH,absoluteD,relative,32,maximumDepth,order,raw,rawError,accepted,surfaceGradient); }
      catch(const std::runtime_error&) { line = true; }
    }
    if(computeDelta && raw[count] == 0.0L) line = true;
  }
  if(fastOnly && (line || (computeDelta && raw[count] == 0.0L))) throw std::runtime_error("TRI graph integration requires a parent-coordinate chart.");
  if(line) {
    accepted = 0;
    LineIntegration::IntegrateFullGaussKronrod(coefficient,degree,true,order,computeDelta,fullBoundaryDelta,absoluteH,absoluteD,relative,maximumDepth,maximumEvaluations,raw,rawError,evaluations,accepted,derivativeBound > 0.0L ? order : 0U,zeroBoundary == nullptr ? nullptr : &boundary,surfaceGradient);
    if(order == 2) {
      for(unsigned i = 0; i <= degree; ++i) {
        for(unsigned j = i+1U; j <= degree-i; ++j) {
          const unsigned a = LineSliceUtilities::MomentIndex(degree,i,j), b = LineSliceUtilities::MomentIndex(degree,j,i);
          for(unsigned pair = 0; pair < 2; ++pair) {
            std::swap(raw[pair*count+a],raw[pair*count+b]);
            std::swap(rawError[pair*count+a],rawError[pair*count+b]);
          }
        }
      }
    }
  }
  if(direction == 3) Tri7SelectedRestoreDiagonal(degree,raw,rawError,values,errors);
  else {
    values.swap(raw);
    errors.swap(rawError);
  }
  for(unsigned k = 0; k < values.size(); ++k) {
    const long double tolerance = absoluteTolerance*(k < count ? 1.0L : scale)+relativeTolerance*std::abs(values[k]);
    if(!std::isfinite(values[k]) || !std::isfinite(errors[k]) || errors[k] < 0.0L || errors[k] > tolerance) throw std::runtime_error("TRI7 transformed directional error exceeds tolerance.");
  }
}

inline void Tri7SurfaceCoordinatePolynomials(const std::vector<std::vector<double>>& geometryCoordinates, std::vector<std::vector<std::vector<long double>>>& coordinatePolynomials) {
  ValidateSurfaceGeometryCoordinates(geometryCoordinates,2,7);
  coordinatePolynomials.assign(2,std::vector<std::vector<long double>>(3,std::vector<long double>(3,0.0L)));
  for(unsigned component = 0; component < 2; ++component) {
    std::vector<double> shifted(7);
    for(unsigned node = 0; node < 7; ++node) shifted[node] = geometryCoordinates[component][node]-geometryCoordinates[component][0];
    std::vector<std::vector<double>> coefficient;
    Tri7ReferenceCoefficientsFromNodalValues(shifted,coefficient);
    for(unsigned i = 0; i < 3; ++i) {
      for(unsigned j = 0; j < 3; ++j) coordinatePolynomials[component][i][j] = coefficient[i][j];
    }
  }
}

inline void IntegrateTri7ReferenceSelected(const std::vector<std::vector<double>>& coefficient, unsigned degree, std::vector<std::vector<double>>& positive, std::vector<std::vector<double>>& delta, double absoluteTolerance = 1.0e-12, double relativeTolerance = 1.0e-12, bool computeDelta = true, bool fullBoundaryDelta = false, unsigned maximumDepth = 28, unsigned integrationDirection = 0, unsigned maximumEvaluations = 200000, std::vector<double>* estimatedError = nullptr, unsigned validationMode = 1, const std::vector<unsigned>* zeroBoundary = nullptr, bool useTaylor = false, bool referenceSurface = false, const std::vector<std::vector<double>>* inverseJacobian = nullptr, const std::vector<std::vector<double>>* geometryCoordinates = nullptr) {
  if(geometryCoordinates != nullptr) {
    if(!referenceSurface || inverseJacobian != nullptr) throw std::invalid_argument("Surface geometry requires surface moments without an inverse Jacobian.");
    ValidateSurfaceGeometryCoordinates(*geometryCoordinates,2,7);
  }
  if(inverseJacobian != nullptr) {
    if(!referenceSurface) throw std::invalid_argument("A surface inverse Jacobian requires surface moments.");
    ValidateSurfaceInverseJacobian(*inverseJacobian,2);
  }
  if(&positive == &delta || coefficient.size() != 3 || degree > 64 || !(absoluteTolerance > 0.0) || !(relativeTolerance >= 0.0) || !std::isfinite(absoluteTolerance) || !std::isfinite(relativeTolerance) || maximumDepth > 60 || integrationDirection > 3 || validationMode > 2 || maximumEvaluations < 15) throw std::invalid_argument("Invalid TRI7 selected-direction input.");
  for(unsigned i = 0; i < 3; ++i) if(coefficient[i].size() != 3) throw std::invalid_argument("The TRI7 coefficient matrix must be 3 by 3.");
  if(coefficient[2][2] != 0.0 || coefficient[2][1] != coefficient[1][2]) throw std::invalid_argument("The TRI7 coefficients must represent a quadratic polynomial with a barycentric bubble.");
  if(zeroBoundary != nullptr && zeroBoundary->size() != 3) throw std::invalid_argument("TRI7 requires three boundary flags.");
  const long double amplification = ReferencePower(2.0L,degree);
  const long double chartAbsoluteTolerance = absoluteTolerance/amplification;
  const long double chartRelativeTolerance = relativeTolerance/amplification;
  std::vector<std::vector<long double>> chartCoefficient;
  ReferencePolynomialAffine2D(coefficient,1.0L,-1.0L,0.0L,1.0L,chartCoefficient);
  std::vector<std::vector<long double>> normalized;
  long double scale;
  LineSliceUtilities::NormalizeCanonicalPolynomial(chartCoefficient,normalized,scale);
  long double norm = 0.0L;
  for(unsigned i = 0; i < 3; ++i) {
    for(unsigned j = 0; j < 3; ++j) norm += std::abs(normalized[i][j]);
  }
  const bool zero = norm == 0.0L;
  if(zero && computeDelta) throw std::domain_error("The Dirac delta of an identically zero TRI7 field is undefined.");
  std::vector<std::vector<std::vector<long double>>> referenceGradient;
  std::vector<std::vector<std::vector<long double>>> standardGradient;
  std::vector<std::vector<std::vector<long double>>> standardCoordinates;
  if(referenceSurface) {
    SurfaceGradientFromPolynomial(coefficient,1.0L,1.0L,referenceGradient);
    if(inverseJacobian != nullptr) SurfaceGradientPhysical(referenceGradient,*inverseJacobian);
    if(geometryCoordinates != nullptr) {
      Tri7SurfaceCoordinatePolynomials(*geometryCoordinates,standardCoordinates);
      SurfaceGradientCofactor2D(referenceGradient,standardCoordinates,1.0L,1.0L);
    }
    standardGradient = referenceGradient;
    for(unsigned component = 0; component < referenceGradient.size(); ++component) {
      std::vector<std::vector<long double>> chartComponent;
      ReferencePolynomialAffine2D(referenceGradient[component],1.0L,-1.0L,0.0L,1.0L,chartComponent);
      referenceGradient[component].swap(chartComponent);
    }
  }
  const unsigned count = (degree+1U)*(degree+2U)/2U;
  std::vector<long double> values(2U*count,0.0L), errors(2U*count,0.0L);
  unsigned evaluations = 0, accepted = 0, selected = 0, compared = 0, attempts = 0;
  bool standardMoments = false;
  int phase = Tri7SelectedPhase(normalized,computeDelta);
  if(computeDelta && zeroBoundary != nullptr) {
    for(unsigned side = 0; side < 3; ++side) {
      if((*zeroBoundary)[side] != 0) phase = 0;
    }
  }
  if(phase != 0 || zero) {
    if(phase > 0 || zero) {
      for(unsigned i = 0; i <= degree; ++i) {
        for(unsigned j = 0; j <= degree-i; ++j) {
          const unsigned k = LineSliceUtilities::MomentIndex(degree,i,j);
          values[k] = (zero ? 0.5L : 1.0L)*ReferenceTriangleMoment(i,j);
          errors[k] = (4.0L+2.0L*j)*std::numeric_limits<long double>::epsilon()*std::abs(values[k]);
        }
      }
    }
  }
  else {
    long double lower[3], upper[3];
    unsigned ranked[3] = {1,2,3};
    for(unsigned k = 0; k < 3; ++k) Tri7SelectedDerivativeBounds(normalized,k+1U,lower[k],upper[k]);
    for(unsigned i = 0; i < 3; ++i) {
      for(unsigned j = i+1U; j < 3; ++j) {
        if(lower[ranked[j]-1U] > lower[ranked[i]-1U] || (lower[ranked[j]-1U] == lower[ranked[i]-1U] && upper[ranked[j]-1U] > upper[ranked[i]-1U])) std::swap(ranked[i],ranked[j]);
      }
    }
    bool useSplitIntegration = integrationDirection == 0 && !fullBoundaryDelta && lower[ranked[0]-1U] == 0.0L;
    if(zeroBoundary != nullptr) {
      for(unsigned side = 0; side < 3; ++side) useSplitIntegration = useSplitIntegration && (*zeroBoundary)[side] == 0;
    }
    bool useStandardIntegration = useTaylor && validationMode != 2 && integrationDirection == 0 && !fullBoundaryDelta;
    if(zeroBoundary != nullptr) for(unsigned side = 0; side < 3; ++side) useStandardIntegration = useStandardIntegration && (*zeroBoundary)[side] == 0;
    if(useStandardIntegration) {
      {
        std::vector<std::vector<long double>> standardCoefficient(3,std::vector<long double>(3,0.0L));
        for(unsigned i = 0; i < 3; ++i) for(unsigned j = 0; j < 3; ++j) standardCoefficient[i][j] = static_cast<long double>(coefficient[i][j])/scale;
        const std::vector<TriangleDomainVertex> vertices = {{0.0L,0.0L},{1.0L,0.0L},{0.0L,1.0L}};
        accepted = 0;
        ++attempts;
        bool tangencyIntegrated = false;
        if(lower[ranked[0]-1U] == 0 && !TangencyStripIntegration::IsGraph(standardCoefficient,true)) tangencyIntegrated = TangencyStripIntegration::Try(standardCoefficient,degree,true,computeDelta,absoluteTolerance,absoluteTolerance*scale,relativeTolerance,maximumDepth,maximumEvaluations,values,errors,evaluations,accepted,selected,referenceSurface ? &standardGradient : nullptr,scale,standardCoordinates.empty() ? nullptr : &standardCoordinates,inverseJacobian);
        if(!tangencyIntegrated) TriangleIntegration::Integrate(standardCoefficient,degree,vertices,computeDelta,absoluteTolerance,absoluteTolerance*scale,relativeTolerance,maximumDepth,maximumEvaluations,0,10,values,errors,evaluations,accepted,selected,referenceSurface ? &standardGradient : nullptr);
        standardMoments = true;
      }
    }
    else if(useSplitIntegration) {
      const std::vector<TriangleDomainVertex> vertices = {{0.0L,0.0L},{1.0L,0.0L},{1.0L,1.0L}};
      TriangleIntegration::Integrate(normalized,degree,vertices,computeDelta,chartAbsoluteTolerance,chartAbsoluteTolerance*scale,chartRelativeTolerance,maximumDepth,maximumEvaluations,0,10,values,errors,evaluations,accepted,selected,referenceSurface ? &referenceGradient : nullptr);
      if(validationMode == 2) {
        bool verified = false;
        std::string failure = "TRI split and independent direction did not agree.";
        for(unsigned candidate = 0; candidate < 3 && !verified; ++candidate) {
          std::vector<long double> independent;
          std::vector<long double> independentError;
          unsigned independentAccepted = 0;
          unsigned independentEvaluations = 0;
          ++attempts;
          try {
            Tri7SelectedOneDirection(normalized,scale,degree,ranked[candidate],lower[ranked[candidate]-1U],computeDelta,fullBoundaryDelta,chartAbsoluteTolerance,chartRelativeTolerance,maximumDepth,maximumEvaluations,zeroBoundary,false,independent,independentError,independentEvaluations,independentAccepted,referenceSurface ? &referenceGradient : nullptr);
            verified = true;
            for(unsigned k = 0; k < values.size(); ++k) {
              const long double estimate = std::max(std::abs(values[k]-independent[k]),std::max(errors[k],independentError[k]));
              const long double tolerance = chartAbsoluteTolerance*(k < count ? 1.0L : scale)+chartRelativeTolerance*std::abs(values[k]);
              if(!std::isfinite(estimate) || estimate > tolerance) verified = false;
            }
            if(verified) {
              for(unsigned k = 0; k < values.size(); ++k) errors[k] = std::max(std::abs(values[k]-independent[k]),std::max(errors[k],independentError[k]));
              compared = ranked[candidate];
              evaluations += independentEvaluations;
            }
          }
          catch(const std::runtime_error& exception) { failure = exception.what(); }
        }
        if(!verified) throw std::runtime_error(failure);
      }
    }
    else {
    if(integrationDirection != 0) {
      for(unsigned k = 0; k < 3; ++k) {
        if(ranked[k] == integrationDirection) std::swap(ranked[0],ranked[k]);
      }
    }
    const bool verify = validationMode == 2;
    bool completed = false;
    std::string failure = "TRI7 integration directions did not agree.";
    for(unsigned retry = 0; retry < (verify ? 2U : 1U) && !completed; ++retry) {
      const long double factor = retry == 0 ? 1.0L : 0.125L;
      std::vector<std::vector<long double>> trial(3), trialError(3);
      unsigned trialAccepted[3] = {0,0,0};
      bool success[3] = {false,false,false};
      const unsigned number = integrationDirection != 0 && !verify ? 1U : 3U;
      for(unsigned candidate = 0; candidate < number && !completed; ++candidate) {
        ++attempts;
        try {
          Tri7SelectedOneDirection(normalized,scale,degree,ranked[candidate],lower[ranked[candidate]-1U],computeDelta,fullBoundaryDelta,chartAbsoluteTolerance*factor,chartRelativeTolerance*factor,maximumDepth,maximumEvaluations,zeroBoundary,useTaylor,trial[candidate],trialError[candidate],evaluations,trialAccepted[candidate],referenceSurface ? &referenceGradient : nullptr);
          success[candidate] = true;
          if(!verify) {
            values.swap(trial[candidate]);
            errors.swap(trialError[candidate]);
            accepted = trialAccepted[candidate];
            selected = ranked[candidate];
            completed = true;
            break;
          }
          for(unsigned previous = 0; previous < candidate && !completed; ++previous) {
            if(!success[previous]) continue;
            bool agree = true;
            std::vector<long double> guarded(errors.size());
            for(unsigned k = 0; k < guarded.size(); ++k) {
              guarded[k] = std::max(std::abs(trial[previous][k]-trial[candidate][k]),std::max(trialError[previous][k],trialError[candidate][k]));
              const long double tolerance = chartAbsoluteTolerance*(k < count ? 1.0L : scale)+chartRelativeTolerance*std::abs(trial[previous][k]);
              agree = agree && guarded[k] <= tolerance;
            }
            if(agree) {
              values.swap(trial[previous]);
              errors.swap(guarded);
              accepted = trialAccepted[previous]+trialAccepted[candidate];
              selected = ranked[previous];
              compared = ranked[candidate];
              completed = true;
            }
          }
        }
        catch(const std::runtime_error& exception) { failure = exception.what(); }
      }
    }
    if(!completed) throw std::runtime_error("TRI7 direction selection failed: "+failure);
    }
    if(!standardMoments) {
      std::vector<long double> referenceValues;
      std::vector<long double> referenceErrors;
      ReferenceAffineMoments2D(degree,values,errors,1.0L,-1.0L,0.0L,1.0L,1.0L,referenceValues,referenceErrors);
      values.swap(referenceValues);
      errors.swap(referenceErrors);
    }
  }
  std::vector<std::vector<double>> resultH(degree+1U), resultD(degree+1U);
  long double errorH = 0.0L, errorD = 0.0L;
  for(unsigned i = 0; i <= degree; ++i) {
    resultH[i].resize(degree-i+1U);
    resultD[i].resize(degree-i+1U);
    for(unsigned j = 0; j <= degree-i; ++j) {
      const unsigned k = LineSliceUtilities::MomentIndex(degree,i,j);
      resultH[i][j] = static_cast<double>(values[k]);
      resultD[i][j] = computeDelta ? static_cast<double>(values[count+k]/scale) : 0.0;
      if(!std::isfinite(resultH[i][j]) || !std::isfinite(resultD[i][j])) throw std::overflow_error("A TRI7 selected-direction moment exceeds double range.");
      const long double momentErrorH = errors[k]+std::abs(static_cast<long double>(resultH[i][j])-values[k]);
      const long double momentErrorD = computeDelta ? errors[count+k]/scale+std::abs(static_cast<long double>(resultD[i][j])-values[count+k]/scale) : 0.0L;
      if(!std::isfinite(momentErrorH) || !std::isfinite(momentErrorD) || momentErrorH > absoluteTolerance+relativeTolerance*std::abs(static_cast<long double>(resultH[i][j])) || momentErrorD > absoluteTolerance+relativeTolerance*std::abs(static_cast<long double>(resultD[i][j]))) throw std::runtime_error("TRI7 reference moment error exceeds tolerance.");
      errorH = std::max(errorH,momentErrorH);
      errorD = std::max(errorD,momentErrorD);
    }
  }
  positive.swap(resultH);
  delta.swap(resultD);
  if(estimatedError != nullptr) *estimatedError = {static_cast<double>(errorH),static_cast<double>(errorD),static_cast<double>(accepted),useTaylor ? -1.0 : static_cast<double>(evaluations),static_cast<double>(selected),static_cast<double>(compared),static_cast<double>(attempts)};
}

inline void IntegrateTri7NodalSelected(const std::vector<double>& nodalValues, unsigned degree, std::vector<std::vector<double>>& positive, std::vector<std::vector<double>>& delta, double absoluteTolerance = 1.0e-12, double relativeTolerance = 1.0e-12, bool computeDelta = true, bool fullBoundaryDelta = false, unsigned validationMode = 1, unsigned integrationDirection = 0, std::vector<double>* estimatedError = nullptr, bool useTaylor = false, bool referenceSurface = false, const std::vector<std::vector<double>>* inverseJacobian = nullptr, const std::vector<std::vector<double>>* geometryCoordinates = nullptr) {
  std::vector<std::vector<double>> coefficient;
  Tri7ReferenceCoefficientsFromNodalValues(nodalValues,coefficient);
  const std::vector<unsigned> boundary = {static_cast<unsigned>(nodalValues[0] == 0.0 && nodalValues[1] == 0.0 && nodalValues[3] == 0.0),static_cast<unsigned>(nodalValues[0] == 0.0 && nodalValues[2] == 0.0 && nodalValues[5] == 0.0),static_cast<unsigned>(nodalValues[1] == 0.0 && nodalValues[2] == 0.0 && nodalValues[4] == 0.0)};
  IntegrateTri7ReferenceSelected(coefficient,degree,positive,delta,absoluteTolerance,relativeTolerance,computeDelta,fullBoundaryDelta,useTaylor ? 24U : 28U,integrationDirection,200000,estimatedError,validationMode,&boundary,useTaylor,referenceSurface,inverseJacobian,geometryCoordinates);
}

}

#endif

