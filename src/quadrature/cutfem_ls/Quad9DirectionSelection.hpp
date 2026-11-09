#ifndef FEMUS_QUAD9_DIRECTION_SELECTION_HPP
#define FEMUS_QUAD9_DIRECTION_SELECTION_HPP

#include "LineIntegration.hpp"

// QUAD9 uses square charts internally and returns moments on the standard FEMuS reference square.
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

inline unsigned Quad9SelectedWeakGraphOrder(const std::vector<std::vector<long double>>& coefficient, unsigned requestedOrder = 0) {
  long double bernstein[3][3];
  CanonicalBernstein(coefficient,bernstein);
  long double norm = 0.0L;
  for(unsigned i = 0; i < 3; ++i) for(unsigned j = 0; j < 3; ++j) norm += std::abs(coefficient[i][j]);
  const long double padding = 1024.0L*std::numeric_limits<long double>::epsilon()*std::max(1.0L,norm);
  unsigned selected = 0;
  long double selectedMagnitude = 0.0L;
  for(unsigned order = 1; order <= 2; ++order) {
    if(requestedOrder != 0 && requestedOrder != order) continue;
    long double minimum = std::numeric_limits<long double>::infinity();
    long double maximum = -minimum;
    for(unsigned i = 0; i < 3; ++i) {
      const long double low = order == 1 ? bernstein[i][0] : bernstein[0][i];
      const long double middle = order == 1 ? bernstein[i][1] : bernstein[1][i];
      const long double high = order == 1 ? bernstein[i][2] : bernstein[2][i];
      const long double derivative[3] = {2.0L*(middle-low),high-low,2.0L*(high-middle)};
      for(unsigned j = 0; j < 3; ++j) {
        minimum = std::min(minimum,derivative[j]);
        maximum = std::max(maximum,derivative[j]);
      }
    }
    const bool nonnegative = minimum >= -padding && maximum > padding;
    const bool nonpositive = maximum <= padding && minimum < -padding;
    if(nonnegative || nonpositive) {
      const long double magnitude = std::max(std::abs(minimum),std::abs(maximum));
      if(magnitude > selectedMagnitude) {
        selectedMagnitude = magnitude;
        selected = order;
      }
    }
  }
  return selected;
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
      try { completed = LineIntegration::IntegrateFastGraph(coefficient,degree,order,false,computeDelta,derivativeBound,absoluteTolerance,absoluteTolerance*scale,relativeTolerance,maximumDepth,values,errors,accepted,surfaceGradient); }
      catch(const std::runtime_error&) { completed = false; }
    }
    if(!completed) {
      accepted = 0;
      const std::vector<long double> domain = {0.0L,1.0L,0.0L,1.0L};
      LineIntegration::IntegrateTaylor(coefficient,degree,0,domain,true,computeDelta,fullBoundaryDelta,absoluteTolerance,absoluteTolerance*scale,relativeTolerance,32,maximumDepth,order,values,errors,accepted,surfaceGradient);
    }
    const unsigned count = (degree+1U)*(degree+2U)/2U;
    if(computeDelta && values[count] == 0.0L) useLine = true;
  }
  if(useLine) {
    accepted = 0;
    LineIntegration::IntegrateFullGaussKronrod(coefficient,degree,false,order,computeDelta,fullBoundaryDelta,absoluteTolerance,absoluteTolerance*scale,relativeTolerance,maximumDepth,maximumEvaluations,values,errors,evaluations,accepted,derivativeBound > 0.0L ? order : 0U,zeroBoundary,surfaceGradient);
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

inline int Quad9SelectedReflectionParity(const std::vector<std::vector<long double>>& coefficient, unsigned axis) {
  unsigned columns = 0;
  for(unsigned row = 0; row < coefficient.size(); ++row) columns = std::max(columns,static_cast<unsigned>(coefficient[row].size()));
  const unsigned longitudinal = axis == 0 ? coefficient.size() : columns;
  const unsigned transverse = axis == 0 ? columns : coefficient.size();
  bool even = true;
  bool odd = true;
  for(unsigned k = 0; k < transverse; ++k) {
    long double polynomial[3] = {0.0L,0.0L,0.0L};
    for(unsigned power = 0; power < longitudinal; ++power) {
      const unsigned row = axis == 0 ? power : k;
      const unsigned column = axis == 0 ? k : power;
      const long double value = column < coefficient[row].size() ? coefficient[row][column] : 0.0L;
      if(power > 2 && value != 0.0L) return 0;
      if(power <= 2) polynomial[power] = value;
    }
    even = even && polynomial[1] == -polynomial[2];
    odd = odd && polynomial[2] == 0.0L && polynomial[1] == -2.0L*polynomial[0];
  }
  return even ? 1 : (odd ? -1 : 0);
}

inline void Quad9SelectedReferenceParity(const std::vector<std::vector<long double>>& coefficient, unsigned degree, std::vector<long double>& values, std::vector<long double>& errors, const std::vector<std::vector<std::vector<long double>>>* surfaceGradient) {
  const unsigned count = (degree+1U)*(degree+2U)/2U;
  for(unsigned axis = 0; axis < 2; ++axis) {
    const int parity = Quad9SelectedReflectionParity(coefficient,axis);
    if(parity == 0) continue;
    bool deltaSymmetric = true;
    if(surfaceGradient != nullptr) {
      for(unsigned component = 0; component < surfaceGradient->size(); ++component) deltaSymmetric = deltaSymmetric && Quad9SelectedReflectionParity((*surfaceGradient)[component],axis) != 0;
    }
    for(unsigned i = 0; i <= degree; ++i) {
      for(unsigned j = 0; j <= degree-i; ++j) {
        if((axis == 0 ? i : j)%2U == 0) continue;
        const unsigned index = LineSliceUtilities::MomentIndex(degree,i,j);
        if(parity == 1) {
          values[index] = 0.0L;
          errors[index] = 0.0L;
        }
        if(deltaSymmetric) {
          values[count+index] = 0.0L;
          errors[count+index] = 0.0L;
        }
      }
    }
  }
}

inline bool Quad9SelectedDerivativeSplitCoordinate(const std::vector<std::vector<long double>>& coefficient, unsigned axis, long double& coordinate) {
  std::vector<long double> roots;
  for(unsigned side = 0; side < 2; ++side) {
    long double constant = 0.0L;
    long double linear = 0.0L;
    if(axis == 0) {
      long double sidePower = 1.0L;
      for(unsigned j = 0; j < 3; ++j) {
        constant += coefficient[1][j]*sidePower;
        linear += 2.0L*coefficient[2][j]*sidePower;
        sidePower *= static_cast<long double>(side);
      }
    }
    else {
      long double sidePower = 1.0L;
      for(unsigned i = 0; i < 3; ++i) {
        constant += coefficient[i][1]*sidePower;
        linear += 2.0L*coefficient[i][2]*sidePower;
        sidePower *= static_cast<long double>(side);
      }
    }
    const long double scale = std::abs(constant)+std::abs(linear);
    const long double padding = 1024.0L*std::numeric_limits<long double>::epsilon()*std::max(1.0L,scale);
    if(std::abs(linear) <= padding) continue;
    const long double root = -constant/linear;
    if(root > 0.02L && root < 0.98L) roots.push_back(root);
  }
  if(roots.empty()) return false;
  coordinate = 0.0L;
  for(unsigned k = 0; k < roots.size(); ++k) coordinate += roots[k];
  coordinate /= static_cast<long double>(roots.size());
  long double spread = 0.0L;
  for(unsigned k = 0; k < roots.size(); ++k) spread = std::max(spread,std::abs(roots[k]-coordinate));
  if(roots.size() > 1 && spread > 0.30L) return false;
  coordinate = std::max(0.05L,std::min(0.95L,coordinate));
  return true;
}

inline void Quad9SelectedSplitCanonical(const std::vector<std::vector<long double>>& coefficient, long double scale, unsigned degree, bool computeDelta, long double absoluteTolerance, long double relativeTolerance, unsigned maximumDepth, unsigned maximumEvaluations, const std::vector<unsigned>* zeroBoundary, bool useTaylor, const std::vector<std::vector<std::vector<long double>>>* surfaceGradient, unsigned splitDepth, unsigned maximumSplitDepth, std::vector<long double>& values, std::vector<long double>& errors, unsigned& evaluations, unsigned& accepted) {
  const unsigned count = (degree+1U)*(degree+2U)/2U;
  values.assign(2U*count,0.0L);
  errors.assign(2U*count,0.0L);
  std::vector<long double> correction(2U*count,0.0L);
  long double splitCoordinate[2] = {0.5L,0.5L};
  const bool splitAvailable[2] = {Quad9SelectedDerivativeSplitCoordinate(coefficient,0,splitCoordinate[0]),Quad9SelectedDerivativeSplitCoordinate(coefficient,1,splitCoordinate[1])};
  long double splitLower[2] = {};
  long double splitUpper[2] = {};
  const unsigned preferredOrder = Quad9SelectedOrder(coefficient,splitLower,splitUpper);
  const unsigned preferredAxis = preferredOrder == 2 ? 0U : 1U;
  unsigned splitAxis = splitAvailable[preferredAxis] ? preferredAxis : (splitAvailable[1U-preferredAxis] ? 1U-preferredAxis : 2U);
  const unsigned firstCount = splitAxis == 0 || splitAxis == 2 ? 2U : 1U;
  const unsigned secondCount = splitAxis == 1 || splitAxis == 2 ? 2U : 1U;
  for(unsigned firstSide = 0; firstSide < firstCount; ++firstSide) {
    const long double firstOrigin = firstCount == 1 ? 0.0L : (firstSide == 0 ? 0.0L : splitCoordinate[0]);
    const long double firstScale = firstCount == 1 ? 1.0L : (firstSide == 0 ? splitCoordinate[0] : 1.0L-splitCoordinate[0]);
    for(unsigned secondSide = 0; secondSide < secondCount; ++secondSide) {
      const long double secondOrigin = secondCount == 1 ? 0.0L : (secondSide == 0 ? 0.0L : splitCoordinate[1]);
      const long double secondScale = secondCount == 1 ? 1.0L : (secondSide == 0 ? splitCoordinate[1] : 1.0L-splitCoordinate[1]);
      std::vector<std::vector<long double>> childCoefficient;
      ReferencePolynomialAffine2D(coefficient,firstOrigin,firstScale,secondOrigin,secondScale,childCoefficient);
      std::vector<std::vector<long double>> childNormalized;
      long double childScale = 1.0L;
      LineSliceUtilities::NormalizeCanonicalPolynomial(childCoefficient,childNormalized,childScale);
      std::vector<std::vector<std::vector<long double>>> childGradient;
      if(surfaceGradient != nullptr) {
        childGradient.resize(surfaceGradient->size());
        for(unsigned component = 0; component < surfaceGradient->size(); ++component) ReferencePolynomialAffine2D((*surfaceGradient)[component],firstOrigin,firstScale,secondOrigin,secondScale,childGradient[component]);
      }
      std::vector<unsigned> boundary(4,1);
      for(unsigned i = 0; i < 3; ++i) {
        boundary[0] = boundary[0] != 0 && childNormalized[0][i] == 0.0L;
        boundary[1] = boundary[1] != 0 && (childNormalized[0][i]+childNormalized[1][i])+childNormalized[2][i] == 0.0L;
        boundary[2] = boundary[2] != 0 && childNormalized[i][0] == 0.0L;
        boundary[3] = boundary[3] != 0 && (childNormalized[i][0]+childNormalized[i][1])+childNormalized[i][2] == 0.0L;
      }
      if(zeroBoundary != nullptr) {
        if(firstCount == 1) {
          boundary[0] = boundary[0] != 0 || (*zeroBoundary)[0] != 0;
          boundary[1] = boundary[1] != 0 || (*zeroBoundary)[1] != 0;
        }
        else if(firstSide == 0) boundary[0] = boundary[0] != 0 || (*zeroBoundary)[0] != 0;
        else boundary[1] = boundary[1] != 0 || (*zeroBoundary)[1] != 0;
        if(secondCount == 1) {
          boundary[2] = boundary[2] != 0 || (*zeroBoundary)[2] != 0;
          boundary[3] = boundary[3] != 0 || (*zeroBoundary)[3] != 0;
        }
        else if(secondSide == 0) boundary[2] = boundary[2] != 0 || (*zeroBoundary)[2] != 0;
        else boundary[3] = boundary[3] != 0 || (*zeroBoundary)[3] != 0;
      }
      int phase = 0;
      LineSliceUtilities::CanonicalGraphOrder(childNormalized,0,phase);
      if(computeDelta) {
        for(unsigned side = 0; side < 4; ++side) if(boundary[side] != 0) phase = 0;
      }
      std::vector<long double> childValues(2U*count,0.0L);
      std::vector<long double> childErrors(2U*count,0.0L);
      unsigned childAccepted = 0;
      if(phase != 0) {
        if(phase > 0) {
          for(unsigned i = 0; i <= degree; ++i) {
            for(unsigned j = 0; j <= degree-i; ++j) {
              const unsigned k = LineSliceUtilities::MomentIndex(degree,i,j);
              childValues[k] = 1.0L/((i+1.0L)*(j+1.0L));
              childErrors[k] = 4.0L*std::numeric_limits<long double>::epsilon()*std::abs(childValues[k]);
            }
          }
        }
      }
      else {
        long double lower[2] = {};
        long double upper[2] = {};
        unsigned childOrder = Quad9SelectedOrder(childNormalized,lower,upper);
        long double derivativeBound = lower[2U-childOrder];
        const unsigned weakOrder = derivativeBound == 0.0L ? Quad9SelectedWeakGraphOrder(childNormalized) : 0U;
        if(derivativeBound == 0.0L && weakOrder != 0) childOrder = weakOrder;
        if(derivativeBound == 0.0L && weakOrder == 0 && splitDepth < maximumSplitDepth) {
          Quad9SelectedSplitCanonical(childNormalized,scale*childScale,degree,computeDelta,absoluteTolerance,0.25L*relativeTolerance,maximumDepth,maximumEvaluations,&boundary,useTaylor,surfaceGradient == nullptr ? nullptr : &childGradient,splitDepth+1U,maximumSplitDepth,childValues,childErrors,evaluations,childAccepted);
        }
        else {
          Quad9SelectedOneDirection(childNormalized,scale*childScale,degree,childOrder,derivativeBound,computeDelta,false,absoluteTolerance,0.25L*relativeTolerance,maximumDepth,maximumEvaluations,&boundary,childValues,childErrors,evaluations,childAccepted,useTaylor,surfaceGradient == nullptr ? nullptr : &childGradient);
        }
      }
      accepted += childAccepted;
      for(unsigned k = count; k < 2U*count; ++k) {
        childValues[k] /= childScale;
        childErrors[k] /= childScale;
      }
      std::vector<long double> mappedValues;
      std::vector<long double> mappedErrors;
      ReferenceAffineMoments2D(degree,childValues,childErrors,firstOrigin,firstScale,secondOrigin,secondScale,firstScale*secondScale,mappedValues,mappedErrors);
      for(unsigned k = 0; k < values.size(); ++k) {
        const long double adjusted = mappedValues[k]-correction[k];
        const long double next = values[k]+adjusted;
        correction[k] = (next-values[k])-adjusted;
        values[k] = next;
        errors[k] += mappedErrors[k];
      }
    }
  }
}

inline void Quad9SelectedReferenceOneDirection(const std::vector<std::vector<long double>>& coefficient, long double scale, unsigned degree, unsigned order, long double derivativeBound, bool computeDelta, bool fullBoundaryDelta, long double absoluteTolerance, long double relativeTolerance, unsigned maximumDepth, unsigned maximumEvaluations, const std::vector<unsigned>* zeroBoundary, std::vector<long double>& values, std::vector<long double>& errors, unsigned& evaluations, unsigned& accepted, bool useTaylor, const std::vector<std::vector<std::vector<long double>>>* surfaceGradient = nullptr) {
  const unsigned count = (degree+1U)*(degree+2U)/2U;
  const long double amplification = 4.0L*ReferencePower(3.0L,degree);
  if(absoluteTolerance/amplification > 64.0L*std::numeric_limits<long double>::epsilon()) {
    try {
      std::vector<long double> chartValues;
      std::vector<long double> chartErrors;
      if(derivativeBound == 0.0L && !fullBoundaryDelta) Quad9SelectedSplitCanonical(coefficient,scale,degree,computeDelta,absoluteTolerance/amplification,relativeTolerance/amplification,maximumDepth,maximumEvaluations,zeroBoundary,useTaylor,surfaceGradient,0,5,chartValues,chartErrors,evaluations,accepted);
      else Quad9SelectedOneDirection(coefficient,scale,degree,order,derivativeBound,computeDelta,fullBoundaryDelta,absoluteTolerance/amplification,relativeTolerance/amplification,maximumDepth,maximumEvaluations,zeroBoundary,chartValues,chartErrors,evaluations,accepted,useTaylor,surfaceGradient);
      ReferenceTensorMoments2D(degree,chartValues,chartErrors,values,errors);
      Quad9SelectedReferenceParity(coefficient,degree,values,errors,surfaceGradient);
      bool completed = true;
      for(unsigned k = 0; k < values.size(); ++k) {
        const long double tolerance = absoluteTolerance*(k < count ? 1.0L : scale)+relativeTolerance*std::abs(values[k]);
        completed = completed && std::isfinite(values[k]) && std::isfinite(errors[k]) && errors[k] >= 0.0L && errors[k] <= tolerance;
      }
      if(completed) return;
    }
    catch(const std::runtime_error&) {}
  }
  values.assign(2U*count,0.0L);
  errors.assign(2U*count,0.0L);
  std::vector<long double> correction(2U*count,0.0L);
  std::vector<long double> absoluteSum(2U*count,0.0L);
  accepted = 0;
  const long double quadrantTolerance = absoluteTolerance/16.0L;
  const long double quadrantRelativeTolerance = relativeTolerance/16.0L;
  for(unsigned firstSide = 0; firstSide < 2; ++firstSide) {
    const long double firstSign = firstSide == 0 ? -1.0L : 1.0L;
    for(unsigned secondSide = 0; secondSide < 2; ++secondSide) {
      const long double secondSign = secondSide == 0 ? -1.0L : 1.0L;
      std::vector<std::vector<long double>> quadrantCoefficient;
      ReferencePolynomialAffine2D(coefficient,0.5L,0.5L*firstSign,0.5L,0.5L*secondSign,quadrantCoefficient);
      std::vector<std::vector<long double>> normalized;
      long double quadrantScale;
      LineSliceUtilities::NormalizeCanonicalPolynomial(quadrantCoefficient,normalized,quadrantScale);
      std::vector<std::vector<std::vector<long double>>> quadrantGradient;
      if(surfaceGradient != nullptr) {
        quadrantGradient.resize(surfaceGradient->size());
        for(unsigned component = 0; component < surfaceGradient->size(); ++component) ReferencePolynomialAffine2D((*surfaceGradient)[component],0.5L,0.5L*firstSign,0.5L,0.5L*secondSign,quadrantGradient[component]);
      }
      std::vector<unsigned> boundary(4,1);
      for(unsigned i = 0; i < 3; ++i) {
        boundary[0] = boundary[0] != 0 && normalized[0][i] == 0.0L;
        boundary[1] = boundary[1] != 0 && (normalized[0][i]+normalized[1][i])+normalized[2][i] == 0.0L;
        boundary[2] = boundary[2] != 0 && normalized[i][0] == 0.0L;
        boundary[3] = boundary[3] != 0 && (normalized[i][0]+normalized[i][1])+normalized[i][2] == 0.0L;
      }
      if(zeroBoundary != nullptr) {
        boundary[1] = boundary[1] != 0 || (*zeroBoundary)[firstSide] != 0;
        boundary[3] = boundary[3] != 0 || (*zeroBoundary)[2U+secondSide] != 0;
      }
      std::vector<long double> quadrantValues(2U*count,0.0L);
      std::vector<long double> quadrantErrors(2U*count,0.0L);
      int phase = 0;
      LineSliceUtilities::CanonicalGraphOrder(normalized,0,phase);
      for(unsigned edge = 0; edge < 4; ++edge) if(computeDelta && boundary[edge] != 0) phase = 0;
      if(phase != 0) {
        if(phase > 0) {
          for(unsigned i = 0; i <= degree; ++i) {
            for(unsigned j = 0; j <= degree-i; ++j) {
              const unsigned k = LineSliceUtilities::MomentIndex(degree,i,j);
              quadrantValues[k] = 1.0L/((i+1.0L)*(j+1.0L));
              quadrantErrors[k] = 4.0L*std::numeric_limits<long double>::epsilon()*std::abs(quadrantValues[k]);
            }
          }
        }
      }
      else {
        long double lower[2];
        long double upper[2];
        const unsigned quadrantOrder = Quad9SelectedOrder(normalized,lower,upper);
        unsigned quadrantAccepted = 0;
        Quad9SelectedOneDirection(normalized,scale*quadrantScale,degree,quadrantOrder,lower[2U-quadrantOrder],computeDelta,false,quadrantTolerance,quadrantRelativeTolerance,maximumDepth,maximumEvaluations,&boundary,quadrantValues,quadrantErrors,evaluations,quadrantAccepted,useTaylor,surfaceGradient == nullptr ? nullptr : &quadrantGradient);
        accepted += quadrantAccepted;
      }
      if(computeDelta && fullBoundaryDelta) {
        for(unsigned axis = 0; axis < 2; ++axis) {
          const unsigned edge = 2U*axis+1U;
          if(boundary[edge] == 0) continue;
          std::vector<std::vector<long double>> edgeCoefficient(3,std::vector<long double>(3,0.0L));
          for(unsigned k = 0; k < 3; ++k) {
            const long double derivative = axis == 0 ? normalized[1][k]+2.0L*normalized[2][k] : normalized[k][1]+2.0L*normalized[k][2];
            if(axis == 0) {
              edgeCoefficient[0][k] = -derivative;
              edgeCoefficient[1][k] = derivative;
            }
            else {
              edgeCoefficient[k][0] = -derivative;
              edgeCoefficient[k][1] = derivative;
            }
          }
          std::vector<std::vector<long double>> normalizedEdge;
          long double edgeScale;
          LineSliceUtilities::NormalizeCanonicalPolynomial(edgeCoefficient,normalizedEdge,edgeScale);
          long double lower[2];
          long double upper[2];
          Quad9SelectedOrder(normalizedEdge,lower,upper);
          std::vector<unsigned> edgeBoundary(4,0);
          edgeBoundary[edge] = 1;
          std::vector<long double> edgeValues;
          std::vector<long double> edgeErrors;
          unsigned edgeAccepted = 0;
          const unsigned edgeOrder = 2U-axis;
          Quad9SelectedOneDirection(normalizedEdge,scale*quadrantScale*edgeScale,degree,edgeOrder,lower[axis],true,false,quadrantTolerance,quadrantRelativeTolerance,maximumDepth,maximumEvaluations,&edgeBoundary,edgeValues,edgeErrors,evaluations,edgeAccepted,useTaylor,surfaceGradient == nullptr ? nullptr : &quadrantGradient);
          accepted += edgeAccepted;
          for(unsigned k = count; k < 2U*count; ++k) {
            const long double contribution = edgeValues[k]/edgeScale;
            quadrantValues[k] += contribution;
            quadrantErrors[k] += edgeErrors[k]/edgeScale+4.0L*std::numeric_limits<long double>::epsilon()*std::abs(contribution);
          }
        }
      }
      for(unsigned i = 0; i <= degree; ++i) {
        for(unsigned j = 0; j <= degree-i; ++j) {
          const unsigned k = LineSliceUtilities::MomentIndex(degree,i,j);
          const long double sign = ReferencePower(firstSign,i)*ReferencePower(secondSign,j);
          for(unsigned pair = 0; pair < 2; ++pair) {
            const unsigned index = pair*count+k;
            const long double multiplier = pair == 0 ? sign : sign/quadrantScale;
            const long double term = multiplier*quadrantValues[index];
            const long double adjusted = term-correction[index];
            const long double next = values[index]+adjusted;
            correction[index] = (next-values[index])-adjusted;
            values[index] = next;
            errors[index] += std::abs(multiplier)*quadrantErrors[index];
            absoluteSum[index] += std::abs(term);
          }
        }
      }
    }
  }
  for(unsigned k = 0; k < values.size(); ++k) errors[k] += 16.0L*std::numeric_limits<long double>::epsilon()*absoluteSum[k];
  Quad9SelectedReferenceParity(coefficient,degree,values,errors,surfaceGradient);
  for(unsigned k = 0; k < values.size(); ++k) {
    const long double tolerance = absoluteTolerance*(k < count ? 1.0L : scale)+relativeTolerance*std::abs(values[k]);
    if(!std::isfinite(values[k]) || !std::isfinite(errors[k]) || errors[k] < 0.0L || errors[k] > tolerance) throw std::runtime_error("QUAD9 native quadrant error exceeds tolerance.");
  }
}

inline void Quad9SurfaceCoordinatePolynomials(const std::vector<std::vector<double>>& geometryCoordinates, std::vector<std::vector<std::vector<long double>>>& coordinatePolynomials) {
  ValidateSurfaceGeometryCoordinates(geometryCoordinates,2,9);
  coordinatePolynomials.assign(2,std::vector<std::vector<long double>>(3,std::vector<long double>(3,0.0L)));
  for(unsigned component = 0; component < 2; ++component) {
    std::vector<double> shifted(9);
    for(unsigned node = 0; node < 9; ++node) shifted[node] = geometryCoordinates[component][node]-geometryCoordinates[component][0];
    std::vector<std::vector<double>> coefficient;
    Quad9ReferenceCoefficientsFromNodalValues(shifted,coefficient);
    for(unsigned i = 0; i < 3; ++i) {
      for(unsigned j = 0; j < 3; ++j) coordinatePolynomials[component][i][j] = coefficient[i][j];
    }
  }
}

inline void IntegrateQuad9ReferenceSelected(const std::vector<std::vector<double>>& coefficient, unsigned degree, std::vector<std::vector<double>>& positive, std::vector<std::vector<double>>& delta, double absoluteTolerance = 1.0e-12, double relativeTolerance = 1.0e-12, bool computeDelta = true, bool fullBoundaryDelta = false, unsigned maximumDepth = 28, unsigned integrationOrder = 0, unsigned maximumEvaluations = 200000, std::vector<double>* estimatedError = nullptr, unsigned validationMode = 1, const std::vector<unsigned>* zeroBoundary = nullptr, bool useTaylor = false, bool referenceSurface = false, const std::vector<std::vector<std::vector<long double>>>* surfaceGradient = nullptr, const std::vector<std::vector<double>>* inverseJacobian = nullptr, const std::vector<std::vector<double>>* geometryCoordinates = nullptr) {
  if(geometryCoordinates != nullptr) {
    if(!referenceSurface || inverseJacobian != nullptr) throw std::invalid_argument("Surface geometry requires surface moments without an inverse Jacobian.");
    ValidateSurfaceGeometryCoordinates(*geometryCoordinates,2,9);
  }
  if(inverseJacobian != nullptr) {
    if(!referenceSurface) throw std::invalid_argument("A surface inverse Jacobian requires surface moments.");
    ValidateSurfaceInverseJacobian(*inverseJacobian,2);
  }
  if(&positive == &delta || coefficient.size() != 3 || degree > 64 || !(absoluteTolerance > 0.0) || !(relativeTolerance >= 0.0) || !std::isfinite(absoluteTolerance) || !std::isfinite(relativeTolerance) || maximumDepth > 60 || integrationOrder > 2 || validationMode > 2 || maximumEvaluations < 15) throw std::invalid_argument("Invalid QUAD9 selected-direction input.");
  for(unsigned i = 0; i < 3; ++i) if(coefficient[i].size() != 3) throw std::invalid_argument("The QUAD9 coefficient matrix must be 3 by 3.");
  if(zeroBoundary != nullptr && zeroBoundary->size() != 4) throw std::invalid_argument("QUAD9 requires four boundary flags.");
  std::vector<std::vector<long double>> chartCoefficient;
  ReferencePolynomialAffine2D(coefficient,-1.0L,2.0L,-1.0L,2.0L,chartCoefficient);
  std::vector<std::vector<long double>> normalized;
  long double scale;
  LineSliceUtilities::NormalizeCanonicalPolynomial(chartCoefficient,normalized,scale);
  std::vector<std::vector<std::vector<long double>>> referenceGradient;
  if(referenceSurface && surfaceGradient == nullptr) {
    SurfaceGradientFromPolynomial(coefficient,1.0L,1.0L,referenceGradient);
    surfaceGradient = &referenceGradient;
  }
  if(inverseJacobian != nullptr) {
    if(surfaceGradient != &referenceGradient) referenceGradient = *surfaceGradient;
    SurfaceGradientPhysical(referenceGradient,*inverseJacobian);
    surfaceGradient = &referenceGradient;
  }
  if(geometryCoordinates != nullptr) {
    if(surfaceGradient != &referenceGradient) referenceGradient = *surfaceGradient;
    std::vector<std::vector<std::vector<long double>>> coordinatePolynomials;
    Quad9SurfaceCoordinatePolynomials(*geometryCoordinates,coordinatePolynomials);
    SurfaceGradientCofactor2D(referenceGradient,coordinatePolynomials,1.0L,1.0L);
    surfaceGradient = &referenceGradient;
  }
  std::vector<std::vector<std::vector<long double>>> chartGradient;
  if(surfaceGradient != nullptr) {
    chartGradient.resize(surfaceGradient->size());
    for(unsigned component = 0; component < surfaceGradient->size(); ++component) ReferencePolynomialAffine2D((*surfaceGradient)[component],-1.0L,2.0L,-1.0L,2.0L,chartGradient[component]);
    surfaceGradient = &chartGradient;
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
        for(unsigned j = 0; j <= degree-i; ++j) {
          const unsigned k = LineSliceUtilities::MomentIndex(degree,i,j);
          values[k] = (norm == 0.0L ? 0.5L : 1.0L)*ReferenceIntervalMoment(i)*ReferenceIntervalMoment(j);
          errors[k] = 4.0L*std::numeric_limits<long double>::epsilon()*std::abs(values[k]);
        }
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
          Quad9SelectedReferenceOneDirection(normalized,scale,degree,order,lower[2U-order],computeDelta,fullBoundaryDelta,absoluteTolerance*factor,relativeTolerance*factor,maximumDepth,maximumEvaluations,zeroBoundary,trial[candidate],trialError[candidate],evaluations,trialAccepted[candidate],useTaylor,surfaceGradient);
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
      const long double momentErrorH = errors[k]+std::abs(static_cast<long double>(resultH[i][j])-values[k]);
      const long double momentErrorD = computeDelta ? errors[count+k]/scale+std::abs(static_cast<long double>(resultD[i][j])-values[count+k]/scale) : 0.0L;
      if(!std::isfinite(momentErrorH) || !std::isfinite(momentErrorD) || momentErrorH > absoluteTolerance+relativeTolerance*std::abs(static_cast<long double>(resultH[i][j])) || momentErrorD > absoluteTolerance+relativeTolerance*std::abs(static_cast<long double>(resultD[i][j]))) throw std::runtime_error("QUAD9 reference moment error exceeds tolerance.");
      errorH = std::max(errorH,momentErrorH);
      errorD = std::max(errorD,momentErrorD);
    }
  }
  positive.swap(resultH);
  delta.swap(resultD);
  if(estimatedError != nullptr) *estimatedError = {static_cast<double>(errorH),static_cast<double>(errorD),static_cast<double>(accepted),useTaylor ? -1.0 : static_cast<double>(evaluations),static_cast<double>(selected),static_cast<double>(compared),static_cast<double>(attempts)};
}

inline void IntegrateQuad9NodalSelected(const std::vector<double>& nodalValues, unsigned degree, std::vector<std::vector<double>>& positive, std::vector<std::vector<double>>& delta, double absoluteTolerance = 1.0e-12, double relativeTolerance = 1.0e-12, bool computeDelta = true, bool fullBoundaryDelta = false, unsigned validationMode = 1, unsigned integrationOrder = 0, std::vector<double>* estimatedError = nullptr, bool useTaylor = false, bool referenceSurface = false, const std::vector<std::vector<double>>* inverseJacobian = nullptr, const std::vector<std::vector<double>>* geometryCoordinates = nullptr) {
  std::vector<std::vector<double>> coefficient;
  Quad9ReferenceCoefficientsFromNodalValues(nodalValues,coefficient);
  std::vector<unsigned> boundary = {static_cast<unsigned>(nodalValues[0] == 0.0 && nodalValues[3] == 0.0 && nodalValues[7] == 0.0),static_cast<unsigned>(nodalValues[1] == 0.0 && nodalValues[2] == 0.0 && nodalValues[5] == 0.0),static_cast<unsigned>(nodalValues[0] == 0.0 && nodalValues[1] == 0.0 && nodalValues[4] == 0.0),static_cast<unsigned>(nodalValues[3] == 0.0 && nodalValues[2] == 0.0 && nodalValues[6] == 0.0)};
  IntegrateQuad9ReferenceSelected(coefficient,degree,positive,delta,absoluteTolerance,relativeTolerance,computeDelta,fullBoundaryDelta,useTaylor ? 24U : 28U,integrationOrder,200000,estimatedError,validationMode,&boundary,useTaylor,referenceSurface,nullptr,inverseJacobian,geometryCoordinates);
}

}

#endif

