#ifndef FEMUS_HEX27_DIRECTION_SELECTION_HPP
#define FEMUS_HEX27_DIRECTION_SELECTION_HPP

#include "Hex27LevelSetMoments.hpp"
#include "Quad9DirectionSelection.hpp"

// Unit-cube HEX27 moments: order 0 is automatic; orders 1..6 are xyz, xzy, yxz, yzx, zxy, zyx, from inner to outer.
namespace FELevelSetMoments {

inline void Hex27SelectedAxes(unsigned order, unsigned& inner, unsigned& middle, unsigned& outer) {
  if(order == 0 || order > 6) throw std::invalid_argument("Invalid explicit HEX27 integration order.");
  const unsigned axes[6][3] = {{0,1,2},{0,2,1},{1,0,2},{1,2,0},{2,0,1},{2,1,0}};
  inner = axes[order-1U][0];
  middle = axes[order-1U][1];
  outer = axes[order-1U][2];
}

inline void Hex27SelectedDerivativeBounds(const std::vector<long double>& coefficient, const std::vector<long double>& bernstein, unsigned axis, long double& minimum, long double& maximum, long double& energy) {
  if(coefficient.size() != 27 || bernstein.size() != 27 || axis > 2) throw std::invalid_argument("Invalid HEX27 derivative input.");
  long double norm = 0.0L;
  for(unsigned n = 0; n < 27; ++n) {
    if(!std::isfinite(coefficient[n]) || !std::isfinite(bernstein[n])) throw std::invalid_argument("A HEX27 derivative coefficient is not finite.");
    norm += std::abs(coefficient[n]);
  }
  unsigned first;
  unsigned second;
  Hex27OtherAxes(axis,first,second);
  minimum = std::numeric_limits<long double>::infinity();
  maximum = -minimum;
  energy = 0.0L;
  for(unsigned i = 0; i < 3; ++i) {
    for(unsigned j = 0; j < 3; ++j) {
      unsigned index[3] = {0,0,0};
      index[first] = i;
      index[second] = j;
      index[axis] = 0;
      long double previous = bernstein[Hex27Index(index[0],index[1],index[2])];
      for(unsigned side = 1; side < 3; ++side) {
        index[axis] = side;
        const long double next = bernstein[Hex27Index(index[0],index[1],index[2])];
        const long double derivative = 2.0L*(next-previous);
        minimum = std::min(minimum,derivative);
        maximum = std::max(maximum,derivative);
        energy += derivative*derivative;
        previous = next;
      }
    }
  }
  const long double padding = 2048.0L*std::numeric_limits<long double>::epsilon()*norm;
  minimum -= padding;
  maximum += padding;
  energy /= 18.0L;
}

inline void Hex27SelectedRank(const std::vector<long double>& coefficient, const std::vector<long double>& bernstein, long double lower[3], long double upper[3], long double energy[3], bool independent[3], unsigned axes[3]) {
  independent[0] = independent[1] = independent[2] = true;
  for(unsigned i = 0; i < 3; ++i) {
    for(unsigned j = 0; j < 3; ++j) {
      for(unsigned k = 0; k < 3; ++k) {
        if(coefficient[Hex27Index(i,j,k)] != 0.0L) {
          if(i != 0) independent[0] = false;
          if(j != 0) independent[1] = false;
          if(k != 0) independent[2] = false;
        }
      }
    }
  }
  long double transverse[3];
  for(unsigned axis = 0; axis < 3; ++axis) {
    long double minimum;
    long double maximum;
    Hex27SelectedDerivativeBounds(coefficient,bernstein,axis,minimum,maximum,energy[axis]);
    upper[axis] = std::max(std::abs(minimum),std::abs(maximum));
    lower[axis] = minimum > 0.0L ? minimum : (maximum < 0.0L ? -maximum : 0.0L);
    if(lower[axis] <= 32.0L*std::sqrt(std::numeric_limits<long double>::epsilon())*upper[axis]) lower[axis] = 0.0L;
    axes[axis] = axis;
  }
  for(unsigned axis = 0; axis < 3; ++axis) transverse[axis] = std::max(lower[(axis+1U)%3U],lower[(axis+2U)%3U]);
  for(unsigned i = 0; i < 3; ++i) {
    for(unsigned j = i+1U; j < 3; ++j) {
      const unsigned first = axes[i];
      const unsigned second = axes[j];
      const bool better = independent[first] != independent[second] ? independent[second] : (transverse[first] != transverse[second] ? transverse[second] > transverse[first] : energy[second] < energy[first]);
      if(better) std::swap(axes[i],axes[j]);
    }
  }
}

inline void Hex27SelectedSlice(const std::vector<long double>& coefficient, unsigned axis, long double point, std::vector<std::vector<long double>>& slice) {
  std::vector<std::vector<double>> original;
  Hex27Slice(coefficient,axis,point,original);
  slice.assign(3,std::vector<long double>(3,0.0L));
  for(unsigned i = 0; i < 3; ++i) {
    for(unsigned j = 0; j < 3; ++j) slice[i][j] = original[i][j];
  }
}

inline void Hex27SelectedSample(const std::vector<long double>& coefficient, unsigned degree, unsigned axis, unsigned innerAxis, long double point, long double jacobian, long double absoluteTolerance, long double deltaTolerance, long double relativeTolerance, bool computeDelta, std::vector<long double>& value, std::vector<long double>& error, bool& cut, unsigned& evaluations, unsigned maximumEvaluations, bool useTaylor, const std::vector<unsigned>* zeroFaces, std::vector<unsigned>& innerCounts) {
  if(evaluations >= maximumEvaluations) throw std::runtime_error("HEX27 slice evaluation budget exceeded.");
  ++evaluations;
  unsigned first;
  unsigned second;
  Hex27OtherAxes(axis,first,second);
  std::vector<std::vector<long double>> slice;
  Hex27SelectedSlice(coefficient,axis,point,slice);
  std::vector<unsigned> boundary(4,0);
  bool anyBoundary = false;
  for(unsigned edge = 0; edge < 4; ++edge) {
    const unsigned localAxis = edge/2U;
    const unsigned globalAxis = localAxis == 0 ? first : second;
    bool zero = true;
    for(unsigned i = 0; i < 3; ++i) {
      const long double edgeValue = localAxis == 0 ? (edge%2U == 0 ? slice[0][i] : (slice[0][i]+slice[1][i])+slice[2][i]) : (edge%2U == 0 ? slice[i][0] : (slice[i][0]+slice[i][1])+slice[i][2]);
      zero = zero && edgeValue == 0.0L;
    }
    boundary[edge] = zero || (zeroFaces != nullptr && (*zeroFaces)[2U*globalAxis+edge%2U] != 0);
    anyBoundary = anyBoundary || boundary[edge] != 0;
  }
  const unsigned count = (degree+1U)*(degree+2U)*(degree+3U)/6U;
  const unsigned sliceCount = (degree+1U)*(degree+2U)/2U;
  std::vector<long double> moments(2U*sliceCount,0.0L);
  std::vector<long double> sliceError(2U*sliceCount,0.0L);
  long double norm = 0.0L;
  for(unsigned i = 0; i < 3; ++i) {
    for(unsigned j = 0; j < 3; ++j) norm = std::max(norm,std::abs(slice[i][j]));
  }
  if(norm == 0.0L && computeDelta) throw std::domain_error("A HEX27 slice is identically zero.");
  int phase = 0;
  LineSliceUtilities::CanonicalGraphOrder(slice,0,phase);
  if(computeDelta && anyBoundary) phase = 0;
  if(phase != 0 || norm == 0.0L) {
    if(phase > 0 || norm == 0.0L) {
      for(unsigned i = 0; i <= degree; ++i) {
        for(unsigned j = 0; j <= degree-i; ++j) moments[MomentIndex(degree,i,j)] = (norm == 0.0L ? 0.5L : 1.0L)/((i+1.0L)*(j+1.0L));
      }
    }
  }
  else {
    for(unsigned i = 0; i < 3; ++i) {
      for(unsigned j = 0; j < 3; ++j) slice[i][j] /= norm;
    }
    long double lower[2];
    long double upper[2];
    const unsigned automatic = Quad9SelectedOrder(slice,lower,upper);
    const unsigned preferred = innerAxis == 3 ? automatic : (innerAxis == first ? 2U : 1U);
    const long double innerTolerance = std::max(static_cast<long double>(std::numeric_limits<double>::min()),0.02L*std::min(absoluteTolerance,computeDelta ? deltaTolerance : absoluteTolerance));
    unsigned lineEvaluations = 0;
    unsigned accepted = 0;
    bool completed = false;
    std::string failure = "No HEX27 slice direction converged.";
    for(unsigned attempt = 0; attempt < (innerAxis == 3 ? 2U : 1U) && !completed; ++attempt) {
      const unsigned order = attempt == 0 ? preferred : 3U-preferred;
      try {
        Quad9SelectedOneDirection(slice,norm,degree,order,lower[2U-order],computeDelta,false,innerTolerance,0.05L*relativeTolerance,28,200000,anyBoundary ? &boundary : nullptr,moments,sliceError,lineEvaluations,accepted,useTaylor);
        ++innerCounts[order == 2 ? first : second];
        completed = true;
      }
      catch(const std::runtime_error& exception) { failure = exception.what(); }
    }
    if(!completed) throw std::runtime_error(failure);
    for(unsigned k = sliceCount; k < 2U*sliceCount; ++k) {
      moments[k] /= norm;
      sliceError[k] /= norm;
    }
  }
  for(unsigned k = 0; k < moments.size(); ++k) {
    const double rounded = static_cast<double>(moments[k]);
    if(!std::isfinite(rounded)) throw std::overflow_error("A HEX27 slice moment exceeds double range.");
    sliceError[k] += std::numeric_limits<double>::epsilon()*std::abs(rounded);
    moments[k] = rounded;
  }
  value.assign(2U*count,0.0L);
  error.assign(2U*count,0.0L);
  for(unsigned i = 0; i <= degree; ++i) {
    for(unsigned j = 0; j <= degree-i; ++j) {
      const unsigned source = MomentIndex(degree,i,j);
      long double power = 1.0L;
      for(unsigned k = 0; k <= degree-i-j; ++k) {
        unsigned index[3] = {0,0,0};
        index[first] = i;
        index[second] = j;
        index[axis] = k;
        const unsigned target = HexMomentIndex(index[0],index[1],index[2]);
        const long double factor = jacobian*power;
        value[target] = factor*moments[source];
        error[target] = factor*sliceError[source]+8.0L*std::numeric_limits<long double>::epsilon()*std::abs(value[target]);
        if(computeDelta) {
          value[count+target] = factor*moments[sliceCount+source];
          error[count+target] = factor*sliceError[sliceCount+source]+8.0L*std::numeric_limits<long double>::epsilon()*std::abs(value[count+target]);
        }
        power *= point;
      }
    }
  }
  cut = (moments[0] > 0.0L && moments[0] < 1.0L) || (computeDelta && moments[sliceCount] != 0.0L);
}

inline void Hex27SelectedInterval(const std::vector<long double>& coefficient, const std::vector<long double>& bernstein, unsigned degree, unsigned axis, unsigned innerAxis, long double left, long double right, long double thetaLeft, long double thetaRight, long double absoluteTolerance, long double deltaTolerance, long double relativeTolerance, bool computeDelta, unsigned depth, unsigned maximumDepth, unsigned& evaluations, unsigned maximumEvaluations, unsigned& accepted, std::vector<long double>& total, std::vector<long double>& totalError, bool regularGraph, bool useTaylor, const std::vector<unsigned>* zeroFaces, std::vector<unsigned>& innerCounts) {
  const long double nodes[8] = {0.991455371120812639206854697526329L,0.949107912342758524526189684047851L,0.864864423359769072789712788640926L,0.741531185599394439863864773280788L,0.586087235467691130294144838258730L,0.405845151377397166906606412076961L,0.207784955007898467600689403773245L,0.0L};
  const long double kronrod[8] = {0.022935322010529224963732008058970L,0.063092092629978553290700663189204L,0.104790010322250183839876322541518L,0.140653259715525918745189590510238L,0.169004726639267902826583426598550L,0.190350578064785409913256402421014L,0.204432940075298892414161999234649L,0.209482141084727828012999174891714L};
  const long double gauss[4] = {0.129484966168869693270611432679082L,0.279705391489276667901467771423780L,0.381830050505118944950369775488975L,0.417959183673469387755102040816327L};
  long double width = right-left;
  long double zLeft = left + width*(regularGraph ? thetaLeft : Power(std::sin(thetaLeft),2));
  long double zRight = left + width*(regularGraph ? thetaRight : Power(std::sin(thetaRight),2));
  if(!(zRight > zLeft)) throw std::runtime_error("HEX27 slice interval reached floating point resolution.");
  unsigned count = (degree+1U)*(degree+2U)*(degree+3U)/6U;
  int phase = regularGraph ? 0 : Hex27SlabPhase(bernstein,axis,zLeft,zRight,depth,computeDelta);
  if(phase != 0) {
    if(phase > 0) {
      for(unsigned i = 0; i <= degree; ++i) {
        for(unsigned j = 0; j <= degree-i; ++j) {
          for(unsigned k = 0; k <= degree-i-j; ++k) {
            unsigned powers[3] = {i,j,k};
            long double value = MonomialIntegral(zLeft,zRight,powers[axis]);
            for(unsigned d = 0; d < 3; ++d) if(d != axis) value /= powers[d]+1U;
            total[HexMomentIndex(i,j,k)] += value;
          }
        }
      }
    }
    ++accepted;
    return;
  }
  long double middle = 0.5L*(thetaLeft+thetaRight);
  long double halfWidth = 0.5L*(thetaRight-thetaLeft);
  std::vector<long double> high(2U*count,0.0L), low(2U*count,0.0L), innerError(2U*count,0.0L), absolute(2U*count,0.0L);
  std::vector<std::vector<long double>> samples;
  std::vector<long double> sampleWeights;
  bool anyCut = false;
  bool failed = false;
  long double split = middle;
  for(unsigned node = 0; node < 8 && !failed; ++node) {
    unsigned sides = node == 7 ? 1 : 2;
    for(unsigned side = 0; side < sides; ++side) {
      long double theta = middle + (side == 0 ? -1.0L : 1.0L)*halfWidth*nodes[node];
      long double sine = std::sin(theta);
      long double point = left + width*(regularGraph ? theta : sine*sine);
      long double jacobian = regularGraph ? width : 2.0L*width*sine*std::cos(theta);
      std::vector<long double> value, error;
      bool sampleCut = false;
      try {
        Hex27SelectedSample(coefficient,degree,axis,innerAxis,point,jacobian,absoluteTolerance,deltaTolerance,relativeTolerance,computeDelta,value,error,sampleCut,evaluations,maximumEvaluations,useTaylor,zeroFaces,innerCounts);
      }
      catch(const std::domain_error&) { failed = true; split = theta; break; }
      anyCut = anyCut || sampleCut;
      samples.push_back(value);
      sampleWeights.push_back(kronrod[node]);
      for(unsigned k = 0; k < high.size(); ++k) {
        high[k] += kronrod[node]*value[k];
        absolute[k] += kronrod[node]*std::abs(value[k]);
        innerError[k] += kronrod[node]*error[k];
        if((node & 1U) != 0) low[k] += gauss[node/2U]*value[k];
      }
    }
  }
  bool converged = !failed && anyCut;
  std::vector<long double> errors(2U*count,0.0L);
  long double intervalFraction = std::max(1.0e-2L,width*(thetaRight-thetaLeft)/(regularGraph ? 1.0L : 0.5L*std::acos(-1.0L)));
  if(failed && zRight-zLeft <= 0.25L*absoluteTolerance*intervalFraction) {
    long double bound = computeDelta ? Hex27SlabDeltaBound(bernstein,axis,zLeft,zRight,0.25L*deltaTolerance*intervalFraction) : 0.0L;
    if(bound <= 0.25L*deltaTolerance*intervalFraction) {
      for(unsigned k = 0; k < count; ++k) { totalError[k] += zRight-zLeft; totalError[count+k] += bound; }
      ++accepted;
      return;
    }
  }
  if(!failed) {
    long double fraction = intervalFraction;
    if(!anyCut && zRight-zLeft <= absoluteTolerance*fraction) {
      long double bound = computeDelta ? Hex27SlabDeltaBound(bernstein,axis,zLeft,zRight,deltaTolerance*fraction) : 0.0L;
      if(bound <= deltaTolerance*fraction) {
        converged = true;
        for(unsigned k = 0; k < count; ++k) { innerError[k] += (zRight-zLeft)/halfWidth; innerError[count+k] += bound/halfWidth; }
      }
    }
    for(unsigned k = 0; k < high.size(); ++k) {
      long double variation = 0.0L;
      for(unsigned node = 0; node < samples.size(); ++node) variation += sampleWeights[node]*std::abs(samples[node][k]-0.5L*high[k]);
      long double error = std::abs(high[k]-low[k]);
      if(variation > 0.0L && error > 0.0L) error = variation*std::min(1.0L,std::pow(200.0L*error/variation,1.5L));
      error = std::max(error,50.0L*std::numeric_limits<long double>::epsilon()*absolute[k]);
      errors[k] = halfWidth*(error+innerError[k]);
      high[k] *= halfWidth;
      long double tolerance = (k < count ? absoluteTolerance : deltaTolerance)*fraction + relativeTolerance*std::abs(high[k]);
      if(errors[k] > tolerance) converged = false;
    }
  }
  if(converged) {
    for(unsigned k = 0; k < high.size(); ++k) { total[k] += high[k]; totalError[k] += errors[k]; }
    ++accepted;
    return;
  }
  if(depth >= maximumDepth || !(split > thetaLeft && split < thetaRight)) throw std::runtime_error("HEX27 integration did not resolve the requested tolerance or slice topology.");
  Hex27SelectedInterval(coefficient,bernstein,degree,axis,innerAxis,left,right,thetaLeft,split,absoluteTolerance,deltaTolerance,relativeTolerance,computeDelta,depth+1U,maximumDepth,evaluations,maximumEvaluations,accepted,total,totalError,regularGraph,useTaylor,zeroFaces,innerCounts);
  Hex27SelectedInterval(coefficient,bernstein,degree,axis,innerAxis,left,right,split,thetaRight,absoluteTolerance,deltaTolerance,relativeTolerance,computeDelta,depth+1U,maximumDepth,evaluations,maximumEvaluations,accepted,total,totalError,regularGraph,useTaylor,zeroFaces,innerCounts);
}

inline void Hex27SelectedOneAxis(const std::vector<long double>& normalized, const std::vector<long double>& bernstein, unsigned degree, unsigned axis, unsigned innerAxis, bool independent, unsigned graphAxes, long double scale, long double absoluteTolerance, long double relativeTolerance, bool computeDelta, unsigned maximumDepth, unsigned maximumEvaluations, unsigned& evaluations, unsigned& accepted, std::vector<long double>& values, std::vector<long double>& errors, bool useTaylor, const std::vector<unsigned>* zeroFaces, std::vector<unsigned>& innerCounts) {
  const unsigned count = (degree+1U)*(degree+2U)*(degree+3U)/6U;
  const bool regularGraph = (graphAxes & ~(1U<<axis)) != 0;
  if(computeDelta && !regularGraph && Hex27FlatSlice(normalized,axis)) throw std::runtime_error("The requested HEX27 height contains an entire zero slice.");
  std::vector<long double> events;
  if(!independent) events = regularGraph ? std::vector<long double>{0.0L,1.0L} : Hex27Events(normalized,axis);
  long double positiveWorkingTolerance = absoluteTolerance;
  long double deltaWorkingTolerance = std::max(absoluteTolerance*scale,relativeTolerance);
  for(unsigned refinement = 0; refinement < 4; ++refinement) {
    values.assign(2U*count,0.0L);
    errors.assign(2U*count,0.0L);
    if(independent) {
      bool cut = false;
      Hex27SelectedSample(normalized,degree,axis,innerAxis,0.5L,1.0L,positiveWorkingTolerance,deltaWorkingTolerance,relativeTolerance,computeDelta,values,errors,cut,evaluations,maximumEvaluations,useTaylor,zeroFaces,innerCounts);
      for(unsigned i = 0; i <= degree; ++i) {
        for(unsigned j = 0; j <= degree-i; ++j) {
          for(unsigned k = 0; k <= degree-i-j; ++k) {
            const unsigned powers[3] = {i,j,k};
            const long double factor = 1.0L/(Power(0.5L,powers[axis])*(powers[axis]+1U));
            const unsigned index = HexMomentIndex(i,j,k);
            values[index] *= factor;
            values[count+index] *= factor;
            errors[index] *= factor;
            errors[count+index] *= factor;
          }
        }
      }
      ++accepted;
    }
    else {
      for(unsigned event = 1; event < events.size(); ++event) {
        if(events[event] > events[event-1U]) Hex27SelectedInterval(normalized,bernstein,degree,axis,innerAxis,events[event-1U],events[event],0.0L,regularGraph ? 1.0L : 0.5L*std::acos(-1.0L),positiveWorkingTolerance,deltaWorkingTolerance,relativeTolerance,computeDelta,0,maximumDepth,evaluations,maximumEvaluations,accepted,values,errors,regularGraph,useTaylor,zeroFaces,innerCounts);
      }
    }
    bool positiveConverged = true;
    bool deltaConverged = true;
    long double nextPositiveTolerance = 0.125L*positiveWorkingTolerance;
    long double nextDeltaTolerance = 0.125L*deltaWorkingTolerance;
    for(unsigned k = 0; k < values.size(); ++k) {
      if(!std::isfinite(values[k]) || !std::isfinite(errors[k]) || errors[k] < 0.0L) throw std::runtime_error("A HEX27 directional result is not finite.");
      const long double tolerance = absoluteTolerance*(k < count ? 1.0L : scale)+relativeTolerance*std::abs(values[k]);
      if(k < count) {
        positiveConverged = positiveConverged && errors[k] <= tolerance;
        nextPositiveTolerance = std::min(nextPositiveTolerance,0.25L*tolerance);
      }
      else if(computeDelta) {
        deltaConverged = deltaConverged && errors[k] <= tolerance;
        nextDeltaTolerance = std::min(nextDeltaTolerance,0.25L*tolerance);
      }
    }
    if(positiveConverged && deltaConverged) return;
    if(!positiveConverged) positiveWorkingTolerance = nextPositiveTolerance;
    if(!deltaConverged) deltaWorkingTolerance = nextDeltaTolerance;
  }
  throw std::runtime_error("HEX27 accumulated directional error exceeds tolerance.");
}

inline bool Hex27SelectedAgreement(const std::vector<long double>& first, const std::vector<long double>& firstError, const std::vector<long double>& second, const std::vector<long double>& secondError, unsigned count, long double scale, long double absoluteTolerance, long double relativeTolerance, bool computeDelta, std::vector<long double>& errors) {
  errors = firstError;
  bool agree = true;
  for(unsigned k = 0; k < (computeDelta ? 2U*count : count); ++k) {
    errors[k] = std::abs(first[k]-second[k])+firstError[k]+secondError[k];
    const long double tolerance = absoluteTolerance*(k < count ? 1.0L : scale)+relativeTolerance*std::abs(first[k]);
    agree = agree && std::isfinite(errors[k]) && errors[k] <= tolerance;
  }
  return agree;
}

inline void IntegrateHex27UnitSelected(const std::vector<std::vector<std::vector<double>>>& coefficient, unsigned degree, std::vector<std::vector<std::vector<double>>>& positive, std::vector<std::vector<std::vector<double>>>& delta, double absoluteTolerance = 1.0e-10, double relativeTolerance = 1.0e-10, bool computeDelta = true, unsigned maximumDepth = 24, unsigned integrationOrder = 0, unsigned maximumEvaluations = 200000, std::vector<double>* estimatedError = nullptr, unsigned validationMode = 1, bool useTaylor = true, const std::vector<unsigned>* zeroFaces = nullptr, std::vector<double>* directionInfo = nullptr) {
  if(&positive == &delta || degree > 64 || !(absoluteTolerance > 0.0) || !(relativeTolerance >= 0.0) || !std::isfinite(absoluteTolerance) || !std::isfinite(relativeTolerance) || maximumDepth > 48 || integrationOrder > 6 || maximumEvaluations == 0 || validationMode > 2 || (directionInfo != nullptr && directionInfo == estimatedError)) throw std::invalid_argument("Invalid selected-direction HEX27 integration input.");
  if(zeroFaces != nullptr) {
    if(zeroFaces->size() != 6) throw std::invalid_argument("HEX27 requires six zero-face flags.");
    for(unsigned face = 0; face < 6; ++face) {
      if((*zeroFaces)[face] > 1) throw std::invalid_argument("Invalid HEX27 zero-face flag.");
    }
  }
  std::vector<long double> normalized;
  long double scale;
  Hex27Normalize(coefficient,normalized,scale);
  const unsigned count = (degree+1U)*(degree+2U)*(degree+3U)/6U;
  std::vector<long double> values(2U*count,0.0L);
  std::vector<long double> errors(2U*count,0.0L);
  std::vector<long double> bernstein = Hex27Bernstein(normalized);
  bool planar = true;
  bool identicallyZero = true;
  for(unsigned i = 0; i < 3; ++i) {
    for(unsigned j = 0; j < 3; ++j) {
      for(unsigned k = 0; k < 3; ++k) {
        const long double value = normalized[Hex27Index(i,j,k)];
        if(value != 0.0L) {
          identicallyZero = false;
          if(i+j+k > 1) planar = false;
        }
      }
    }
  }
  if(identicallyZero && computeDelta) throw std::domain_error("The Dirac delta of an identically zero HEX27 field is undefined.");
  unsigned phaseBudget = 0;
  int phase = Hex27BernsteinPhase(bernstein,0,1,0,0,phaseBudget,computeDelta);
  if(computeDelta && zeroFaces != nullptr) {
    for(unsigned face = 0; face < 6; ++face) {
      if((*zeroFaces)[face] != 0) phase = 0;
    }
  }
  unsigned evaluations = 0;
  unsigned accepted = 0;
  unsigned attempts = 0;
  unsigned successes = 0;
  unsigned selectedAxis = 0;
  unsigned compared = 0;
  std::vector<unsigned> innerCounts(3,0);
  if(phase != 0) {
    if(phase > 0) {
      for(unsigned i = 0; i <= degree; ++i) {
        for(unsigned j = 0; j <= degree-i; ++j) {
          for(unsigned k = 0; k <= degree-i-j; ++k) values[HexMomentIndex(i,j,k)] = 1.0L/((i+1.0L)*(j+1.0L)*(k+1.0L));
        }
      }
    }
  }
  else if(planar) Hex27PlanarMoments(normalized,degree,computeDelta,absoluteTolerance,relativeTolerance,values);
  else {
    long double lower[3];
    long double upper[3];
    long double energy[3];
    bool independent[3];
    unsigned axes[3];
    Hex27SelectedRank(normalized,bernstein,lower,upper,energy,independent,axes);
    unsigned requestedInner = 3;
    unsigned requestedOuter = 3;
    if(integrationOrder != 0) {
      unsigned middle;
      Hex27SelectedAxes(integrationOrder,requestedInner,middle,requestedOuter);
      for(unsigned n = 0; n < 3; ++n) {
        if(axes[n] == requestedOuter) std::swap(axes[0],axes[n]);
      }
    }
    const long double transverse = requestedInner != 3 ? lower[requestedInner] : std::max(lower[(axes[0]+1U)%3U],lower[(axes[0]+2U)%3U]);
    const bool validate = validationMode == 2 || (validationMode == 1 && computeDelta && transverse == 0.0L);
    const unsigned graphAxes = Hex27GraphAxes(bernstein);
    const unsigned number = validate || integrationOrder == 0 ? 3U : 1U;
    bool completed = false;
    std::string failure = "No HEX27 direction converged.";
    for(unsigned retry = 0; retry < (validate ? 2U : 1U) && !completed; ++retry) {
      const long double factor = retry == 0 ? 1.0L : 0.125L;
      std::vector<std::vector<long double>> previousValues;
      std::vector<std::vector<long double>> previousErrors;
      std::vector<unsigned> previousAxes;
      for(unsigned candidate = 0; candidate < number && !completed; ++candidate) {
        const unsigned axis = axes[candidate];
        const unsigned inner = integrationOrder != 0 && candidate == 0 ? requestedInner : 3U;
        const long double candidateLower = inner == 3 ? std::max(lower[(axis+1U)%3U],lower[(axis+2U)%3U]) : lower[inner];
        if(!validate && validationMode == 1 && computeDelta && candidateLower == 0.0L) continue;
        std::vector<long double> trial;
        std::vector<long double> trialError;
        ++attempts;
        const std::vector<unsigned> countsBefore = innerCounts;
        try {
          Hex27SelectedOneAxis(normalized,bernstein,degree,axis,inner,independent[axis],graphAxes,scale,absoluteTolerance*factor,relativeTolerance*factor,computeDelta,maximumDepth,maximumEvaluations,evaluations,accepted,trial,trialError,useTaylor,zeroFaces,innerCounts);
          ++successes;
          if(!validate) {
            values.swap(trial);
            errors.swap(trialError);
            selectedAxis = axis+1U;
            completed = true;
          }
          else {
            if(independent[axis]) {
              unsigned first;
              unsigned second;
              Hex27OtherAxes(axis,first,second);
              const unsigned otherInner = innerCounts[first] > countsBefore[first] ? second : first;
              std::vector<long double> other;
              std::vector<long double> otherError;
              ++attempts;
              try {
                Hex27SelectedOneAxis(normalized,bernstein,degree,axis,otherInner,true,graphAxes,scale,absoluteTolerance*factor,relativeTolerance*factor,computeDelta,maximumDepth,maximumEvaluations,evaluations,accepted,other,otherError,useTaylor,zeroFaces,innerCounts);
                ++successes;
                std::vector<long double> checkedError;
                if((integrationOrder == 0 || axis == requestedOuter) && Hex27SelectedAgreement(trial,trialError,other,otherError,count,scale,absoluteTolerance,relativeTolerance,computeDelta,checkedError)) {
                  values = trial;
                  errors.swap(checkedError);
                  selectedAxis = axis+1U;
                  compared = 1;
                  completed = true;
                }
              }
              catch(const std::runtime_error& exception) { failure = exception.what(); }
              catch(const std::domain_error& exception) { failure = exception.what(); }
            }
            for(unsigned previous = 0; previous < previousValues.size() && !completed; ++previous) {
              if(integrationOrder != 0 && previousAxes[previous] != requestedOuter) continue;
              std::vector<long double> checkedError;
              if(Hex27SelectedAgreement(previousValues[previous],previousErrors[previous],trial,trialError,count,scale,absoluteTolerance,relativeTolerance,computeDelta,checkedError)) {
                values = previousValues[previous];
                errors.swap(checkedError);
                selectedAxis = previousAxes[previous]+1U;
                compared = 1;
                completed = true;
              }
            }
            previousValues.push_back(trial);
            previousErrors.push_back(trialError);
            previousAxes.push_back(axis);
            if(!completed) failure = "Independent HEX27 slicing directions did not agree within tolerance.";
          }
        }
        catch(const std::runtime_error& exception) { failure = exception.what(); }
        catch(const std::domain_error& exception) { failure = exception.what(); }
        if(evaluations >= maximumEvaluations && !completed) throw std::runtime_error("HEX27 total slice evaluation budget exhausted.");
        if(integrationOrder != 0 && candidate == 0 && previousValues.empty() && !completed) break;
      }
    }
    if(!completed) throw std::runtime_error("HEX27 selected-direction integration failed: "+failure);
  }
  std::vector<std::vector<std::vector<double>>> resultH(degree+1U);
  std::vector<std::vector<std::vector<double>>> resultD(degree+1U);
  long double maxH = 0.0L;
  long double maxD = 0.0L;
  for(unsigned i = 0; i <= degree; ++i) {
    resultH[i].resize(degree-i+1U);
    resultD[i].resize(degree-i+1U);
    for(unsigned j = 0; j <= degree-i; ++j) {
      resultH[i][j].resize(degree-i-j+1U);
      resultD[i][j].resize(degree-i-j+1U);
      for(unsigned k = 0; k <= degree-i-j; ++k) {
        const unsigned index = HexMomentIndex(i,j,k);
        const long double h = values[index];
        const long double d = computeDelta ? values[count+index]/scale : 0.0L;
        resultH[i][j][k] = static_cast<double>(h);
        resultD[i][j][k] = static_cast<double>(d);
        const long double hError = errors[index]+std::abs(h-static_cast<long double>(resultH[i][j][k]));
        const long double dError = computeDelta ? errors[count+index]/scale+std::abs(d-static_cast<long double>(resultD[i][j][k])) : 0.0L;
        if(!std::isfinite(resultH[i][j][k]) || !std::isfinite(resultD[i][j][k])) throw std::overflow_error("A HEX27 selected-direction moment exceeds double range.");
        if(hError > absoluteTolerance+relativeTolerance*std::abs(h) || dError > absoluteTolerance+relativeTolerance*std::abs(d)) throw std::runtime_error("A HEX27 output moment exceeds its estimated error budget.");
        maxH = std::max(maxH,hError);
        maxD = std::max(maxD,dError);
      }
    }
  }
  positive.swap(resultH);
  delta.swap(resultD);
  if(estimatedError != nullptr) *estimatedError = {static_cast<double>(maxH),static_cast<double>(maxD),static_cast<double>(accepted),static_cast<double>(evaluations),static_cast<double>(selectedAxis)};
  if(directionInfo != nullptr) *directionInfo = {static_cast<double>(selectedAxis),static_cast<double>(compared),static_cast<double>(attempts),static_cast<double>(successes),static_cast<double>(innerCounts[0]),static_cast<double>(innerCounts[1]),static_cast<double>(innerCounts[2])};
}

inline void IntegrateHex27NodalSelected(const std::vector<double>& nodalValues, unsigned degree, std::vector<std::vector<std::vector<double>>>& positive, std::vector<std::vector<std::vector<double>>>& delta, double absoluteTolerance = 1.0e-10, double relativeTolerance = 1.0e-10, bool computeDelta = true, unsigned validationMode = 1, unsigned integrationOrder = 0, std::vector<double>* estimatedError = nullptr, bool useTaylor = true, std::vector<double>* directionInfo = nullptr) {
  std::vector<std::vector<std::vector<double>>> coefficient;
  Hex27UnitCoefficientsFromNodalValues(nodalValues,coefficient);
  std::vector<unsigned> zeroFaces(6,0);
  for(unsigned axis = 0; axis < 3; ++axis) {
    for(unsigned side = 0; side < 2; ++side) {
      std::vector<double> face;
      Hex27FaceNodalValues(nodalValues,axis,side,face);
      bool zero = true;
      for(unsigned node = 0; node < 9; ++node) zero = zero && face[node] == 0.0;
      zeroFaces[2U*axis+side] = zero;
    }
  }
  IntegrateHex27UnitSelected(coefficient,degree,positive,delta,absoluteTolerance,relativeTolerance,computeDelta,24,integrationOrder,200000,estimatedError,validationMode,useTaylor,&zeroFaces,directionInfo);
}

}

#endif
