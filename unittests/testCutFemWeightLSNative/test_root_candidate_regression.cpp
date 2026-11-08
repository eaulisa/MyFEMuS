#include "AdaptiveLineSliceMoments.hpp"
#include <cmath>
#include <iomanip>
#include <iostream>
#include <vector>

unsigned candidateChecks = 0;
unsigned candidateFailures = 0;

void CheckCandidate(bool success, const char* description) {
  ++candidateChecks;
  if(success) return;
  ++candidateFailures;
  std::cerr << description << '\n';
}

void CheckDuplicateCandidate() {
  const std::vector<long double> polynomial = {3.66848065654445018278e-17L,-1.16299486905685379029e-15L,1.12631006249140947052e-15L,-2.56245806002129715809e-31L,7.42051998969827515356e-32L};
  std::vector<long double> candidates;
  LineSliceUtilities::PolynomialRoots(polynomial,0.0L,1.0L,candidates);
  CheckCandidate(candidates.size() == 3,"The duplicate-candidate fixture changed.");
  std::vector<std::pair<long double,long double>> events;
  std::vector<unsigned> sources;
  std::vector<long double> uncertainties;
  try {
    FELevelSetMoments::LineSlicePolynomialEvents(polynomial,1U,events,sources,uncertainties);
    CheckCandidate(events.size() == 2,"Duplicate candidates were not consolidated.");
    if(events.size() != 2) return;
    CheckCandidate(events[0] < events[1] && sources[0] == 1U && sources[1] == 1U,"Root provenance or ordering changed.");
    CheckCandidate(events[1].first == 1.0L && events[1].second < 0.0L,"The refined interior root was moved to the boundary.");
    FELevelSetMoments::LineSliceInsertEvent(1.0L,0.0L,128U,0.0L,events,sources,uncertainties);
    CheckCandidate(events.size() == 3 && sources.back() == 128U,"The independent boundary event was lost.");
    const long double width = (events[2].first-events[1].first)+(events[2].second-events[1].second);
    CheckCandidate(width > 4.0L*(uncertainties[2]+uncertainties[1]),"The positive interval next to the boundary was lost.");
    std::cout << "duplicate root boundary width=" << width << '\n';
  }
  catch(const std::exception& exception) {
    CheckCandidate(false,"A duplicated simple-root candidate was rejected.");
    std::cout << "duplicate candidate rejection=" << exception.what() << '\n';
  }
}

void CheckDistinctCandidates(const std::vector<long double>& polynomial) {
  std::vector<std::pair<long double,long double>> events;
  std::vector<unsigned> sources;
  std::vector<long double> uncertainties;
  try {
    FELevelSetMoments::LineSlicePolynomialEvents(polynomial,2U,events,sources,uncertainties);
    CheckCandidate(events.size() == 2,"Two genuinely distinct roots were consolidated.");
    if(events.size() != 2) return;
    const long double width = (events[1].first-events[0].first)+(events[1].second-events[0].second);
    CheckCandidate(width > 0.0L && sources[0] == 2U && sources[1] == 2U,"A genuine narrow root interval lost its ordering or provenance.");
    std::cout << "distinct root width=" << width << '\n';
  }
  catch(const std::exception& exception) {
    CheckCandidate(false,"A genuinely distinct close-root control was rejected.");
    std::cout << "distinct candidate rejection=" << exception.what() << '\n';
  }
}

void CheckDifferentPolynomials() {
  std::vector<std::pair<long double,long double>> events;
  std::vector<unsigned> sources;
  std::vector<long double> uncertainties;
  FELevelSetMoments::LineSlicePolynomialEvents({-0.5L,1.0L},16U,events,sources,uncertainties);
  FELevelSetMoments::LineSlicePolynomialEvents({-0.5L-std::ldexp(1.0L,-60),1.0L},32U,events,sources,uncertainties);
  CheckCandidate(events.size() == 2 && sources[0] == 16U && sources[1] == 32U,"Different-polynomial root identities were merged.");
}

int main() {
  std::cout << std::setprecision(21);
  CheckDuplicateCandidate();
  CheckDistinctCandidates({0.0L,-1.0e-24L,1.0L});
  CheckDistinctCandidates({0.25L-std::ldexp(1.0L,-60),-1.0L,1.0L});
  CheckDifferentPolynomials();
  std::cout << "candidate checks=" << candidateChecks << " failures=" << candidateFailures << '\n';
  return candidateFailures == 0 ? 0 : 1;
}
