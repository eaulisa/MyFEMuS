#include "Quad9DirectionSelection.hpp"
#include "Tri7DirectionSelection.hpp"
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

// Benchmark stage: derivative-aligned TRI splitting with mapped Fast GK plus optimized QUAD.

struct BenchmarkCase {
  std::string element;
  std::string name;
  unsigned degree;
  std::vector<std::vector<double>> coefficient;
  bool nodal;
  std::vector<double> nodalValues;
};

struct CaseResult {
  bool success;
  double medianSeconds;
  unsigned calls;
  long double checksum;
  long double h00;
  long double d00;
  long double maximumMagnitude;
  std::string error;
};

volatile long double benchmarkSink = 0.0L;

long double WeightedChecksum(const std::vector<std::vector<double>>& positive, const std::vector<std::vector<double>>& delta) {
  long double value = 0.0L;
  unsigned counter = 1;
  for(unsigned i = 0; i < positive.size(); ++i) {
    for(unsigned j = 0; j < positive[i].size(); ++j) {
      value += counter*static_cast<long double>(positive[i][j]);
      value -= (counter+17U)*static_cast<long double>(delta[i][j]);
      ++counter;
    }
  }
  return value;
}

long double MaximumMagnitude(const std::vector<std::vector<double>>& positive, const std::vector<std::vector<double>>& delta) {
  long double value = 0.0L;
  for(unsigned i = 0; i < positive.size(); ++i) {
    for(unsigned j = 0; j < positive[i].size(); ++j) {
      value = std::max(value,std::abs(static_cast<long double>(positive[i][j])));
      value = std::max(value,std::abs(static_cast<long double>(delta[i][j])));
    }
  }
  return value;
}

double PolynomialValue(const std::vector<std::vector<double>>& coefficient, double x, double y) {
  long double value = 0.0L;
  long double xPower = 1.0L;
  for(unsigned i = 0; i < coefficient.size(); ++i) {
    long double yPower = 1.0L;
    for(unsigned j = 0; j < coefficient[i].size(); ++j) {
      value += static_cast<long double>(coefficient[i][j])*xPower*yPower;
      yPower *= y;
    }
    xPower *= x;
  }
  return static_cast<double>(value);
}

std::vector<double> NodalValues(const BenchmarkCase& test) {
  std::vector<double> values;
  if(test.element == "QUAD") {
    static const double point[9][2] = {{-1.0,-1.0},{1.0,-1.0},{1.0,1.0},{-1.0,1.0},{0.0,-1.0},{1.0,0.0},{0.0,1.0},{-1.0,0.0},{0.0,0.0}};
    values.resize(9);
    for(unsigned k = 0; k < 9; ++k) values[k] = PolynomialValue(test.coefficient,point[k][0],point[k][1]);
  }
  else {
    static const double point[7][2] = {{0.0,0.0},{1.0,0.0},{0.0,1.0},{0.5,0.0},{0.5,0.5},{0.0,0.5},{1.0/3.0,1.0/3.0}};
    values.resize(7);
    for(unsigned k = 0; k < 7; ++k) values[k] = PolynomialValue(test.coefficient,point[k][0],point[k][1]);
  }
  return values;
}

void RunIntegration(const BenchmarkCase& test, std::vector<std::vector<double>>& positive, std::vector<std::vector<double>>& delta) {
  if(!test.nodal) {
    if(test.element == "QUAD") FELevelSetMoments::IntegrateQuad9ReferenceSelected(test.coefficient,test.degree,positive,delta,1.0e-10,1.0e-10,true,false,24,0,200000,nullptr,1,nullptr,true);
    else FELevelSetMoments::IntegrateTri7ReferenceSelected(test.coefficient,test.degree,positive,delta,1.0e-10,1.0e-10,true,false,24,0,200000,nullptr,1,nullptr,true);
    return;
  }
  const std::vector<double>& nodalValues = test.nodalValues;
  if(test.element == "QUAD") FELevelSetMoments::IntegrateQuad9NodalSelected(nodalValues,test.degree,positive,delta,1.0e-10,1.0e-10,true,false,1,0,nullptr,true);
  else FELevelSetMoments::IntegrateTri7NodalSelected(nodalValues,test.degree,positive,delta,1.0e-10,1.0e-10,true,false,1,0,nullptr,true);
}

CaseResult Measure(const BenchmarkCase& test) {
  CaseResult result;
  result.success = false;
  result.medianSeconds = 0.0;
  result.calls = 0;
  result.checksum = 0.0L;
  result.h00 = 0.0L;
  result.d00 = 0.0L;
  result.maximumMagnitude = 0.0L;
  try {
    std::vector<std::vector<double>> positive;
    std::vector<std::vector<double>> delta;
    RunIntegration(test,positive,delta);
    result.checksum = WeightedChecksum(positive,delta);
    result.h00 = positive[0][0];
    result.d00 = delta[0][0];
    result.maximumMagnitude = MaximumMagnitude(positive,delta);
    benchmarkSink += result.checksum*1.0e-300L;
    unsigned calls = 1;
    const double target = 0.025;
    while(calls < 16384U) {
      const std::chrono::steady_clock::time_point start = std::chrono::steady_clock::now();
      for(unsigned k = 0; k < calls; ++k) {
        RunIntegration(test,positive,delta);
        benchmarkSink += WeightedChecksum(positive,delta)*1.0e-300L;
      }
      const double seconds = std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
      if(seconds >= target) break;
      calls *= 2U;
    }
    std::vector<double> samples;
    for(unsigned round = 0; round < 7; ++round) {
      const std::chrono::steady_clock::time_point start = std::chrono::steady_clock::now();
      for(unsigned k = 0; k < calls; ++k) {
        RunIntegration(test,positive,delta);
        benchmarkSink += WeightedChecksum(positive,delta)*1.0e-300L;
      }
      samples.push_back(std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count()/calls);
    }
    std::sort(samples.begin(),samples.end());
    result.medianSeconds = samples[samples.size()/2U];
    result.calls = calls;
    result.success = true;
  }
  catch(const std::exception& exception) {
    result.error = exception.what();
  }
  return result;
}

std::vector<std::vector<double>> Quad(double c00, double c10, double c01, double c20, double c11, double c02, double c21 = 0.0, double c12 = 0.0, double c22 = 0.0) {
  std::vector<std::vector<double>> c(3,std::vector<double>(3,0.0));
  c[0][0]=c00; c[1][0]=c10; c[0][1]=c01; c[2][0]=c20; c[1][1]=c11; c[0][2]=c02; c[2][1]=c21; c[1][2]=c12; c[2][2]=c22;
  return c;
}

std::vector<std::vector<double>> Tri(double c00, double c10, double c01, double c20, double c11, double c02, double bubble = 0.0) {
  std::vector<std::vector<double>> c(3,std::vector<double>(3,0.0));
  c[0][0]=c00; c[1][0]=c10; c[0][1]=c01; c[2][0]=c20; c[1][1]=c11; c[0][2]=c02; c[2][1]=bubble; c[1][2]=bubble;
  return c;
}

void AddTimingCases(std::vector<BenchmarkCase>& tests) {
  for(unsigned degreeIndex = 0; degreeIndex < 2; ++degreeIndex) {
    const unsigned degree = degreeIndex == 0 ? 2U : 6U;
    tests.push_back({"QUAD","plane_x_d"+std::to_string(degree),degree,Quad(-0.23,1.0,0.0,0.0,0.0,0.0)});
    tests.push_back({"QUAD","tilted_plane_d"+std::to_string(degree),degree,Quad(-0.15,1.0,0.6,0.0,0.0,0.0)});
    tests.push_back({"QUAD","curved_graph_d"+std::to_string(degree),degree,Quad(-0.12,1.0,0.0,0.0,0.0,-0.28)});
    tests.push_back({"QUAD","circle_d"+std::to_string(degree),degree,Quad(-0.36,0.0,0.0,1.0,0.0,1.0)});
    tests.push_back({"QUAD","small_circle_d"+std::to_string(degree),degree,Quad(0.68*0.68+0.14*0.14-0.07*0.07,-1.36,-0.28,1.0,0.0,1.0)});
    tests.push_back({"QUAD","near_plane_d"+std::to_string(degree),degree,Quad(-0.02,1.0,0.03,0.0,0.0,0.001)});
    tests.push_back({"TRI","plane_x_d"+std::to_string(degree),degree,Tri(-0.28,1.0,0.0,0.0,0.0,0.0)});
    tests.push_back({"TRI","tilted_plane_d"+std::to_string(degree),degree,Tri(-0.38,1.0,0.7,0.0,0.0,0.0)});
    tests.push_back({"TRI","curved_graph_d"+std::to_string(degree),degree,Tri(-0.18,1.0,0.0,0.0,0.0,-0.35)});
    tests.push_back({"TRI","circle_d"+std::to_string(degree),degree,Tri(0.35*0.35+0.35*0.35-0.15*0.15,-0.70,-0.70,1.0,0.0,1.0)});
    tests.push_back({"TRI","small_circle_d"+std::to_string(degree),degree,Tri(0.72*0.72+0.18*0.18-0.045*0.045,-1.44,-0.36,1.0,0.0,1.0)});
    tests.push_back({"TRI","bubble_rich_d"+std::to_string(degree),degree,Tri(-0.08,0.9,0.4,-0.25,0.35,-0.15,-0.22)});
  }
}

void DumpMoments(const BenchmarkCase& test, const std::vector<std::vector<double>>& positive, const std::vector<std::vector<double>>& delta) {
  for(unsigned i = 0; i < positive.size(); ++i) for(unsigned j = 0; j < positive[i].size(); ++j) std::cout << "MOMENT," << test.element << "," << test.name << "," << test.degree << "," << i << "," << j << "," << positive[i][j] << "," << delta[i][j] << "\n";
}

void RunStress(const std::string& element) {
  std::mt19937_64 generator(element == "QUAD" ? 918273645ULL : 564738291ULL);
  std::uniform_real_distribution<double> distribution(-0.45,0.45);
  for(unsigned testIndex = 0; testIndex < 120; ++testIndex) {
    std::vector<std::vector<double>> coefficient(3,std::vector<double>(3,0.0));
    coefficient[0][0] = distribution(generator);
    coefficient[1][0] = distribution(generator);
    coefficient[0][1] = distribution(generator);
    coefficient[2][0] = 0.35*distribution(generator);
    coefficient[1][1] = 0.35*distribution(generator);
    coefficient[0][2] = 0.35*distribution(generator);
    if(element == "QUAD") {
      coefficient[2][1] = 0.12*distribution(generator);
      coefficient[1][2] = 0.12*distribution(generator);
      coefficient[2][2] = 0.04*distribution(generator);
    }
    else {
      const double bubble = 0.12*distribution(generator);
      coefficient[2][1] = bubble;
      coefficient[1][2] = bubble;
    }
    BenchmarkCase test = {element,"stress_"+std::to_string(testIndex),4,coefficient,false};
    try {
      std::vector<std::vector<double>> positive;
      std::vector<std::vector<double>> delta;
      RunIntegration(test,positive,delta);
      DumpMoments(test,positive,delta);
      std::cout << "STRESS," << element << "," << testIndex << ",OK," << std::setprecision(18) << static_cast<double>(WeightedChecksum(positive,delta)) << "," << positive[0][0] << "," << delta[0][0] << "," << static_cast<double>(MaximumMagnitude(positive,delta)) << "\n";
    }
    catch(const std::exception& exception) {
      std::string error = exception.what();
      std::replace(error.begin(),error.end(),',',';');
      std::cout << "STRESS," << element << "," << testIndex << ",FAIL,0,0,0,0," << error << "\n";
    }
  }
}

int main() {
  std::cout << std::setprecision(18);
  std::vector<BenchmarkCase> tests;
  AddTimingCases(tests);
  const unsigned coefficientCases = static_cast<unsigned>(tests.size());
  for(unsigned k = 0; k < coefficientCases; ++k) {
    BenchmarkCase nodal = tests[k];
    nodal.name = "nodal_"+nodal.name;
    nodal.nodal = true;
    nodal.nodalValues = NodalValues(nodal);
    tests.push_back(nodal);
    BenchmarkCase matched = tests[k];
    matched.name = "reconstructed_"+matched.name;
    if(matched.element == "QUAD") FELevelSetMoments::Quad9ReferenceCoefficientsFromNodalValues(nodal.nodalValues,matched.coefficient);
    else FELevelSetMoments::Tri7ReferenceCoefficientsFromNodalValues(nodal.nodalValues,matched.coefficient);
    tests.push_back(matched);
  }
  for(unsigned k = 0; k < tests.size(); ++k) {
    const CaseResult result = Measure(tests[k]);
    if(result.success) {
      std::vector<std::vector<double>> positive;
      std::vector<std::vector<double>> delta;
      RunIntegration(tests[k],positive,delta);
      DumpMoments(tests[k],positive,delta);
    }
    if(result.success) std::cout << "TIMING," << tests[k].element << "," << tests[k].name << "," << tests[k].degree << ",OK," << result.medianSeconds << "," << result.calls << "," << static_cast<double>(result.checksum) << "," << static_cast<double>(result.h00) << "," << static_cast<double>(result.d00) << "," << static_cast<double>(result.maximumMagnitude) << "\n";
    else {
      std::string error = result.error;
      std::replace(error.begin(),error.end(),',',';');
      std::cout << "TIMING," << tests[k].element << "," << tests[k].name << "," << tests[k].degree << ",FAIL,0,0,0,0,0,0," << error << "\n";
    }
  }
  RunStress("QUAD");
  RunStress("TRI");
  return benchmarkSink == 123456789.0L ? 1 : 0;
}
