from pathlib import Path

root = Path(__file__).resolve().parents[2]
folder = root / 'src/quadrature/cutfem_ls'
suite = Path(__file__).resolve().parent


def once(text, old, new):
    if old not in text:
        raise RuntimeError('Missing follow-up anchor: ' + old[:100])
    return text.replace(old, new, 1)


path = folder / 'LineDomainIntegration.hpp'
text = path.read_text()
helper = r'''  static bool HasVerticalFiber(const std::vector<std::vector<long double>>& coefficient, unsigned order) {
    std::vector<std::vector<long double>> fibers(3,std::vector<long double>(3,0.0L));
    for(unsigned j = 0; j < 3; ++j) {
      for(unsigned i = 0; i < 3; ++i) fibers[j][i] = order == 1 ? coefficient[i][j] : coefficient[j][i];
      LineSliceUtilities::Trim(fibers[j]);
      if(fibers[j].size() == 1 && fibers[j][0] != 0.0L) return false;
    }
    for(unsigned j = 0; j < 3; ++j) {
      if(fibers[j].size() < 2) continue;
      std::vector<long double> roots;
      LineSliceUtilities::PolynomialRoots(fibers[j],0.0L,1.0L,roots);
      for(unsigned r = 0; r < roots.size(); ++r) {
        bool common = true;
        for(unsigned k = 0; k < 3; ++k) {
          const long double residual = std::abs(LineSliceUtilities::PolynomialValue(fibers[k],roots[r]));
          const long double rounding = 128.0L*std::numeric_limits<long double>::epsilon()*LineSliceUtilities::PolynomialBound(fibers[k],roots[r]);
          common = common && residual <= rounding;
        }
        if(common) return true;
      }
      return false;
    }
    return false;
  }

'''
text = once(text, '  static void IntegrateCell(', helper+'  static void IntegrateCell(')
needle = '    std::vector<std::vector<long double>> working(coefficient);'
replacement = r'''    if(computeDelta && HasVerticalFiber(coefficient,order)) {
      LineSliceDirection(coefficient,degree,triangle,order,computeDelta,false,absoluteH,absoluteD,relative,maximumDepth,maximumEvaluations,values,errors,evaluations,accepted,0,nullptr,surfaceGradient);
      if(order == 2) {
        const unsigned count = (degree+1U)*(degree+2U)/2U;
        for(unsigned i = 0; i <= degree; ++i) for(unsigned j = i+1U; j <= degree-i; ++j) {
          const unsigned a = LineSliceUtilities::MomentIndex(degree,i,j);
          const unsigned b = LineSliceUtilities::MomentIndex(degree,j,i);
          for(unsigned pair = 0; pair < 2; ++pair) {
            std::swap(values[pair*count+a],values[pair*count+b]);
            std::swap(errors[pair*count+a],errors[pair*count+b]);
          }
        }
      }
      return;
    }
    std::vector<std::vector<long double>> working(coefficient);'''
text = once(text, needle, replacement)
path.write_text(text)
for name, old, new in [('Quad9DirectionSelection.hpp','bool useLine = !useTaylor || zeroEdge;','bool useLine = !useTaylor || zeroEdge || fullBoundaryDelta;'),('Tri7DirectionSelection.hpp','bool line = !useTaylor || zeroEdge;','bool line = !useTaylor || zeroEdge || fullBoundaryDelta;')]:
    path = folder / name
    path.write_text(once(path.read_text(), old, new))

# Emit every moment, and separate input reconstruction from timed nodal evaluation.
path = suite / 'benchmark_quad_tri.cpp'
text = path.read_text()
text = once(text, '  bool nodal;\n};', '  bool nodal;\n  std::vector<double> nodalValues;\n};')
text = once(text, '  const std::vector<double> nodalValues = NodalValues(test);', '  const std::vector<double>& nodalValues = test.nodalValues;')
text = text.replace('true,false,28,0,200000,nullptr,0,nullptr,true);', 'true,false,24,0,200000,nullptr,1,nullptr,true);')
helper = r'''void DumpMoments(const BenchmarkCase& test, const std::vector<std::vector<double>>& positive, const std::vector<std::vector<double>>& delta) {
  for(unsigned i = 0; i < positive.size(); ++i) for(unsigned j = 0; j < positive[i].size(); ++j) std::cout << "MOMENT," << test.element << "," << test.name << "," << test.degree << "," << i << "," << j << "," << positive[i][j] << "," << delta[i][j] << "\n";
}

'''
text = once(text, 'void RunStress(', helper+'void RunStress(')
text = once(text, '      RunIntegration(test,positive,delta);\n      std::cout << "STRESS,', '      RunIntegration(test,positive,delta);\n      DumpMoments(test,positive,delta);\n      std::cout << "STRESS,')
text = once(text, '    nodal.nodal = true;\n    tests.push_back(nodal);', '    nodal.nodal = true;\n    nodal.nodalValues = NodalValues(nodal);\n    tests.push_back(nodal);\n    BenchmarkCase matched = tests[k];\n    matched.name = "reconstructed_"+matched.name;\n    if(matched.element == "QUAD") FELevelSetMoments::Quad9ReferenceCoefficientsFromNodalValues(nodal.nodalValues,matched.coefficient);\n    else FELevelSetMoments::Tri7ReferenceCoefficientsFromNodalValues(nodal.nodalValues,matched.coefficient);\n    tests.push_back(matched);')
text = once(text, '    const CaseResult result = Measure(tests[k]);', '    const CaseResult result = Measure(tests[k]);\n    if(result.success) {\n      std::vector<std::vector<double>> positive;\n      std::vector<std::vector<double>> delta;\n      RunIntegration(tests[k],positive,delta);\n      DumpMoments(tests[k],positive,delta);\n    }')
path.write_text(text)

path = suite / 'build.py'
text = path.read_text()
text = once(text, "'geometry': ['line_integration'", "'geometry': ['review_regressions', 'line_integration'")
path.write_text(text)
print('Degenerate fibers retained; benchmark now records every moment and matched input paths.')
