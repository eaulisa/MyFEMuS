#!/usr/bin/env python3
import argparse
from pathlib import Path
import subprocess

BASELINE = 'f5f81a8a088cb990780f326c4a8692ac98265536'
FIELDS = ('Regular15','Mapped15','Regular31','Mapped31','Subdivisions','Collapses','Emergency7','Emergency3')


def guarded(statement):
    return '#ifdef FEMUS_BASELINE_OUTER_COUNTERS\n'+statement+'\n#endif\n'


def main():
    parser = argparse.ArgumentParser(description='Add counting-only instrumentation to the pinned baseline checkout.')
    parser.add_argument('checkout',type=Path)
    args = parser.parse_args()
    root = args.checkout.resolve()
    if subprocess.check_output(['git','-C',str(root),'rev-parse','HEAD'],text=True).strip() != BASELINE:
        raise RuntimeError('The checkout must be at the accepted f5f81a8 baseline.')
    for filename,structure,signature in [('LineDomainIntegration.hpp','LineIntegrationStatistics','  static bool GaussKronrodPanel('),('EndpointStripQuadrature.hpp','EndpointStripStatistics','  template<unsigned Rule> static void Panel('),('SurfaceEndpointIntegration.hpp','SurfaceEndpointStatistics','  template<unsigned Rule> static void Panel(')]:
        relative = 'src/quadrature/cutfem_ls/'+filename
        original = subprocess.check_output(['git','-C',str(root),'show',BASELINE+':'+relative],text=True)
        text = original
        field_text = ''.join('  unsigned long long benchmark'+name+' = 0;\n' for name in FIELDS)
        text = text.replace('struct '+structure+' {\n','struct '+structure+' {\n'+guarded(field_text.rstrip()))
        offset = text.index('\n',text.index(signature))+1
        if filename == 'LineDomainIntegration.hpp':
            increment = '    if(mapped) ++Statistics().benchmarkMapped15; else ++Statistics().benchmarkRegular15;'
            lines = text.splitlines(keepends=True)
            text = ''.join((guarded('    ++Statistics().benchmarkSubdivisions;') if '      Adaptive(coefficient,degree,strip,left,right,parameterLeft,middle,' in line or '    Adaptive(coefficient,degree,strip,left,right,parameterLeft,middle,' in line else '')+line for line in lines)
            offset = text.index('\n',text.index(signature))+1
        elif filename == 'EndpointStripQuadrature.hpp':
            increment = '    if(Rule == 15) { if(interval.mapping == 0) ++Statistics().benchmarkRegular15; else ++Statistics().benchmarkMapped15; }\n    if(Rule == 31) { if(interval.mapping == 0) ++Statistics().benchmarkRegular31; else ++Statistics().benchmarkMapped31; }'
        else:
            increment = '    if(Rule == 15) { if(mapping == 0) ++Statistics().benchmarkRegular15; else ++Statistics().benchmarkMapped15; }\n    if(Rule == 31) { if(mapping == 0) ++Statistics().benchmarkRegular31; else ++Statistics().benchmarkMapped31; }\n    if(Rule == 7) ++Statistics().benchmarkEmergency7;\n    if(Rule == 3) ++Statistics().benchmarkEmergency3;'
            text = text.replace('    if(!ResolvedPanel<15>(left,right,mapping,endpoint,delta,exponent)) {\n','    if(!ResolvedPanel<15>(left,right,mapping,endpoint,delta,exponent)) {\n'+guarded('      ++Statistics().benchmarkCollapses;'))
            text = text.replace('    Adaptive(data,cell,outer,left,middle,',guarded('    ++Statistics().benchmarkSubdivisions;')+'    Adaptive(data,cell,outer,left,middle,')
        text = text[:offset]+guarded(increment)+text[offset:]
        path = root/relative
        if path.read_text() not in (original,text):
            raise RuntimeError('Refusing to overwrite unrelated changes: '+relative)
        path.write_text(text)


if __name__ == '__main__':
    main()
