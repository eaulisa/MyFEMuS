import argparse
import concurrent.futures
import hashlib
import json
import os
from pathlib import Path
import re
import subprocess
import time

here = Path(__file__).resolve().parent
repo = here.parents[3]
suite = repo / 'unittests/testCutFemWeightLSNative'
parser = argparse.ArgumentParser()
parser.add_argument('--build', type=Path, required=True)
parser.add_argument('--tests', nargs='*')
parser.add_argument('--jobs', type=int, default=2)
parser.add_argument('--backend', choices=['old','prepared','both'], default='both')
parser.add_argument('--boost', type=Path, required=True)
parser.add_argument('--eigen', type=Path, required=True)
args = parser.parse_args()
args.build.mkdir(parents=True,exist_ok=True)
support = repo / 'unittests/testCutFemWeightLSNative/build'
objects = [support / (name+'.o') for name in ['Line','Quadrilateral','Triangle','Hexaedron','Tetrahedral','Wedge']]
objects += [support/(p.stem+'.o') for p in sorted((repo/'src/quadrature').glob('*Points.cpp'))]
includes = [here,repo/'src/quadrature/cutfem_ls',repo/'src/quadrature',repo/'src/quadrature/cutfem',repo/'src/fe',repo/'src/enums',args.boost,args.eigen]
flags = ['-std=c++14','-O2','-DNDEBUG',*['-I'+str(p) for p in includes]]
if not all(p.exists() for p in objects):
    support=args.build/'support'
    support.mkdir(exist_ok=True)
    support_sources=[repo/'src/fe'/(name+'.cpp') for name in ['Line','Quadrilateral','Triangle','Hexaedron','Tetrahedral','Wedge']]+sorted((repo/'src/quadrature').glob('*Points.cpp'))
    def build_support(source):
        target=support/(source.stem+'.o')
        command=[os.environ.get('CXX','c++'),*flags,'-c',str(source),'-o',str(target)]
        subprocess.run(command,check=True,capture_output=True,text=True,timeout=300)
        return target
    with concurrent.futures.ThreadPoolExecutor(max_workers=args.jobs) as pool:
        objects=list(pool.map(build_support,support_sources))
all_names = sorted(p.stem for p in suite.glob('*.cpp'))
names = args.tests if args.tests else all_names
assert set(names)<=set(all_names+['native_suite_speed','curved_public_regression','benchmark_graph_accuracy','benchmark_surface_accuracy'])
backends = ['old','prepared'] if args.backend=='both' else [args.backend]

def generate(name,backend):
    original={'native_suite_speed':'native_suite','benchmark_graph_accuracy':'benchmark_graph_dimensions','benchmark_surface_accuracy':'benchmark_surface_reuse'}.get(name,name)
    source_path=here/(name+'.cpp') if name=='curved_public_regression' else suite/(original+'.cpp')
    source=source_path.read_text()
    if name=='benchmark_graph_accuracy':
        source=source.replace('double DimensionIntegrate(','std::vector<long double> dimensionLastMoments;\n\ndouble DimensionIntegrate(',1)
        capture='  dimensionLastMoments.clear();\n  for(unsigned pair = 0; pair < 2; ++pair) for(unsigned i = 0; i <= degree; ++i) for(unsigned j = 0; j <= degree-i; ++j) { if(element < 2) dimensionLastMoments.push_back(pair == 0 ? positive2[i][j] : surface2[i][j]); else for(unsigned k = 0; k <= degree-i-j; ++k) dimensionLastMoments.push_back(pair == 0 ? positive3[i][j][k] : surface3[i][j][k]); }\n'
        source=source.replace('  double checksum = 0.0;',capture+'  double checksum = 0.0;',1)
        source=source[:source.index('int main(int argc')]+'''int main(int argc, char** argv) {
  if(argc != 2) return 64;
  const unsigned degree = std::strtoul(argv[1],nullptr,10);
  std::cout << std::setprecision(21);
  for(unsigned index = 0; index < 10; ++index) for(unsigned element = 0; element < 4; ++element) { DimensionIntegrate(element,DimensionNodal(element,index),degree); std::cout << element << ',' << index; for(long double value : dimensionLastMoments) std::cout << ',' << value; std::cout << '\\n'; }
}
'''
    if name=='benchmark_surface_accuracy':
        helper='''template<class Weight> void OutputSurfaceMoments(Weight& rule, unsigned mode, unsigned cell, const std::vector<double>& positive, const std::vector<double>& tangential, const std::vector<std::vector<double>>& normal) {
  const unsigned count = rule.GetGaussQuadraturePointNumber();
  std::vector<std::array<long double,14>> first(count), second(count), third(count);
  for(unsigned n = 0; n < count; ++n) { first[n][0] = second[n][0] = third[n][0] = 1; for(unsigned power = 1; power <= 13; ++power) { first[n][power] = first[n][power-1]*(static_cast<long double>(rule.GetGaussCoordinatePointer(0)[n])+1)/2; second[n][power] = second[n][power-1]*(static_cast<long double>(rule.GetGaussCoordinatePointer(1)[n])+1)/2; third[n][power] = third[n][power-1]*(static_cast<long double>(rule.GetGaussCoordinatePointer(2)[n])+1)/2; } }
  std::cout << std::setprecision(21);
  for(unsigned i = 0; i <= 13; ++i) for(unsigned j = 0; j <= 13-i; ++j) for(unsigned k = 0; k <= 13-i-j; ++k) { long double sum[5] = {0,0,0,0,0}; for(unsigned n = 0; n < count; ++n) { const long double factor = rule.GetGaussWeightPointer()[n]*first[n][i]*second[n][j]*third[n][k]; sum[0] += factor*positive[n]; sum[1] += factor*tangential[n]; for(unsigned component = 0; component < 3; ++component) sum[2+component] += factor*normal[component][n]; } std::cout << "M," << mode << ',' << cell << ',' << i << ',' << j << ',' << k; for(long double value : sum) std::cout << ',' << value; std::cout << '\\n'; }
}

'''
        source=source.replace('int main() {',helper+'int main() {',1)
        call='      rule.GetNormalSurfaceWeightsFromNodalGeometry(nodal,coordinates,positive,normal,1e-8,1e-8);'
        source=source.replace(call,call+'\n      OutputSurfaceMoments(rule,mode,cell,positive,tangential,normal);',1)
    if name=='native_suite_speed':
        source=source.replace('if(argc != 8) return 64;','if(argc != 9) return 64;\n  const unsigned repetitions = std::strtoul(argv[8],nullptr,10);\n  if(repetitions == 0 || repetitions > 10000) return 64;')
        source=re.sub(r'\brule\.(?:GetWeightsFromNodalValues|GetPhysicalSurfaceWeightsFromNodalGeometry)\([^;]+;',r'{ for(unsigned repeat = 0; repeat < repetitions; ++repeat) { rule.ClearCache(); \g<0> } }',source)
    source=re.sub(r'\bCutFemWeight_LS\b(?=\s*<)','Ex47NativeTests::'+('OriginalWeight' if backend=='old' else 'PreparedWeight'),source)
    source=re.sub(r'\bCutFemWeightRules_LS\b(?=\s*<)','Ex47NativeTests::'+('OriginalRules' if backend=='old' else 'PreparedRules'),source)
    preface='#include "NativeCompatibility.hpp"\n'
    if name=='test_order_planar_first': preface='#include "CutFemWeight.hpp"\n'+preface
    declarations=(here/'NativeCompatibility.hpp').read_text()
    functions=re.findall(r'^(?:void|bool|int) (\w+)\(',declarations,re.M)
    for symbol in functions+['LineIntegration','TriangleIntegration','PlanarEndpointIntegration','SurfaceEndpointIntegration']:
        source=re.sub(r'(?<![\w:])(?:FELevelSetMoments::)?'+symbol+r'\b', 'Ex47NativeTests::'+symbol,source)
    source=source.replace('FELevelSetMoments::LineDomainIntegration::','Ex47NativeTests::LineIntegration::')
    source=re.sub(r'(#include\s+")Ex47NativeTests::',r'\1',source)
    if backend=='prepared':
        for symbol in ['SurfaceEndpointStatistics','LineIntegrationStatistics','EndpointStripStatistics','EndpointStripQuadrature','EndpointProjection']:
            source=source.replace('FELevelSetMoments::'+symbol,'FELevelSetMoments::Ex47Prepared::'+symbol)
    return preface+source

def compatibility(backend):
    source=here/('OriginalCompatibility.cpp' if backend=='old' else 'NativeCompatibility.cpp')
    target=args.build/(backend+'_compatibility.o')
    command=[os.environ.get('CXX','c++'),*flags,'-DEX47_NATIVE_USE_OLD' if backend=='old' else '-DEX47_NATIVE_PREPARED','-c',str(source),'-o',str(target)]
    dependencies=[source,*here.glob('*.hpp'),*(repo/'applications/MultiPhase/ex47/include').rglob('*.hpp')]
    fingerprint=hashlib.sha256((json.dumps(command)+''.join(hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(dependencies))).encode()).hexdigest()
    stamp=target.with_suffix('.sha256')
    if not target.exists() or not stamp.exists() or stamp.read_text()!=fingerprint:
        result=subprocess.run(command,capture_output=True,text=True,timeout=600)
        target.with_suffix('.build.log').write_text(result.stdout+result.stderr)
        if result.returncode: raise RuntimeError(target.with_suffix('.build.log').read_text()[-12000:])
        stamp.write_text(fingerprint)
    return target

compatibility_objects={backend:compatibility(backend) for backend in backends}

def build(item):
    name,backend=item
    folder=args.build/backend
    folder.mkdir(exist_ok=True)
    source=folder/(name+'.cpp')
    text=generate(name,backend)
    if not source.exists() or source.read_text()!=text: source.write_text(text)
    binary=folder/name
    command=[os.environ.get('CXX','c++'),*flags,'-DEX47_NATIVE_USE_OLD' if backend=='old' else '-DEX47_NATIVE_PREPARED']
    if name=='endpoint_comparison': command.append('-DFEMUS_ENDPOINT_ONLY')
    object_file=folder/(name+'.o')
    compile_command=command+['-c',str(source),'-o',str(object_file)]
    command=[os.environ.get('CXX','c++'),str(object_file),str(compatibility_objects[backend]),*[str(p) for p in objects],'-o',str(binary)]
    dependencies=[*here.glob('*.hpp'),*(repo/'applications/MultiPhase/ex47/include').rglob('*.hpp')]
    source_fingerprint=hashlib.sha256((text+json.dumps(compile_command)+''.join(hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(dependencies))).encode()).hexdigest()
    fingerprint=hashlib.sha256((source_fingerprint+json.dumps(command)+hashlib.sha256(compatibility_objects[backend].read_bytes()).hexdigest()).encode()).hexdigest()
    stamp=folder/(name+'.sha256')
    if binary.exists() and stamp.exists() and stamp.read_text()==fingerprint: return {'test':name,'backend':backend,'exit':0,'cached':True}
    started=time.monotonic()
    try:
        object_stamp=object_file.with_suffix('.object.sha256')
        if not object_file.exists() or not object_stamp.exists() or object_stamp.read_text()!=source_fingerprint:
            result=subprocess.run(compile_command,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True,timeout=600)
            (folder/(name+'.build.log')).write_text(result.stdout)
            if result.returncode:
                record={'test':name,'backend':backend,'exit':result.returncode,'seconds':time.monotonic()-started}
                print(json.dumps(record),flush=True)
                return record
            object_stamp.write_text(source_fingerprint)
        result=subprocess.run(command,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True,timeout=600)
        (folder/(name+'.build.log')).write_text(result.stdout)
        if result.returncode==0: stamp.write_text(fingerprint)
        record={'test':name,'backend':backend,'exit':result.returncode,'seconds':time.monotonic()-started}
    except subprocess.TimeoutExpired:
        record={'test':name,'backend':backend,'exit':124,'seconds':time.monotonic()-started}
    print(json.dumps(record),flush=True)
    return record

with concurrent.futures.ThreadPoolExecutor(max_workers=args.jobs) as pool:
    records=list(pool.map(build,[(name,backend) for name in names for backend in backends]))
(args.build/'build_results.json').write_text(json.dumps(records,indent=2)+'\n')
(args.build/'source_manifest.json').write_text(json.dumps({str(p.relative_to(repo)):hashlib.sha256(p.read_bytes()).hexdigest() for p in [*suite.glob('*.cpp'),*here.glob('*.hpp'),*here.glob('*.cpp'),*(repo/'applications/MultiPhase/ex47/include').rglob('*.hpp'),*(repo/'src/quadrature/cutfem_ls').glob('*.hpp')]},indent=2)+'\n')
raise SystemExit(any(record['exit'] for record in records))
