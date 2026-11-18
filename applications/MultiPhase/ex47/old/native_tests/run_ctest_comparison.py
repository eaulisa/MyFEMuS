import argparse
import json
from pathlib import Path
import re
import subprocess
import xml.etree.ElementTree as ET

parser = argparse.ArgumentParser(description='Run the unchanged native CTest commands against both comparison backends.')
parser.add_argument('--build', type=Path, required=True)
parser.add_argument('--output', type=Path, required=True)
parser.add_argument('--ctest', type=Path, default=Path('ctest'))
parser.add_argument('--jobs', type=int, default=2)
parser.add_argument('--detection-depth', type=int, default=12)
args = parser.parse_args()
args.output.mkdir(parents=True, exist_ok=True)
here = Path(__file__).resolve().parent
suite = here.parents[3] / 'unittests/testCutFemWeightLSNative'
cmake = (suite / 'CMakeLists.txt').read_text()
names = re.search(r'set\(FEMUS_NATIVE_TESTS\s+([^)]*)\)', cmake).group(1).split()
commands = [(name, [name]) for name in names]
commands += [('hex_sphere_mesh', ['hex_tet_validation','mesh','HEX','2','4','1e-9']), ('tet_sphere_mesh', ['hex_tet_validation','mesh','TET','2','4','1e-9'])]
summary = {'commands': len(commands), 'detection_depth': args.detection_depth, 'backends': {}}
for backend in ['old', 'prepared']:
    folder = args.build.resolve() / backend
    lines = []
    for name, arguments in commands:
        invocation = [str(folder / arguments[0]), *arguments[1:]]
        lines.append('add_test(cutfem_'+name+' '+ ' '.join(json.dumps(value) for value in invocation)+')')
        lines.append('set_tests_properties(cutfem_'+name+' PROPERTIES WORKING_DIRECTORY '+json.dumps(str(suite))+')')
        if name == 'nodal_cut_detection':
            lines.append('set_tests_properties(cutfem_'+name+' PROPERTIES ENVIRONMENT "EX47_NATIVE_DETECTION_DEPTH='+str(args.detection_depth)+'")')
    (folder / 'CTestTestfile.cmake').write_text('\n'.join(lines)+'\n')
    xml = args.output.resolve() / (backend+'.xml')
    result = subprocess.run([str(args.ctest),'--test-dir',str(folder),'--parallel',str(args.jobs),'--output-on-failure','--output-junit',str(xml)],capture_output=True,text=True)
    (args.output / (backend+'.log')).write_text(result.stdout+result.stderr)
    tests = ET.parse(xml).getroot().findall('testcase')
    failures = [test.attrib['name'] for test in tests if test.find('failure') is not None]
    summary['backends'][backend] = {'exit': result.returncode, 'tests':len(tests), 'failures':failures}
    print(backend,summary['backends'][backend],flush=True)
(args.output / 'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
raise SystemExit(any(row['exit'] or row['tests'] != len(commands) for row in summary['backends'].values()))
