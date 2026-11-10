from pathlib import Path
import hashlib
import json
import subprocess

suite=Path(__file__).resolve().parent
repo=suite.parents[1]
manifest={str(p.relative_to(repo)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted((repo/'src/quadrature/cutfem_ls').glob('*.hpp'))}
reports={}
for name in ['tangency_strip_experiment','tangency_preservation','tangency_sanitizers','tangent_quad_summary','tangent_tri_summary','tangent_quad_mesh_summary','tangent_tri_mesh_summary','tangent_double_summary']:
    report=json.loads((suite/'review_results'/(name+'.json')).read_text())
    assert report['accepted'] and report['source_sha256']==manifest,name
    reports[name]=report
original=json.loads((suite/'candidate_original_summary.json').read_text())
supplemental=json.loads((suite/'candidate_supplemental_summary.json').read_text())
focused=json.loads((suite/'extended_results.json').read_text())
planar=json.loads((suite/'planar_comparison_results.json').read_text())['summary']
build=json.loads((suite/'build_results.json').read_text())
assert json.loads((suite/'build/source_manifest.json').read_text())==manifest
assert original['cases']==386 and original['flagged_cases']==0
assert supplemental['cases']==8 and supplemental['flagged_cases']==0
assert len(focused)==19 and all(r['exit']==0 for r in focused)
assert planar['cases']==147 and planar['flagged_cases']==0
assert len(build)==30 and all(r['exit']==0 for r in build)
trial=reports['tangency_strip_experiment']
assert trial['cases']==405 and trial['graph_cases_bitwise']==180
preservation=reports['tangency_preservation']
assert len(preservation['three_dimensional_cases'])==16 and all(r['unchanged'] for r in preservation['three_dimensional_cases'])
result={'accepted':True,'input_sha':subprocess.check_output(['git','rev-parse','HEAD'],cwd=repo,text=True).strip(),'baseline_sha':trial['baseline_sha'],'branch':'rifat-cutfem-2d-tangent-strips','maximum_requested_degree':10,'experiment_cases':405,'independent_experiment_checks':trial['checks'],'bitwise_graph_cases':180,'three_dimensional_cases_unchanged':16,'three_dimensional_checks':preservation['three_dimensional_checks'],'wide_cases':sum(reports[n]['cases'] for n in ['tangent_quad_summary','tangent_tri_summary']),'wide_checks':sum(reports[n]['checks'] for n in ['tangent_quad_summary','tangent_tri_summary']),'mesh_levels':sum(reports[n]['mesh_levels'] for n in ['tangent_quad_mesh_summary','tangent_tri_mesh_summary']),'sanitizer_cases':reports['tangency_sanitizers']['cases'],'double_output_cases':reports['tangent_double_summary']['cases'],'native_cases':original['cases']+supplemental['cases'],'native_checks':original['assertions']+supplemental['assertions'],'focused_executions':19,'planar_comparison':planar,'new_path_calls':sum(r['stats']['attempts'] for r in trial['records']),'new_path_completed':sum(r['stats']['completed'] for r in trial['records']),'legacy_fallbacks':sum(r['stats']['rejected'] for r in trial['records']),'source_sha256':manifest,'source_digest':hashlib.sha256(json.dumps(manifest,sort_keys=True).encode()).hexdigest(),'scope':'2D non-graph integration experiment; existing graph algorithms and 3D algorithms retained; full MPI/PETSc ex44 was not run.','performance_decision':'Experimental branch; accuracy passes, but the quadratic-strip method is not a uniform speed improvement.'}
(suite/'review_results/tangency_acceptance.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({k:v for k,v in result.items() if k!='source_sha256'},indent=2))
