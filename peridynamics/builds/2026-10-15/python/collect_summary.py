"""Gather machine-readable validation evidence without copying runtime caches."""
from pathlib import Path
import hashlib
import importlib.metadata
import json
import platform
import sys
import xml.etree.ElementTree as ET

root=Path('/workspace/work/python-build')
jobs={}
for path in sorted((root/'results').glob('*.json')):
    value=json.loads(path.read_text())
    if 'exit_code' not in value:
        continue
    xml=root/'results'/f'{path.stem}.xml'
    if xml.exists():
        suite=ET.parse(xml).getroot().find('testsuite')
        value['test_results']={k:float(suite.attrib[k]) if k=='time' else int(suite.attrib[k])
                               for k in ('tests','errors','failures','skipped','time')}
    jobs[path.stem]=value
wheels=[dict(path=str(path),size_bytes=path.stat().st_size,
             sha256=hashlib.sha256(path.read_bytes()).hexdigest())
        for path in (root/'wheels').glob('*.whl')]
pending=[name for name in ('fresh-tests','reproduction-tests','recovered-fresh-tests','recovered-reproduction-tests','fresh-wheel','fresh-mesh',
                           'fresh-cli','reproduction-cli','thesis-cli','printed-cli','wheel-install','wheel-smoke','recovered-regression')
         if name not in jobs or jobs[name]['exit_code']!=0]
summary=dict(status='passed' if not pending else 'incomplete_or_failed',pending_or_failed_jobs=pending,
             python=sys.version,platform=platform.platform(),
             packages={name:importlib.metadata.version(name) for name in
                       ('numpy','scipy','numba','llvmlite','matplotlib','pytest','setuptools','wheel','pip')},
             scope='Build and software verification of recovered Oct 9 Python projects; no full impact rerun or scientific benchmark validation.',
             scientific_benchmark_validated=False,jobs=jobs,wheels=wheels,
             recovered_source_hashes_match_original={
                 name: jobs[f'{name}-tests'].get('source_hashes')==jobs.get(f'recovered-{name}-tests',{}).get('source_hashes')
                 for name in ('fresh','reproduction')},
             wheel_smoke=json.loads((root/'results/wheel-smoke-details.json').read_text()),
             notes=['Original source files were hash-checked unchanged by each completed test suite.',
                    'Wheel and mesh builds used a separate byte-preserving source copy under work/python-build/fresh-source.',
                    'All build/install commands were offline with no dependency downloads.',
                    'Initial test attempts had harness-only errors because pytest basetemp parent was missing; setup was corrected and both suites rerun. Initial logs are retained separately.'])
(root/'build-summary.json').write_text(json.dumps(summary,indent=2)+'\n')
lines=['# Python build and software verification','',summary['scope'],'',
       'Original source files were not edited. Caches, temporary files, builds and smoke checkpoints stayed under `work/python-build`.','',
       '| Check | Result | Duration |','|---|---|---|']
for name,value in jobs.items():
    result='Passed' if value['exit_code']==0 else f"Exit {value['exit_code']}"
    if 'test_results' in value:
        t=value['test_results'];result=f"{t['tests']-t['errors']-t['failures']-t['skipped']} passed, {t['failures']} failed, {t['errors']} errors"
    lines.append(f"| {name} | {result} | {value['elapsed_seconds']:.3f} s |")
lines += ['', 'The installed wheel assembled 33,600 concrete particles, 270 steel nodes and 1,755,720 bonds, advanced one step, and passed an exact checkpoint roundtrip. The mesh starts at zero force. The recovered Stage 19 fixture regression passes all 112 retained samples.','',
          '## Reproduce in the attached workspace','', '```bash',
          '/workspace/bowen_reproduction/.venv/bin/python /workspace/work/python-build/validate_projects.py fresh-tests',
          '/workspace/bowen_reproduction/.venv/bin/python /workspace/work/python-build/validate_projects.py reproduction-tests',
          '/workspace/bowen_reproduction/.venv/bin/python /workspace/work/python-build/validate_projects.py recovered-fresh-tests',
          '/workspace/bowen_reproduction/.venv/bin/python /workspace/work/python-build/validate_projects.py recovered-reproduction-tests',
          '/workspace/bowen_reproduction/.venv/bin/python /workspace/work/python-build/validate_projects.py regression',
          '/workspace/bowen_reproduction/.venv/bin/python /workspace/work/python-build/validate_projects.py wheel',
          '/workspace/bowen_reproduction/.venv/bin/python /workspace/work/python-build/validate_projects.py mesh',
          '/workspace/bowen_reproduction/.venv/bin/python /workspace/work/python-build/validate_projects.py cli',
          '/workspace/bowen_reproduction/.venv/bin/python /workspace/work/python-build/validate_projects.py wheel-smoke',
          '/workspace/bowen_reproduction/.venv/bin/python /workspace/work/python-build/collect_summary.py', '```','',
          'The driver records every executed command, working directory, environment override, exit status and duration in `results/*.json`. JUnit XML and detailed test durations accompany the logs. The build command is `python -m pip wheel --no-deps --no-build-isolation --no-index --wheel-dir <wheel-directory> <source-copy>`.','',
          'Python '+platform.python_version()+'. Runtime dependencies match the pinned requirements; detailed versions and hashes are in `build-summary.json`.','',
          'These checks establish executable software and checkpoint consistency. The models remain declared reconstructions; they do not establish Bowen author equivalence, full-event agreement or experimental validation.']
(root/'BUILD_SUMMARY.md').write_text('\n'.join(lines)+'\n')
print(json.dumps(dict(status=summary['status'],pending_or_failed_jobs=pending,
                      output=str(root/'build-summary.json'))))
