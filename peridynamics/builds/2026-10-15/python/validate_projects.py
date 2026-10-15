"""Validate recovered Python sources without writing into source directories."""
from __future__ import annotations
import argparse
import hashlib
import importlib.metadata
import json
import os
from pathlib import Path
import shutil
import subprocess
import time

ROOT = Path('/workspace/work/python-build')
PYTHON = Path('/workspace/bowen_reproduction/.venv/bin/python')
FRESH = Path('/workspace/bowen_from_scratch')
REPRODUCTION = Path('/workspace/bowen_reproduction')
ENV = dict(os.environ, PYTHONDONTWRITEBYTECODE='1', NUMBA_CACHE_DIR=str(ROOT/'cache/numba'),
           MPLCONFIGDIR=str(ROOT/'cache/matplotlib'), XDG_CACHE_HOME=str(ROOT/'cache'),
           PIP_CACHE_DIR=str(ROOT/'cache/pip'), OPENBLAS_NUM_THREADS='1', OMP_NUM_THREADS='1',
           MKL_NUM_THREADS='1', MPLBACKEND='Agg')
for folder in ('logs', 'results', 'cache', 'wheels', 'tmp'):
    (ROOT/folder).mkdir(parents=True, exist_ok=True)

def source_hashes(root):
    return {str(p.relative_to(root)): hashlib.sha256(p.read_bytes()).hexdigest()
            for p in sorted(root.rglob('*'))
            if p.is_file() and p.suffix in ('.py', '.toml', '.txt', '.md')
            and not any(x in p.relative_to(root).parts for x in ('.venv','outputs','__pycache__','.pytest_cache'))}

def run(name, command, cwd, timeout=1200):
    started=time.time()
    with (ROOT/'logs'/f'{name}.log').open('w') as stream:
        stream.write('$ '+__import__('shlex').join(map(str,command))+'\n')
        stream.flush()
        try:
            proc=subprocess.run(list(map(str,command)),cwd=cwd,env=ENV,stdout=stream,stderr=subprocess.STDOUT,timeout=timeout)
            code=proc.returncode
        except subprocess.TimeoutExpired:
            code=124
            stream.write(f'\nValidation timeout after {timeout} seconds.\n')
    result=dict(name=name,command=list(map(str,command)),cwd=str(cwd),exit_code=code,
                elapsed_seconds=round(time.time()-started,3),log=str(ROOT/'logs'/f'{name}.log'),
                environment_overrides={k:ENV[k] for k in ENV if k not in os.environ or ENV[k]!=os.environ[k]})
    (ROOT/'results'/f'{name}.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result),flush=True)
    return result

def tests(name, source):
    before=source_hashes(source)
    result=run(name,[PYTHON,'-m','pytest','-q','-p','no:cacheprovider','--durations=12',
                     f'--junitxml={ROOT}/results/{name}.xml',f'--basetemp={ROOT}/tmp/{name}'],source)
    result['source_files_unchanged']=before==source_hashes(source)
    result['source_hashes']=before
    (ROOT/'results'/f'{name}.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(dict(name=name,source_files_unchanged=result['source_files_unchanged'])),flush=True)

def stage_fresh():
    stage=ROOT/'fresh-source'
    if not stage.exists():
        shutil.copytree(FRESH,stage,ignore=shutil.ignore_patterns('outputs','__pycache__','.pytest_cache','*.pyc','.venv'))
    return stage

def wheel():
    stage=stage_fresh()
    run('fresh-wheel',[PYTHON,'-m','pip','wheel','--no-deps','--no-build-isolation','--no-index',
                       '--wheel-dir',ROOT/'wheels',stage],ROOT)

def mesh():
    stage=stage_fresh()
    run('fresh-mesh',[PYTHON,stage/'build.py'],stage)

def cli():
    for name,path in [('fresh-cli',FRESH/'run.py'),('reproduction-cli',REPRODUCTION/'run.py'),
                      ('thesis-cli',REPRODUCTION/'run_thesis.py'),('printed-cli',REPRODUCTION/'run_printed.py')]:
        run(name,[PYTHON,path,'--help'],path.parent,timeout=120)

def wheel_smoke():
    wheel=next((ROOT/'wheels').glob('*.whl'))
    installed=run('wheel-install',[PYTHON,'-m','pip','install','--no-deps','--no-index',
                                  '--target',ROOT/'wheel-site',wheel],ROOT)
    if installed['exit_code']==0:
        run('wheel-smoke',[PYTHON,'-I',ROOT/'smoke_installed_wheel.py'],ROOT)

def regression():
    path=Path('/workspace/Khan/peridynamics/versions/workspace/bowen_reproduction')
    run('recovered-regression',[PYTHON,path/'verify_regression.py'],path,timeout=120)

if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('job',choices=['fresh-tests','reproduction-tests','recovered-fresh-tests',
                                     'recovered-reproduction-tests','wheel','mesh','cli','wheel-smoke','regression'])
    job=parser.parse_args().job
    if job=='fresh-tests': tests(job,FRESH)
    elif job=='reproduction-tests': tests(job,REPRODUCTION)
    elif job=='recovered-fresh-tests': tests(job,Path('/workspace/Khan/peridynamics/versions/workspace/bowen_from_scratch'))
    elif job=='recovered-reproduction-tests': tests(job,Path('/workspace/Khan/peridynamics/versions/workspace/bowen_reproduction'))
    elif job=='wheel': wheel()
    elif job=='mesh': mesh()
    elif job=='cli': cli()
    elif job=='wheel-smoke': wheel_smoke()
    elif job=='regression': regression()
