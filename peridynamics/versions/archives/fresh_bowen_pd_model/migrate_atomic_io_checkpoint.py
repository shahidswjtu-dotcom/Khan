"""Explicit checkpoint migration for the verified atomic-I/O-only patch.

Reject every other source/parameter change. Never overwrite the original.
"""
import argparse,hashlib,json
from pathlib import Path
import numpy as np
from pdbeam.core import Benchmark
from pdbeam.assembly import Assembly,Reconstruction
from pdbeam.checkpoints import atomic_npz

root=Path(__file__).parent
parser=argparse.ArgumentParser()
parser.add_argument('checkpoint'); parser.add_argument('destination')
a=parser.parse_args(); source=Path(a.checkpoint); destination=Path(a.destination)
if source.resolve()==destination.resolve(): raise ValueError('Preserve original checkpoint')
old=root/'versions/fresh1/pdbeam'; new=root/'pdbeam'
for f in old.glob('*.py'):
    expected=f.read_text()
    if f.name in ('core.py','solver.py'):
        expected=expected.replace('import numpy as np','import numpy as np\nfrom .checkpoints import atomic_npz',1).replace('np.savez_compressed(path,','atomic_npz(path,')
    if expected!=(new/f.name).read_text(): raise ValueError(f'Non-I/O source change: {f.name}')
if {f.name for f in new.glob('*.py')}!={f.name for f in old.glob('*.py')}|{'checkpoints.py'}:
    raise ValueError('Unexpected package module change')
with np.load(source,allow_pickle=False) as d: payload={k:d[k].copy() for k in d.files}
old_manifest=json.loads(str(payload['manifest']))
for name,sha in old_manifest['source_sha256'].items():
    if hashlib.sha256((old/name).read_bytes()).hexdigest()!=sha: raise ValueError('Archived source does not match checkpoint')
p=Benchmark(**old_manifest['benchmark'])
choices=dict(old_manifest['reconstruction']); choices['bar_y']=tuple(choices['bar_y'])
m=Assembly(p,Reconstruction(**choices))
new_manifest=json.loads(m.manifest)
if {k:v for k,v in old_manifest.items() if k!='source_sha256'}!={k:v for k,v in new_manifest.items() if k!='source_sha256'}:
    raise ValueError('Configuration/protocol changed')
payload['manifest']=np.array(m.manifest)
atomic_npz(destination,**payload)
record=dict(source=str(source),destination=str(destination),step=int(payload['step']),
            original_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
            migrated_sha256=hashlib.sha256(destination.read_bytes()).hexdigest(),
            correction='Atomic checkpoint writer only; deterministic source transformation verified; all state arrays and scalars unchanged')
destination.with_suffix('.migration.json').write_text(json.dumps(record,indent=2)+'\n')
print(json.dumps(record,indent=2))
