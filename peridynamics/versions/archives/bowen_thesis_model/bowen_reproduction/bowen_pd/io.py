"""Atomic checkpoints, checksums, strict continuation compatibility."""
import json,hashlib,os
from pathlib import Path
import numpy as np

def source_hash():
    h=hashlib.sha256()
    for p in sorted(Path(__file__).parent.glob('*.py')):
        h.update(p.name.encode());h.update(p.read_bytes())
    return h.hexdigest()

def save(path,arrays,metadata):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    temp=path.with_suffix('.tmp')
    with temp.open('wb') as f:
        np.savez(f,metadata=np.array(json.dumps(metadata)),**arrays)
        f.flush();os.fsync(f.fileno())
    os.replace(temp,path)

def load(path,expected):
    with np.load(path,allow_pickle=False) as z:
        meta=json.loads(str(z['metadata']))
        for key,val in expected.items():
            if meta.get(key)!=val: raise ValueError('Checkpoint mismatch: '+key)
        return {k:z[k].copy() for k in z.files if k!='metadata'},meta
