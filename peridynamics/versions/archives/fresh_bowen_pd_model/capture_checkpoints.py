import zipfile
"""Copy completed milestones without changing a running simulation."""
from pathlib import Path
import argparse,time
from pdbeam.checkpoints import atomic_npz
import numpy as np
parser=argparse.ArgumentParser()
parser.add_argument('--folder',default='outputs/hooked_dt3e7')
a=parser.parse_args()
root=Path(__file__).parent/a.folder
seen=set()
while True:
    path=root/'checkpoint.npz'
    if path.exists():
        try:
            with np.load(path,allow_pickle=False) as data:
                payload={k:data[k].copy() for k in data.files}
                step=int(payload['step'])
            if step not in seen:
                atomic_npz(root/f'snapshot_{step:05d}.npz',**payload)
                seen.add(step)
                print(step,flush=True)
            if step>=5000: break
        except (OSError,ValueError,EOFError,zipfile.BadZipFile):
            pass
    time.sleep(2)
