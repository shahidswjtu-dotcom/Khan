import zipfile
"""Extend the response window without changing physics or tuning to targets."""
from pathlib import Path
import subprocess,sys,time
import numpy as np
root=Path(__file__).parent
checkpoint=root/'outputs/hooked_dt3e7/checkpoint.npz'
while True:
    try:
        with np.load(checkpoint,allow_pickle=False) as state: step=int(state['step'])
        if step>=5000: break
    except (OSError,ValueError,EOFError,zipfile.BadZipFile): pass
    time.sleep(2)
subprocess.run([sys.executable,str(root/'run.py'),'--steps','10000','--restart',str(checkpoint),'--out',str(root/'outputs/hooked_continuation_3ms'),'--accept-bulk-interface-hypothesis'],check=True,cwd=root)
