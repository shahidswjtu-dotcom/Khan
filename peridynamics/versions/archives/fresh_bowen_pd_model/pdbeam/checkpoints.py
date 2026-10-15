"""Atomic local checkpoints: readers see either old or complete new state."""
import os
from pathlib import Path
import tempfile
import numpy as np


def atomic_npz(path, **values):
    target=Path(path)
    if target.suffix!='.npz': target=target.with_name(target.name+'.npz')
    target.parent.mkdir(parents=True,exist_ok=True)
    descriptor,temporary=tempfile.mkstemp(prefix='.'+target.name+'.',suffix='.tmp',dir=target.parent)
    try:
        with os.fdopen(descriptor,'wb') as stream:
            np.savez_compressed(stream,**values)
            stream.flush()
            os.fsync(stream.fileno())
        os.replace(temporary,target)
    finally:
        if os.path.exists(temporary): os.unlink(temporary)
