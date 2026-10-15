import numpy as np
import pytest
from pdbeam.checkpoints import atomic_npz

def test_atomic_failure_preserves_previous_checkpoint(tmp_path,monkeypatch):
    path=tmp_path/'checkpoint.npz'
    atomic_npz(path,step=np.int64(5))
    def fail(*args,**kwargs):
        args[0].write(b'partial')
        raise OSError('Synthetic disk-write failure')
    monkeypatch.setattr(np,'savez_compressed',fail)
    with pytest.raises(OSError): atomic_npz(path,step=np.int64(6))
    with np.load(path,allow_pickle=False) as d: assert int(d['step'])==5
    assert list(tmp_path.iterdir())==[path]

def test_atomic_success_replaces_complete_state(tmp_path):
    path=tmp_path/'checkpoint.npz'
    atomic_npz(path,step=np.int64(5))
    atomic_npz(path,step=np.int64(6),u=np.ones((3,3)))
    with np.load(path,allow_pickle=False) as d:
        assert int(d['step'])==6
        np.testing.assert_array_equal(d['u'],np.ones((3,3)))
