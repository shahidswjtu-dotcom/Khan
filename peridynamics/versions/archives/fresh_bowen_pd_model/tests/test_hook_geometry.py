import numpy as np
import pytest
from pdbeam.core import Benchmark
from pdbeam.assembly import Assembly,Reconstruction

def test_hook_topology_lengths_and_lumped_mass():
    p=Benchmark(); choices=Reconstruction.bulk_interface_hypothesis(p)
    m=Assembly(p,choices)
    delta=m.sx[m.rod_ends[:,1]]-m.sx[m.rod_ends[:,0]]
    vertical=np.isclose(abs(delta[:,2]),p.spacing)
    assert vertical.sum()==4*round(choices.hook_rise/p.spacing)
    assert m.rod_length.sum()==pytest.approx(2*(choices.bar_x_end-choices.bar_x_start+2*choices.hook_rise))
    assert m.steel_mass.sum()==pytest.approx(m.rod_length.sum()*m.area*p.steel_density)
    assert m.sx[:,2].max()==pytest.approx(.080)
    assert m.sx[:,2].min()==pytest.approx(.020)
