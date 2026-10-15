"""Regression verifies the original saved sample times, not a newly fitted peak."""
from pathlib import Path
import csv,json
import numpy as np
r=Path(__file__).parent
base=list(csv.DictReader((r/'baselines/stage19/extended_impact_history.csv').open()))
run={int(x['step']):x for x in csv.DictReader((r/'outputs/legacy550/history.csv').open())}
for old,new in [('impact_kN','impact_kN'),('midspan_mm','legacy_half_depth_mm'),('max_concrete_mm','max_concrete_mm'),('max_steel_mm','max_steel_mm'),('broken_concrete','broken_concrete')]:
 np.testing.assert_allclose([float(x[old]) for x in base],[float(run[int(x['step'])][new]) for x in base],rtol=1e-8,atol=1e-7)
print('Stage19 sample regression passed:',len(base),'samples')
