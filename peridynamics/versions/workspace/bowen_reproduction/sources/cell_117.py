from numba import njit
from pathlib import Path
import csv, json, time
STAGE18_DIR=Path('bowen_stage18_outputs')
STAGE18_DIR.mkdir(exist_ok=True)
STEEL_CONTROL_VOLUME = STEEL_MASS / p.rho_steel  # m^3, derived from rod lumped masses
assert np.all(STEEL_CONTROL_VOLUME>0)
print('Steel effective node volume range m^3:',float(STEEL_CONTROL_VOLUME.min()),float(STEEL_CONTROL_VOLUME.max()))
print('Concrete node volume m^3:',float(V))