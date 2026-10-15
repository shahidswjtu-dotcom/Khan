from pathlib import Path
import csv, json, time
import numpy as np
STAGE17_DIR=Path('bowen_stage17_outputs')
STAGE17_DIR.mkdir(exist_ok=True)

def normalized_acceleration(physical_concrete, physical_steel, interface_concrete_density,
                            interface_steel_density):
    """Only the Appendix C entries are force density; all other inputs are nodal N."""
    return (physical_concrete/CONCRETE_NODE_MASS + interface_concrete_density/p.rho_concrete,
            physical_steel/STEEL_MASS[:,None] + interface_steel_density/p.rho_steel)

# Check the normalization independently of contact and fracture.
fq=np.array([[74204.,0.,0.]])
a_source=fq/p.rho_steel
m=STEEL_MASS[0]
a_wrong=fq/m
assert np.linalg.norm(a_wrong)>np.linalg.norm(a_source)
np.testing.assert_allclose(a_source, fq*STEEL_MASS[0]/p.rho_steel/STEEL_MASS[0])
print('PASS: source acceleration = force-density / steel density')
print('Example source acceleration:',a_source[0,0], 'm/s^2, previous lumped mass:',a_wrong[0,0])