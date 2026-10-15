from pathlib import Path
import csv, json
from datetime import datetime

STAGE6_DIR=Path('bowen_stage6_outputs')
STAGE6_DIR.mkdir(exist_ok=True)
TARGET_PEAK_FORCE_KN=359.8
TARGET_PEAK_DEFLECTION_MM=22.70
REPRODUCTION_AUDIT={
    'beam_geometry': ('verified_in_notebook','600 x 70 x 100 mm; 33,600 material points'),
    'bond_force_assembly': ('verified_unit_tests','pairwise momentum conservation'),
    'microelastic_fracture': ('partly_verified','continuous radial fracture energy checked, discrete fracture energy pending'),
    'strain_rate_DIF': ('implemented_not_source_equivalent','DIF_E force-law scaling and update order require original code'),
    'PDROD': ('prototype_only','2 samples per steel node; linear penalty, not thesis nonlinear interface'),
    'contact': ('prototype_only','per-particle spherical penalty; Bowen contact resolution not established'),
    'supports': ('provisional','two point-row vertical masks; fixture agreement not established'),
    'hammer_motion': ('implemented','free rigid hammer with v0 = -10.8 m/s'),
    'energy': ('incomplete','no validated elastic+fracture+rate dissipation balance'),
    'experimental_comparison': ('not_attempted','no full benchmark computation verified')
}
with (STAGE6_DIR/'source_fidelity_audit.csv').open('w',newline='') as f:
    w=csv.writer(f);w.writerow(['component','status','evidence_or_missing_detail'])
    for k,(status,desc) in REPRODUCTION_AUDIT.items():w.writerow([k,status,desc])
print('Audit saved to',STAGE6_DIR/'source_fidelity_audit.csv')
for k,(status,_) in REPRODUCTION_AUDIT.items():print(f'{k:24s}: {status}')
