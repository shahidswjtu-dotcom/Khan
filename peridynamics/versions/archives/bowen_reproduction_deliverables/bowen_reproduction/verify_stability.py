"""Conservative undeformed elastic spectral upper bounds, not nonlinear proof."""
import json
import numpy as np
from bowen_pd.model import create
m=create();p=m.p
k=72*p.E_concrete/(np.pi*m.delta**4)*(1-(m.LEN/m.delta)**2)**2*m.V**2*m.PAIR_WEIGHT/m.LEN
rows=np.bincount(m.I,weights=k,minlength=len(m.xyz))+np.bincount(m.J,weights=k,minlength=len(m.xyz))
# rank-one axial bond block norm k; graph bound lambda_max <= 2 max(sum k/m)
concrete_bound=2*max(rows)/m.CONCRETE_NODE_MASS
ks=p.E_steel*m.steel_A/m.steel_lengths
rs=np.bincount(m.rod_pairs[:,0],weights=ks,minlength=len(m.steel_xyz))+np.bincount(m.rod_pairs[:,1],weights=ks,minlength=len(m.steel_xyz))
rod_bound=2*max(rs/m.STEEL_MASS)
it=m.interface;pr=m.params
ki=pr.bcc*(1-(it.lengths/pr.horizon)**2)**2*pr.quadrature_volume*(m.STEEL_MASS[it.steel_ids]/p.rho_steel)/it.lengths
# sum rank-one mass-normalized interaction norms gives a valid loose bound
interface_bound=float(np.sum(ki*(1/m.STEEL_MASS[it.steel_ids]+np.sum(it.weights**2,axis=1)/m.CONCRETE_NODE_MASS)))
# every point contact is deliberately conservative; relevant contact family smaller
contact_bound=m.CONTACT_STIFFNESS*(1/m.CONCRETE_NODE_MASS+len(m.xyz)/m.HAMMER_MASS)
r={'concrete_lambda_upper_s-2':float(concrete_bound),'rod_lambda_upper_s-2':float(rod_bound),'interface_lambda_upper_s-2_loose':interface_bound,'contact_lambda_upper_s-2_all_particles':float(contact_bound),'concrete_alone_safe_dt_s':float(2/np.sqrt(concrete_bound)),'rod_alone_safe_dt_s':float(2/np.sqrt(rod_bound)),'combined_verlet_sufficient_dt_s':float(2/np.sqrt(concrete_bound+rod_bound+interface_bound+contact_bound)),'qualification':'Sufficient bound for undeformed elastic Verlet only. Rate enhancement, compression geometric stiffness, projection, and damage require separate analysis; exceeding a sufficient bound does not prove instability.'}
print(json.dumps(r,indent=2))
from pathlib import Path
Path('outputs/stability.json').write_text(json.dumps(r,indent=2))
