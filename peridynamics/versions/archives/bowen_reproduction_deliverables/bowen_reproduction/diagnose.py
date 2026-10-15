"""Diagnose detached particles and compression overlap at saved states."""
import json
from pathlib import Path
import numpy as np
from scipy.spatial import cKDTree
root=Path(__file__).parent
for folder in ['legacy_extended','projection_attempt','legacy_half_dt','family_projection']:
 p=root/'outputs'/folder
 if not (p/'final.npz').exists():continue
 with np.load(p/'final.npz') as s,np.load(p/'damage.npz') as d:
  X=d['xyz'];u=s['u'];v=s['v'];Y=X+u;k=int(np.argmax(np.linalg.norm(u,axis=1)))
  # Short-range criteria from Eq2-49, all current neighbours, irrespective of old rupture.
  pair=cKDTree(Y).query_pairs(.00675,output_type='ndarray')
  r0=np.linalg.norm(X[pair[:,1]]-X[pair[:,0]],axis=1);length=np.linalg.norm(Y[pair[:,1]]-Y[pair[:,0]],axis=1)
  threshold=np.minimum(.9*r0,.00675);mask=length<threshold
  result={'diagnostic':'Overlap threshold Eq2-49; pair counts alone are not force/energy or proof of instability',
   'largest_displacement_node':k,'reference_m':X[k].tolist(),'displacement_m':u[k].tolist(),'velocity_mps':v[k].tolist(),'node_damage':float(d['damage'][k]),'largest_speed_mps':float(np.linalg.norm(v,axis=1).max()),'short_range_active_pair_count':int(mask.sum()),'minimum_current_pair_distance_m':float(length.min(initial=np.inf))}
  (p/'diagnosis.json').write_text(json.dumps(result,indent=2));print(folder,json.dumps(result))
