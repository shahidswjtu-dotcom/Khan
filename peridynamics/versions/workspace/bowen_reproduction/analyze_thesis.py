"""Compare the integrated reconstruction without inventing reference curves."""
import csv,json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
p=Path(__file__).parent/'outputs/thesis_integrated'
rows=list(csv.DictReader((p/'history.csv').open()))
continuation=p.parent/'thesis_complete5000'
if (continuation/'summary.json').exists():
    rows+=list(csv.DictReader((continuation/'history.csv').open()))
    p=continuation
arrays={key:np.array([float(r[key]) for r in rows]) for key in rows[0]}
t=arrays['time_ms'];F=arrays['impact_kN'];u=arrays['support_relative_bottom_mm']
result={'observed_peak_force_kN':float(F.max()),'observed_peak_time_ms':float(t[F.argmax()]),'max_observed_bottom_probe_mm':float(u.max()),'force_difference_vs_Bowen_percent':float(100*(F.max()/359.8-1)),'force_difference_vs_experiment_percent':float(100*(F.max()/301-1)),'deflection_difference_vs_Bowen_percent':float(100*(u.max()/22.70-1)),'deflection_difference_vs_experiment_percent':float(100*(u.max()/23.38-1)),'qualification':'Observed over computed interval; sensor unmatched; no full-response or exact author equivalence claim'}
if (p/'damage.npz').exists():
 from thesis_pd.geometry import rod_pairs,steel_lengths,steel_E,p as material
 with np.load(p/'damage.npz') as z:
  x=z['xyz'];damage=z['damage'];xs=z['steel_xyz'];us=z['us']
  unit=(xs[rod_pairs[:,1]]-xs[rod_pairs[:,0]])/steel_lengths[:,None]
  strain=np.sum((us[rod_pairs[:,1]]-us[rod_pairs[:,0]])*unit,axis=1)/steel_lengths
  result['max_abs_rod_strain']=float(max(abs(strain)));result['yield_strain']=material.fy_steel/steel_E;result['elastic_rods_exceed_yield']=bool(max(abs(strain))>material.fy_steel/steel_E)
  mask=abs(x[:,1]-.035)<=.0025+1e-12
  fig,ax=plt.subplots(figsize=(10,3),layout='constrained');sc=ax.scatter(x[mask,0]*1000,x[mask,2]*1000,c=damage[mask],s=8,vmin=0,vmax=1,cmap='inferno');ax.set_aspect('equal');ax.set_xlabel('Reference x (mm)');ax.set_ylabel('Reference z (mm)');ax.set_title('Integrated thesis reconstruction: computed damage, not validated reproduction');fig.colorbar(sc,ax=ax,label='Unweighted failed-family fraction');fig.savefig(p/'damage_reference.png',dpi=200);plt.close(fig)
fig,axes=plt.subplots(3,1,figsize=(8,9),sharex=True,layout='constrained')
axes[0].plot(t,F,label='Computed reaction');axes[0].axhline(359.8,color='k',ls='--',label='Bowen359.8kN (peak reference only)');axes[0].axhline(301,color='grey',ls=':',label='Experimental301kN (peak reference only)');axes[0].set_ylabel('Reaction (kN)');axes[0].legend(fontsize=8)
axes[1].plot(t,u);axes[1].set_ylabel('Bottom layer relative\ndeflection (mm)')
axes[2].plot(t,arrays['envelope_energy_budget_residual_J'],label='Envelope energy residual');axes[2].plot(t,arrays['repulsion_potential_J'],label='Short-range potential');axes[2].set_ylabel('Energy (J)');axes[2].set_xlabel('Time from initial state (ms)');axes[2].legend()
fig.suptitle('Source-based reconstruction with declared unresolved inputs');fig.savefig(p/'comparison.png',dpi=200);fig.savefig(p/'comparison.pdf');plt.close(fig)
(p/'comparison.json').write_text(json.dumps(result,indent=2));print(json.dumps(result,indent=2))
