"""Comparison of diagnostic trajectories; no fitting or digitized curves invented."""
from pathlib import Path
import csv,json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
r=Path(__file__).parent;rows=[];fig,axs=plt.subplots(3,1,figsize=(8,9),sharex=True,layout='constrained')
for name in ['legacy550','legacy_extended','legacy_half_dt','projection_attempt','family_projection']:
 folder=r/'outputs'/name
 if not (folder/'summary.json').exists():continue
 data=[]
 for p in [folder/'history.csv',folder/'history_resume.csv']:
  if p.exists():data+=list(csv.DictReader(p.open()))
 if name=='legacy_extended':data=list(csv.DictReader((r/'outputs/legacy550/history.csv').open()))+data
 if not data:continue
 t=np.array([float(x['time_ms']) for x in data]);force=np.array([float(x['impact_kN']) for x in data]);u=np.array([float(x['support_relative_bottom_mm']) for x in data])
 summary=json.loads((folder/'summary.json').read_text());peak=force.max();deflection=u.max()
 rows.append(dict(run=name,completed_steps=summary['steps_completed'],time_ms=float(t[-1]),peak_force_kN=float(peak),peak_force_time_ms=float(t[force.argmax()]),bottom_deflection_mm=float(deflection),force_error_vs_Bowen_percent=100*(peak/359.8-1),force_error_vs_experiment_percent=100*(peak/301-1),deflection_error_vs_Bowen_percent=100*(deflection/22.70-1),deflection_error_vs_experiment_percent=100*(deflection/23.38-1),status=summary['stop'] or 'requested interval completed; full response not established'))
 axs[0].plot(t,force,label=name);axs[1].plot(t,u,label=name);axs[2].plot(t,[float(x['numerical_work_budget_residual_J']) for x in data],label=name)
axs[0].axhline(359.8,c='black',ls='--',lw=.7,label='Bowen peak (timing unavailable)');axs[0].axhline(301,c='grey',ls=':',lw=.7,label='Experimental peak (timing unavailable)')
axs[0].set_ylabel('Contact reaction (kN)');axs[1].set_ylabel('Bottom centre-layer\ndeflection (mm)');axs[2].set_ylabel('Work budget residual (J)');axs[2].set_xlabel('Time from initialized state (ms)');axs[0].legend(fontsize=7)
fig.suptitle('Diagnostic reconstructions: incomplete response; coupling calibration unverified',fontsize=10)
fig.savefig(r/'outputs/comparison.png',dpi=200);fig.savefig(r/'outputs/comparison.pdf');plt.close(fig)
with (r/'outputs/comparison.csv').open('w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
# Damage in REFERENCE coordinates; don't mistake detached particle motion for crack openings.
for name in ['legacy_extended','legacy_half_dt','projection_attempt','family_projection']:
 p=r/'outputs'/name/'damage.npz'
 if not p.exists():continue
 with np.load(p) as d:
  x=d['xyz'];mask=abs(x[:,1]-.035)<=.0025+1e-12
  fig,ax=plt.subplots(figsize=(10,2.6),layout='constrained')
  q=ax.scatter(x[mask,0]*1000,x[mask,2]*1000,c=d['damage'][mask],s=8,vmin=0,vmax=1,cmap='inferno')
  ax.set_aspect('equal');ax.set_xlabel('Reference x (mm)');ax.set_ylabel('Reference z (mm)');ax.set_title(name+': unweighted failed-family fraction, final diagnostic state');fig.colorbar(q,ax=ax,label='Damage')
  fig.savefig(p.parent/'damage_reference.png',dpi=200);plt.close(fig)
print(json.dumps(rows,indent=2))
