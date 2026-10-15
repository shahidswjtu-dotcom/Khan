"""Compare computed histories and damage with documented scalar targets."""
import argparse,csv,json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

root=Path(__file__).parent
parser=argparse.ArgumentParser()
parser.add_argument('--coarse',default='outputs/hooked_dt3e7')
parser.add_argument('--fine',default='outputs/hooked_dt1p5e7')
parser.add_argument('--report-dir',default='outputs/hooked_report')
args=parser.parse_args()
coarse=root/args.coarse; fine=root/args.fine
report=root/args.report_dir; report.mkdir(parents=True,exist_ok=True)

def history(folder):
    with (folder/'history.csv').open() as f:
        rows=list(csv.DictReader(f))
    return {k:np.array([float(row[k]) for row in rows]) for k in rows[0]}

h=history(coarse); hf=history(fine)
cm=json.loads((coarse/'manifest.json').read_text()); fm=json.loads((fine/'manifest.json').read_text())
cc=json.loads(json.dumps(cm)); ff=json.loads(json.dumps(fm))
cc['benchmark'].pop('dt'); ff['benchmark'].pop('dt')
same_configuration=cc==ff
comparison_key='timestep_comparison' if same_configuration else 'cross_variant_comparison'
coarse_label=('irreversible' if cm['reconstruction']['irreversible'] else 'printed history')+f", dt {cm['benchmark']['dt']*1e6:.2f} µs"
fine_label=('irreversible' if fm['reconstruction']['irreversible'] else 'printed history')+f", dt {fm['benchmark']['dt']*1e6:.2f} µs"
peak=int(np.argmax(h['force_N']))
summary={'time_window_ms':h['time_s'][-1]*1000,
         'peak_force_kN':h['force_N'][peak]/1000,
         'peak_force_time_ms':h['time_s'][peak]*1000,
         'first_nonzero_contact_time_ms':float(h['time_s'][np.flatnonzero(h['force_N']>0)[0]])*1000,
         'max_observed_bottom_midspan_mm':h['bottom_midspan_down_m'].max()*1000,
         'thesis_peak_force_kN':359.8,'thesis_max_deflection_mm':22.70,
         'experiment_peak_force_kN':301.,'experiment_max_deflection_mm':23.38,
         'max_momentum_residual_kg_m_s':h['momentum_residual_kg_m_s'].max(),
         'last_sampled_energy_residual_J':h['energy_balance_residual_J'][np.isfinite(h['energy_balance_residual_J'])][-1],
         'max_abs_sampled_energy_residual_J':np.nanmax(abs(h['energy_balance_residual_J'])),
         'max_abs_closed_system_material_energy_drift_J':float(np.nanmax(abs(h['energy_balance_residual_J']+h['contact_energy_change_J']+h['support_energy_change_J']))),
         'force_error_vs_thesis_percent':100*(h['force_N'][peak]/359800-1),
         'deflection_sensor_matched_to_thesis':False,'full_event_peak_established':False}
# Match physical time; never compare equal step numbers with different dt.
limit=hf['time_s'][-1]
sel=h['time_s']<=limit+1e-12
matched=h['force_N'][sel].max()/1000
summary[comparison_key]={'matched_time_ms':limit*1000,
    'coarse_peak_kN':matched,'fine_peak_kN':hf['force_N'].max()/1000,
    'peak_change_percent':100*(hf['force_N'].max()/1000/matched-1),
    'coarse_bottom_at_end_mm':h['bottom_midspan_down_m'][sel][-1]*1000,
    'fine_bottom_at_end_mm':hf['bottom_midspan_down_m'][-1]*1000,
    'coarse_impulse_N_s':float(np.sum(h['force_N'][sel])*3e-7),
    'fine_impulse_N_s':float(np.sum(hf['force_N'])*1.5e-7)}
if not same_configuration:
    summary[comparison_key]['note']='Different history/configuration and timestep; not a timestep convergence study.'
(report/'comparison.json').write_text(json.dumps(summary,indent=2)+'\n')
fig,axes=plt.subplots(2,2,figsize=(11,8),constrained_layout=True)
axes[0,0].plot(h['time_s']*1000,h['force_N']/1000,label=coarse_label)
axes[0,0].plot(hf['time_s']*1000,hf['force_N']/1000,label=fine_label,alpha=.7)
axes[0,0].axhline(359.8,color='k',ls='--',label='Bowen table peak')
axes[0,0].axhline(301.,color='gray',ls=':',label='Experiment table peak')
axes[0,0].set(xlabel='Time from declared initial state [ms]',ylabel='Impact force [kN]'); axes[0,0].legend(fontsize=8)
axes[0,1].plot(h['time_s']*1000,h['bottom_midspan_down_m']*1000)
axes[0,1].axhline(22.70,color='k',ls='--',label='Bowen table maximum; sensor unmatched')
axes[0,1].set(xlabel='Time [ms]',ylabel='Bottom midspan mean, downward [mm]'); axes[0,1].legend(fontsize=8)
valid=np.isfinite(h['energy_balance_residual_J'])
axes[1,0].plot(h['time_s'][valid]*1000,h['energy_balance_residual_J'][valid],label='Substep ledger residual')
axes[1,0].plot(h['time_s'][valid]*1000,(h['energy_balance_residual_J']+h['contact_energy_change_J']+h['support_energy_change_J'])[valid],label='Closed-system material drift')
axes[1,0].legend(fontsize=8)
axes[1,0].set(xlabel='Time [ms]',ylabel='Energy ledger residual [J]')
axes[1,1].semilogy(h['time_s']*1000,np.maximum(h['momentum_residual_kg_m_s'],1e-18))
axes[1,1].set(xlabel='Time [ms]',ylabel='Momentum residual [kg m/s]')
fig.suptitle('Fresh reconstruction with unresolved interface hypothesis')
fig.savefig(report/'comparison.png',dpi=180)
plt.close(fig)
mesh=np.load(root/'outputs/reference_mesh.npz',allow_pickle=False)
x=mesh['positions']; ends=mesh['ends']; vc=np.ones(len(ends))
# Eq. 2-10: uniform nodal volume cancels, including all reference family bonds.
denominator=np.zeros(len(x))
for side in (0,1): np.add.at(denominator,ends[:,side],vc)
files=sorted(coarse.glob('snapshot_*.npz'))
if files:
    damage_statistics=[]
    fig,axes=plt.subplots(len(files),1,figsize=(10,2.2*len(files)),constrained_layout=True,squeeze=False)
    for ax,path in zip(axes[:,0],files):
        with np.load(path,allow_pickle=False) as state:
            failed=~state['alive']; num=np.zeros(len(x))
            for side in (0,1): np.add.at(num,ends[:,side],vc*failed)
            damage=num/denominator
            damage_statistics.append({'step':int(state['step']),'time_ms':int(state['step'])*cm['benchmark']['dt']*1000,'global_mean_damage':float(damage.mean())})
            mid=abs(x[:,1]-.035)<.005+1e-12
            scat=ax.scatter(x[mid,0]*1000,x[mid,2]*1000,c=damage[mid],s=8,vmin=0,vmax=1,cmap='magma')
            ax.set(xlabel='Reference x [mm]',ylabel='Reference z [mm]',title=f"Family-volume damage (Eq. 2-10), {int(state['step'])*cm['benchmark']['dt']*1000:.2f} ms")
            ax.set_aspect('equal')
    fig.colorbar(scat,ax=axes[:,0],label='Eq. 2-10 damage, reference central slices')
    fig.savefig(report/'damage.png',dpi=180)
    (report/'damage_statistics.json').write_text(json.dumps(damage_statistics,indent=2)+'\n')
print(json.dumps(summary,indent=2))
