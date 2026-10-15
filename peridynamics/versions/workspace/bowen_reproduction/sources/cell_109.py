from time import perf_counter
import csv
from pathlib import Path
OUT14=Path('bowen_stage14_outputs');OUT14.mkdir(exist_ok=True)

def source_c_force(uc,us,mode,commit=False):
    result=appendix_c_beam.evaluate(uc,us,mode=mode,commit=commit)
    return result['concrete_force'], result['steel_force'],int((result['failfail']==0).sum())

def run_stage14(mode='literal',steps=380,dt=p.dt,stride=5,limit_mm=50.):
    if mode not in ('literal','irreversible_corrected'):raise ValueError(mode)
    appendix_c_beam.reset(); fast_concrete.reset()
    uc=np.zeros_like(xyz);vc=np.zeros_like(uc)
    us=np.zeros_like(steel_xyz);vs=np.zeros_like(us)
    zh=float(hammer_z_initial);vh=float(HAMMER_V0)
    mid=(np.abs(xyz[:,0]-p.L/2)<=dx/2+1.e-12)&(xyz[:,2]>=p.H/2)
    def forces(commit):
        fc,n_b=fast_concrete.evaluate(uc,dt,commit)
        fi,fs,n_i=source_c_force(uc,us,mode,commit)
        fc+=fi
        fs+=rod_forces(us)
        ct,fh,cu,nc,pen=sphere_contact(xyz+uc,zh)
        fc+=ct
        return fc,fs,float(fh[2]),int(nc),float(pen),int(n_b),int(n_i)
    fc,fs,fh,nc,pen,nb,ni=forces(False)
    ac=fc/CONCRETE_NODE_MASS;ast=fs/STEEL_MASS[:,None]
    ah=fh/HAMMER_MASS-(9.81 if USE_GRAVITY else 0.)
    rows=[];first_contact=None;stop=None;t0=perf_counter()
    for step in range(1,steps+1):
        vc+=0.5*dt*ac;vs+=0.5*dt*ast;vh+=0.5*dt*ah
        uc+=dt*vc;us+=dt*vs;zh+=dt*vh
        apply_supports(uc,vc)
        fc,fs,fh,nc,pen,nb,ni=forces(True)
        ac=fc/CONCRETE_NODE_MASS;ast=fs/STEEL_MASS[:,None]
        ah=fh/HAMMER_MASS-(9.81 if USE_GRAVITY else 0.)
        vc+=0.5*dt*ac;vs+=0.5*dt*ast;vh+=0.5*dt*ah
        apply_supports(uc,vc)
        if nc and first_contact is None:first_contact=step
        max_c=float(np.max(np.linalg.norm(uc,axis=1)))*1e3
        max_s=float(np.max(np.linalg.norm(us,axis=1)))*1e3
        if (not np.isfinite(uc).all()) or (not np.isfinite(us).all()) or max(max_c,max_s)>limit_mm or pen>0.5*dx:
            stop=f'safety stop at step {step}: max_c={max_c:.4g} mm max_s={max_s:.4g} mm'
        if step==1 or step%stride==0 or nc and step==first_contact or stop:
            rows.append(dict(step=step,time_ms=step*dt*1000,force_kN=max(0.,fh)*1e-3,
                contact_points=nc,midspan_mm=-1000*np.mean(uc[mid,2]),
                max_concrete_mm=max_c,max_steel_mm=max_s,broken_concrete=nb,
                broken_interface=ni,hammer_v_mps=vh,penetration_mm=pen*1000))
        if stop:break
    output=OUT14/f'{mode}_history.csv'
    with output.open('w',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
    summary=dict(mode=mode,steps_executed=step,first_contact_step=first_contact,
                 elapsed_s=perf_counter()-t0,stop=stop,
                 peak_recorded_kN=max(x['force_kN'] for x in rows))
    print('Stage 14:',summary)
    return rows,summary

RUN_STAGE14=False
if RUN_STAGE14:
    # Run separately: each call resets concrete and interface histories.
    result_literal,summary_literal=run_stage14('literal',steps=380)
    result_corrected,summary_corrected=run_stage14('irreversible_corrected',steps=380)
else:
    print('Stage 14 dynamic run is OFF. Set RUN_STAGE14=True to test the unverified benchmark inputs.')
