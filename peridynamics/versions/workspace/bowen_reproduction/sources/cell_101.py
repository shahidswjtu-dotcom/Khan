from pathlib import Path
import csv

STAGE11_DIR = Path('bowen_stage11_outputs')
STAGE11_DIR.mkdir(exist_ok=True)

def diagnostic_run(dt=p.dt, t_end=0.000115, write_prefix='baseline', 
                   max_disp=0.05, max_pen=0.5*dx, verbose=True):
    fast_concrete.reset(); fast_interface.reset()
    uc=np.zeros_like(xyz); vc=np.zeros_like(uc)
    us=np.zeros_like(steel_xyz); vs=np.zeros_like(us)
    zh=float(hammer_z_initial); vh=float(HAMMER_V0)
    fc,fs,fh,cu,nc,pen,nb,ni=fast_system_forces(uc,us,zh,dt,False)
    ac=fc/CONCRETE_NODE_MASS; ast=fs/STEEL_MASS[:,None]
    ah=fh/HAMMER_MASS-(9.81 if USE_GRAVITY else 0.)
    mid=(np.abs(xyz[:,0]-p.L/2)<=dx/2+1.e-12)&(xyz[:,2]>=p.H/2)
    records=[]; failure=None; snapshot=None
    for step in range(1,int(np.ceil(t_end/dt))+1):
        vc+=0.5*dt*ac; vs+=0.5*dt*ast; vh+=0.5*dt*ah
        uc+=dt*vc; us+=dt*vs; zh+=dt*vh
        apply_supports(uc,vc)
        fc,fs,fh,cu,nc,pen,nb,ni=fast_system_forces(uc,us,zh,dt,True)
        ac=fc/CONCRETE_NODE_MASS; ast=fs/STEEL_MASS[:,None]
        ah=fh/HAMMER_MASS-(9.81 if USE_GRAVITY else 0.)
        vc+=0.5*dt*ac; vs+=0.5*dt*ast; vh+=0.5*dt*ah
        apply_supports(uc,vc)
        maxdisp=np.linalg.norm(uc,axis=1)
        maxidx=int(np.argmax(maxdisp)); maxacc=np.linalg.norm(ac,axis=1)
        maxai=int(np.argmax(maxacc));maxsi=int(np.argmax(np.linalg.norm(us,axis=1)))
        kinetic_con=0.5*CONCRETE_NODE_MASS*np.sum(vc*vc)
        kinetic_steel=0.5*np.sum(STEEL_MASS[:,None]*vs*vs)
        kinetic_ham=0.5*HAMMER_MASS*vh*vh
        q=dict(step=step,time_ms=step*dt*1e3,contact_kN=max(0,fh)*1e-3,
               contact_nodes=int(nc),penetration_mm=pen*1e3,
               midspan_mm=-float(np.mean(uc[mid,2]))*1e3,
               max_disp_mm=float(maxdisp[maxidx])*1e3,max_disp_node=maxidx,
               max_disp_x_mm=xyz[maxidx,0]*1e3,max_disp_y_mm=xyz[maxidx,1]*1e3,
               max_disp_z_mm=xyz[maxidx,2]*1e3,
               max_acc_mps2=float(maxacc[maxai]),max_acc_node=maxai,
               max_steel_disp_mm=float(np.linalg.norm(us[maxsi]))*1e3,
               max_steel_node=maxsi,concrete_broken=int(fast_concrete.broken.sum()),
               interface_broken=int(fast_interface.broken.sum()),
               concrete_kinetic_J=float(kinetic_con),steel_kinetic_J=float(kinetic_steel),
               hammer_kinetic_J=float(kinetic_ham),contact_potential_J=float(cu),
               hammer_velocity_mps=float(vh))
        records.append(q)
        finite=(np.isfinite(uc).all() and np.isfinite(us).all() and np.isfinite(zh))
        if not finite or maxdisp[maxidx]>max_disp or pen>max_pen:
            failure=f'At step {step}: max displacement {maxdisp[maxidx]*1e3:.3f} mm, penetration {pen*1e3:.5f} mm'
            break
        if verbose and (step%50==0 or (nc>0 and step>1 and records[-2]['contact_nodes']==0)):
            print(f'step={step} time={q["time_ms"]:.5f}ms F={q["contact_kN"]:.3f}kN maxU={q["max_disp_mm"]:.4f}mm maxA={q["max_acc_mps2"]:.3g}')
    with open(STAGE11_DIR/f'{write_prefix}_history.csv','w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(records[0]));w.writeheader();w.writerows(records)
    np.savez_compressed(STAGE11_DIR/f'{write_prefix}_last_state.npz',
        uc=uc,vc=vc,us=us,vs=vs,zh=zh,vh=vh,
        broken=fast_concrete.broken,interface_broken=fast_interface.broken)
    print('RUN:',write_prefix,'steps',len(records),'failure:',failure or 'none')
    return records,failure

RUN_STAGE11 = False
if RUN_STAGE11:
    stage11_baseline, baseline_failure=diagnostic_run(write_prefix='dt_3e-7')
    stage11_half, half_failure=diagnostic_run(dt=p.dt/2,
       write_prefix='dt_1p5e-7')
    import matplotlib.pyplot as plt
    fig,axs=plt.subplots(3,1,figsize=(9,10),sharex=True)
    for rr,label in [(stage11_baseline,'dt=3e-7 s'),(stage11_half,'dt=1.5e-7 s')]:
        tt=np.array([x['time_ms'] for x in rr]);
        axs[0].plot(tt,[x['contact_kN'] for x in rr],label=label)
        axs[1].plot(tt,[x['max_disp_mm'] for x in rr],label=label)
        axs[2].plot(tt,[x['max_acc_mps2'] for x in rr],label=label)
    for ax,lbl in zip(axs,['Contact force (kN)','Maximum concrete point displacement (mm)','Max concrete acceleration (m/s²)']):
        ax.set_ylabel(lbl);ax.grid(alpha=.3);ax.legend()
    axs[-1].set_xlabel('Physical time (ms)');fig.tight_layout()
    fig.savefig(STAGE11_DIR/'timestep_comparison.png',dpi=180)
    print('Outputs saved to',STAGE11_DIR)
else:
    print('Stage 11 ready. Set RUN_STAGE11=True for both diagnostic runs.')
