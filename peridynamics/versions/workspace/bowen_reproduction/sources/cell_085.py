def run_stage8(steps=5, stride=1, dt=p.dt, verbose=True):
    reset_stage8()
    uc=np.zeros_like(xyz);vc=np.zeros_like(xyz)
    us=np.zeros_like(steel_xyz);vs=np.zeros_like(steel_xyz)
    zh=float(hammer_z_initial);vh=float(HAMMER_V0)
    records=[]
    fc,fs,fh,contact_U,ncontact,pen,ncb,nib=source_shaped_forces(uc,us,zh,dt,commit=False)
    ac=fc/CONCRETE_NODE_MASS
    ast=fs/STEEL_MASS[:,None]
    ah=fh/HAMMER_MASS-(9.81 if USE_GRAVITY else 0.)
    mid=(np.abs(xyz[:,0]-p.L/2)<=dx/2+1e-12)&(xyz[:,2]>=p.H/2)
    for step in range(1,steps+1):
        vc+=0.5*dt*ac;vs+=0.5*dt*ast;vh+=0.5*dt*ah
        uc+=dt*vc;us+=dt*vs;zh+=dt*vh
        apply_supports(uc,vc)
        fc,fs,fh,contact_U,ncontact,pen,ncb,nib=source_shaped_forces(uc,us,zh,dt,commit=True)
        ac=fc/CONCRETE_NODE_MASS;ast=fs/STEEL_MASS[:,None]
        ah=fh/HAMMER_MASS-(9.81 if USE_GRAVITY else 0.)
        vc+=0.5*dt*ac;vs+=0.5*dt*ast;vh+=0.5*dt*ah
        apply_supports(uc,vc)
        if not (np.isfinite(uc).all() and np.isfinite(us).all() and np.isfinite(zh)):
            raise FloatingPointError(f'Nonfinite state at step {step}')
        maxu=np.max(np.abs(uc))
        if maxu>0.05:raise RuntimeError(f'Max concrete displacement {maxu:.3g} m at step {step}: aborting')
        if pen>0.5*dx:raise RuntimeError(f'Contact penetration {pen:.3g} m > 0.5 dx at step {step}')
        if step%stride==0 or step==1:
            rec=dict(step=step,time_ms=step*dt*1e3,
                     impact_force_kN=max(0.,fh)*1e-3,
                     midspan_deflection_mm=-np.mean(uc[mid,2])*1e3,
                     hammer_z_m=zh,hammer_velocity_mps=vh,
                     contact_count=ncontact,penetration_mm=pen*1e3,
                     concrete_broken=int(concrete8.broken.sum()),
                     interface_broken=int(interface8.broken.sum()),
                     kinetic_J=(0.5*CONCRETE_NODE_MASS*np.sum(vc*vc)
                         +0.5*np.sum(STEEL_MASS[:,None]*vs*vs)
                         +0.5*HAMMER_MASS*vh*vh),
                     contact_potential_J=contact_U)
            records.append(rec)
            if verbose:print(f"step={step:5d}, force={rec['impact_force_kN']:9.3f} kN, "
                             f"deflection={rec['midspan_deflection_mm']:.6f} mm, "
                             f"contacts={ncontact}, broken={rec['concrete_broken']}")
    return records,(uc,vc,us,vs,zh,vh)

RUN_STAGE8=False   # TURN ON only after reviewing the model limitations
STAGE8_STEPS=5     # 5-step test checks time stepping; it does not reach contact
STAGE8_STRIDE=1
if RUN_STAGE8:
    stage8_records,stage8_state=run_stage8(STAGE8_STEPS,STAGE8_STRIDE)
    from pathlib import Path
    import csv
    folder=Path('bowen_stage8_outputs');folder.mkdir(exist_ok=True)
    if stage8_records:
        with (folder/'stage8_diagnostics.csv').open('w',newline='') as f:
            w=csv.DictWriter(f,fieldnames=stage8_records[0].keys());w.writeheader();w.writerows(stage8_records)
else:
    print('Stage 8 full dynamics NOT executed by default. Set RUN_STAGE8=True for a controlled test.')
    print('Estimated first-contact step from Stage 6:',first_step)
