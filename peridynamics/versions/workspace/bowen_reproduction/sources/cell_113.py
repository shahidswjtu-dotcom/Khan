def stage17_components(uc,us,zh,dt,commit=False,mode='irreversible_corrected'):
    """Separate physical nodal force contributions from source-shaped interface density entries."""
    fcp,broken_con=fast_concrete.evaluate(uc,dt,commit)
    fsp=rod_forces(us)
    fi_c,fi_s,ni=source_c_force(uc,us,mode,commit)
    contact_c,hammer_f,contact_U,ncontact,maxpen=sphere_contact(xyz+uc,zh)
    fcp+=contact_c
    return fcp,fsp,fi_c,fi_s,float(hammer_f[2]),int(ncontact),float(maxpen),int(broken_con),int(ni)

def stage17_run(steps=420, dt=p.dt, save_every=2,mode='irreversible_corrected',limit_mm=50.,verbose=True):
    assert mode in ('literal','irreversible_corrected')
    appendix_c_beam.reset();fast_concrete.reset()
    uc=np.zeros_like(xyz);vc=np.zeros_like(xyz)
    us=np.zeros_like(steel_xyz);vs=np.zeros_like(steel_xyz)
    zh=float(hammer_z_initial);vh=float(HAMMER_V0)
    mid=(np.abs(xyz[:,0]-p.L/2)<=dx/2+1e-12)&(xyz[:,2]>=p.H/2)
    def force_step(commit):
        physical_c,physical_s,density_c,density_s,fh,nc,pen,nb,ni=stage17_components(uc,us,zh,dt,commit,mode)
        ac,ast=normalized_acceleration(physical_c,physical_s,density_c,density_s)
        ah=fh/HAMMER_MASS-(9.81 if USE_GRAVITY else 0.0)
        physical_interface_sum=density_c.sum(axis=0)*dx**3 + np.sum(density_s*STEEL_MASS[:,None]/p.rho_steel,axis=0)
        return ac,ast,ah,fh,nc,pen,nb,ni,float(np.linalg.norm(physical_interface_sum)),float(np.max(np.linalg.norm(ast,axis=1)))
    ac,ast,ah,fh,nc,pen,nb,ni,mom_error,amax_s=force_step(False)
    rows=[];error=None;first_contact=None;t0=time.perf_counter()
    for step in range(1,steps+1):
        vc+=.5*dt*ac;vs+=.5*dt*ast;vh+=.5*dt*ah
        uc+=dt*vc;us+=dt*vs;zh+=dt*vh
        apply_supports(uc,vc)
        ac,ast,ah,fh,nc,pen,nb,ni,mom_error,amax_s=force_step(True)
        vc+=.5*dt*ac;vs+=.5*dt*ast;vh+=.5*dt*ah
        apply_supports(uc,vc)
        if nc and first_contact is None:first_contact=step
        mc=float(np.max(np.linalg.norm(uc,axis=1)))*1000
        ms=float(np.max(np.linalg.norm(us,axis=1)))*1000
        if not all(np.isfinite(z) for z in (mc,ms,fh,amax_s)) or max(mc,ms)>limit_mm or pen>dx/2:
            error=f'safety stop step {step}: concrete={mc:.4g} mm; steel={ms:.4g} mm; penetration={pen*1000:.4g} mm'
        if step==1 or step%save_every==0 or error or (first_contact==step):
            rows.append(dict(step=step,time_ms=step*dt*1000,impact_kN=max(0,fh)/1000,
                 midspan_mm=-1000*np.mean(uc[mid,2]),max_concrete_mm=mc,max_steel_mm=ms,
                 n_contact=nc,penetration_mm=pen*1000,broken_concrete=nb,
                 broken_interface=ni,max_steel_acceleration=amax_s,
                 interface_physical_momentum_residual_N=mom_error,hammer_velocity=vh))
        if verbose and (step%50==0 or error): print(step,'force kN',round(fh/1000,3),'max concrete mm',round(mc,6),'max steel mm',round(ms,6),'interface failures',ni)
        if error:break
    summary=dict(mode=mode,steps=int(step),first_contact=first_contact,stop=error,
                 max_recorded_contact_kN=max(r['impact_kN'] for r in rows),
                 max_concrete_mm=max(r['max_concrete_mm'] for r in rows),
                 max_steel_mm=max(r['max_steel_mm'] for r in rows),elapsed_s=time.perf_counter()-t0,
                 source_status='source-density convention adapted; interface, contact and geometry provisional')
    with (STAGE17_DIR/'density_normalized_history.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    (STAGE17_DIR/'summary.json').write_text(json.dumps(summary,indent=2))
    np.savez_compressed(STAGE17_DIR/'state.npz',uc=uc,vc=vc,us=us,vs=vs,hammer_z=zh,hammer_velocity=vh)
    print('Stage17 summary:',summary)
    return rows,summary

# Conservation warning: physical interface momentum is not guaranteed by force-density parity
r=np.random.default_rng(17)
u0=r.normal(0,1.e-8,size=xyz.shape);s0=r.normal(0,1.e-8,size=steel_xyz.shape)
fc_d,fs_d,_=source_c_force(u0,s0,'literal',commit=False)
print('Interface force-density vector sum, for regression only:', np.linalg.norm(fc_d.sum(axis=0)+fs_d.sum(axis=0)))
print('Physical force residual under mixed control volumes:',np.linalg.norm(fc_d.sum(axis=0)*dx**3+np.sum(fs_d*STEEL_MASS[:,None]/p.rho_steel,axis=0)))