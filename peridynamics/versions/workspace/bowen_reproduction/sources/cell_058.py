def apply_supports(u,v):
    u[support_mask,2]=0.0
    v[support_mask,2]=0.0
    u[anchor_x,0]=0.0; v[anchor_x,0]=0.0
    u[anchor_y,1]=0.0; v[anchor_y,1]=0.0

def system_forces(uc, us, zh, dt, commit=False):
    fpd,stretch,failed=model.evaluate(uc,dt,commit=commit)
    fc,fs=coupling_forces(uc,us)
    fs+=rod_forces(us)
    cforce,hforce,epot,ncontact,maxpenetration=sphere_contact(xyz+uc,zh)
    fc+=fpd+cforce
    return fc,fs,hforce[2],epot,ncontact,maxpenetration,int(failed.sum())

def reset_concrete_history():
    model.broken[:]=False
    model.prev_stretch[:]=0.

def run_impact(num_steps=100,save_stride=10,dt=p.dt,verbose=True):
    reset_concrete_history()
    uc=np.zeros_like(xyz); vc=np.zeros_like(xyz)
    us=np.zeros_like(steel_xyz); vs=np.zeros_like(steel_xyz)
    zh=hammer_z_initial; vh=HAMMER_V0
    results=[]
    fcon,fsteel,fh,ec,ncontact,pen,nfailed=system_forces(uc,us,zh,dt,commit=True)
    acon=fcon/CONCRETE_NODE_MASS
    asteel=fsteel/STEEL_MASS[:,None]
    ah=(fh/HAMMER_MASS)-(9.81 if USE_GRAVITY else 0.)
    t0=time.perf_counter()
    midspan=(np.abs(xyz[:,0]-p.L/2)<=dx/2+1e-12) & (xyz[:,2]>p.H/2)
    for step in range(1,num_steps+1):
        vc+=0.5*dt*acon;vs+=0.5*dt*asteel;vh+=0.5*dt*ah
        uc+=dt*vc;us+=dt*vs;zh+=dt*vh
        apply_supports(uc,vc)
        # Recalculate forces with state at the NEW step.
        fcon,fsteel,fh,ec,ncontact,pen,nfailed=system_forces(uc,us,zh,dt,commit=True)
        acon=fcon/CONCRETE_NODE_MASS
        asteel=fsteel/STEEL_MASS[:,None]
        ah=(fh/HAMMER_MASS)-(9.81 if USE_GRAVITY else 0.)
        vc+=0.5*dt*acon;vs+=0.5*dt*asteel;vh+=0.5*dt*ah
        apply_supports(uc,vc)
        if step%save_stride==0 or step==1:
            kinetic=(0.5*CONCRETE_NODE_MASS*np.sum(vc**2)+0.5*np.sum(STEEL_MASS[:,None]*vs**2)+0.5*HAMMER_MASS*vh**2)
            rec=dict(step=step,time_s=step*dt,hammer_z=zh,hammer_v=vh,
                     impact_force_kN=max(fh,0.0)/1e3,contact_points=ncontact,
                     penetration_mm=pen*1e3,
                     midspan_deflection_mm=-float(np.mean(uc[midspan,2]))*1e3,
                     max_displacement_mm=float(np.max(np.linalg.norm(uc,axis=1)))*1e3,
                     broken_bonds=int(model.broken.sum()),kinetic_J=kinetic,contact_potential_J=ec)
            results.append(rec)
            if verbose and (step==1 or step%max(save_stride,50)==0):
                print(f"{step:6d}/{num_steps}: force={rec['impact_force_kN']:.3f} kN; "
                      f"deflection={rec['midspan_deflection_mm']:.5f} mm; "
                      f"contacts={ncontact}; broken={rec['broken_bonds']}")
        if not (np.isfinite(uc).all() and np.isfinite(us).all() and np.isfinite(zh)):
            raise FloatingPointError(f'Nonfinite state at step {step}')
    if verbose:print(f'Elapsed {time.perf_counter()-t0:.2f}s')
    return results,(uc,vc,us,vs,zh,vh)
