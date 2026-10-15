# Optimized solver: same loading and physical parameters, compiled force assembly.
# CONTACT LAW AND INTERFACE QUADRATURE REMAIN PROVISIONAL.
def fast_system_forces(uc,us,zh,dt,commit=False):
    Fc,nb=fast_concrete.evaluate(uc,dt,commit)
    Fic,Fis,ni=fast_interface.evaluate(uc,us,commit)
    Fs=rod_forces(us)+Fis
    Fc+=Fic
    C,H,Ucontact,ncontact,penetration=sphere_contact(xyz+uc,zh)
    Fc+=C
    return Fc,Fs,float(H[2]),Ucontact,ncontact,penetration,nb,ni

def run_fast_impact(steps=450,stride=5,dt=p.dt,verbose=True):
    fast_concrete.reset();fast_interface.reset()
    uc=np.zeros_like(xyz); vc=np.zeros_like(uc)
    us=np.zeros_like(steel_xyz);vs=np.zeros_like(us)
    zh=float(hammer_z_initial);vh=float(HAMMER_V0)
    FC,FS,FH,CU,NC,PEN,NB,NI=fast_system_forces(uc,us,zh,dt,False)
    ac=FC/CONCRETE_NODE_MASS;ast=FS/STEEL_MASS[:,None]
    ah=FH/HAMMER_MASS-(9.81 if USE_GRAVITY else 0.)
    mid=(np.abs(xyz[:,0]-p.L/2)<=dx/2+1.e-12)&(xyz[:,2]>=p.H/2)
    recs=[];t0=perf_counter()
    for step in range(1,steps+1):
        vc+=0.5*dt*ac;vs+=0.5*dt*ast;vh+=0.5*dt*ah
        uc+=dt*vc;us+=dt*vs;zh+=dt*vh
        apply_supports(uc,vc)
        FC,FS,FH,CU,NC,PEN,NB,NI=fast_system_forces(uc,us,zh,dt,True)
        ac=FC/CONCRETE_NODE_MASS;ast=FS/STEEL_MASS[:,None]
        ah=FH/HAMMER_MASS-(9.81 if USE_GRAVITY else 0.)
        vc+=0.5*dt*ac;vs+=0.5*dt*ast;vh+=0.5*dt*ah
        apply_supports(uc,vc)
        if not np.isfinite(uc).all() or not np.isfinite(us).all() or not np.isfinite(zh):
            raise FloatingPointError(f'Nonfinite state at step {step}')
        if np.max(np.abs(uc))>0.05 or PEN>0.5*dx:
            raise RuntimeError(f'Unphysical displacement/contact penetration at step {step}, pen={PEN:.5g} m')
        if step%stride==0 or step==1 or (NC and (not recs or recs[-1]['contact_count']==0)):
            q={'step':step,'time_ms':step*dt*1.e3,'force_kN':max(0.,FH)*1.e-3,
               'midspan_deflection_mm':-np.mean(uc[mid,2])*1.e3,'hammer_velocity_mps':vh,
               'contact_count':NC,'penetration_mm':PEN*1.e3,
               'concrete_broken':int(fast_concrete.broken.sum()),
               'interface_broken':int(fast_interface.broken.sum()),
               'kinetic_J':0.5*CONCRETE_NODE_MASS*np.sum(vc**2)+0.5*np.sum(STEEL_MASS[:,None]*vs**2)+0.5*HAMMER_MASS*vh**2,
               'contact_potential_J':CU}
            recs.append(q)
            if verbose and (step%50==0 or NC and len(recs)>1 and recs[-2]['contact_count']==0):
                print(f"step={step} t={q['time_ms']:.4f} ms contact={NC} F={q['force_kN']:.2f} kN defl={q['midspan_deflection_mm']:.4f} mm elapsed={perf_counter()-t0:.1f}s")
    return recs,(uc,vc,us,vs,zh,vh)

# Compilation warm-up, not counted in the solve runtime.
t0=perf_counter()
fast_concrete.evaluate(np.zeros_like(xyz),p.dt,False)
fast_interface.evaluate(np.zeros_like(xyz),np.zeros_like(steel_xyz),False)
print(f'Stage 10 compilation warm-up: {perf_counter()-t0:.2f}s')

# Explicit user-controlled run. Default is intentionally modest; full run uses 450.
RUN_STAGE10=False
STAGE10_STEPS=365  # 365-step first-contact interval verified; 450-step attempt failed at 371
STAGE10_STRIDE=5
if RUN_STAGE10:
    data10,state10=run_fast_impact(STAGE10_STEPS,STAGE10_STRIDE)
    import csv
    out10=Path('bowen_stage10_outputs');out10.mkdir(exist_ok=True)
    with open(out10/'impact_history.csv','w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=data10[0].keys());w.writeheader();w.writerows(data10)
    np.savez_compressed(out10/'final_state.npz', uc=state10[0],vc=state10[1],us=state10[2],vs=state10[3],
                        zh=state10[4],vh=state10[5],broken=fast_concrete.broken,interface_broken=fast_interface.broken)
    print('FIRST CONTACT:',next((x['step'] for x in data10 if x['contact_count']>0),None))
else:
    print('Set RUN_STAGE10=True to execute optimized contact-reaching simulation')