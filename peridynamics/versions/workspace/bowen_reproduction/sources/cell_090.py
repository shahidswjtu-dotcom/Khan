RUN_STAGE9 = False
STAGE9_STEPS = 450       # geometric contact ~341; just 109 steps beyond contact
STAGE9_STRIDE = 5

if RUN_STAGE9:
    started=perf_counter()
    records9, final_state9=run_stage8(STAGE9_STEPS, STAGE9_STRIDE, verbose=False)
    elapsed9=perf_counter()-started
    outdir=Path('bowen_stage9_outputs')
    outdir.mkdir(exist_ok=True)
    with (outdir/'impact_history.csv').open('w',newline='') as fp:
        writer=csv.DictWriter(fp,fieldnames=list(records9[0]))
        writer.writeheader();writer.writerows(records9)
    np.savez_compressed(outdir/'final_state.npz',
        uc=final_state9[0],vc=final_state9[1],us=final_state9[2],
        vs=final_state9[3],hammer_z=final_state9[4],hammer_v=final_state9[5],
        concrete_broken=concrete8.broken, interface_broken=interface8.broken,
        concrete_prev_stretch=concrete8.prev_stretch)
    r=np.array([v['impact_force_kN'] for v in records9]);
    d=np.array([v['midspan_deflection_mm'] for v in records9]);
    nc=np.array([v['contact_count'] for v in records9]);
    print('Wall time:',elapsed9,'s')
    print('First sampled contact step:',next((v['step'] for v in records9 if v['contact_count']>0),None))
    print('Maximum *sampled* impact force:',r.max(),'kN')
    print('Maximum *sampled* deflection magnitude:',np.abs(d).max(),'mm')
    print('Outputs:',outdir.resolve())
else:
    print('Stage 9 full-size contact run not executed. Set RUN_STAGE9=True to run.')