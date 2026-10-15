RUN_DYNAMIC=False      # Set True only when ready to run full model
DYNAMIC_STEPS=400      # Initial trial reaches contact if solver remains stable
SAVE_EVERY=10
if RUN_DYNAMIC:
    impact_records,impact_state=run_impact(num_steps=DYNAMIC_STEPS,save_stride=SAVE_EVERY)
    tt=np.array([r['time_s'] for r in impact_records])*1e3
    fig,ax=plt.subplots(2,1,figsize=(10,7),sharex=True,constrained_layout=True)
    ax[0].plot(tt,[r['impact_force_kN'] for r in impact_records]);ax[0].set_ylabel('Contact force (kN)')
    ax[1].plot(tt,[r['midspan_deflection_mm'] for r in impact_records]);ax[1].set_ylabel('Top-midspan deflection (mm)')
    ax[1].set_xlabel('Time (ms)')
    plt.show()
    print('PROTOTYPE peak force (kN)',max(r['impact_force_kN'] for r in impact_records))
    print('PROTOTYPE peak deflection (mm)',max(r['midspan_deflection_mm'] for r in impact_records))
else:
    print('Full 3-D dynamics NOT RUN. Set RUN_DYNAMIC=True after reviewing assumptions.')
    print('Expected first-contact step:',int(np.ceil(INITIAL_GAP/(abs(HAMMER_V0)*p.dt))))
