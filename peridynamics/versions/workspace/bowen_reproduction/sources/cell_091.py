if RUN_STAGE9:
    import matplotlib.pyplot as plt
    t=np.array([v['time_ms'] for v in records9]);
    fig,axs=plt.subplots(3,1,figsize=(9,9),sharex=True)
    axs[0].plot(t,[v['impact_force_kN'] for v in records9]);axs[0].set_ylabel('Contact force (kN)')
    axs[1].plot(t,[v['midspan_deflection_mm'] for v in records9]);axs[1].set_ylabel('Midspan displacement (mm)')
    axs[2].plot(t,[v['contact_count'] for v in records9]);axs[2].set_ylabel('Contacting particles');axs[2].set_xlabel('Time (ms)')
    fig.suptitle('Stage 9 — provisional contact model: diagnostic curves')
    fig.tight_layout();plt.show()