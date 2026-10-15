if RUN_STAGE19 and stage19_data:
    import pandas as pd
    d=pd.DataFrame(stage19_data)
    import matplotlib.pyplot as plt
    fig,axs=plt.subplots(3,1,figsize=(9,10),sharex=True,constrained_layout=True)
    axs[0].plot(d.time_ms,d.impact_kN);axs[0].set_ylabel('Hammer reaction (kN)')
    axs[1].plot(d.time_ms,d.midspan_mm,label='Recorded midspan');axs[1].plot(d.time_ms,d.max_concrete_mm,label='Max particle magnitude',alpha=.7);axs[1].legend();axs[1].set_ylabel('Displacement (mm)')
    axs[2].plot(d.time_ms,d.broken_concrete,label='Concrete bonds');axs[2].plot(d.time_ms,d.broken_interface,label='PDROD samples');axs[2].legend();axs[2].set_ylabel('Failed count');axs[2].set_xlabel('Time (ms)')
    fig.savefig(STAGE19_DIR/'extended_impact_diagnostics.png',dpi=180)
    plt.show()
    print('Peak recorded force (kN):',round(d.impact_kN.max(),4))
    print('Recorded max midspan deflection magnitude (mm):',round(d.midspan_mm.abs().max(),4))
    print('End hammer velocity (m/s):',round(d.hammer_velocity_mps.iloc[-1],5))
    print('Max concrete particle displacement (mm):',round(d.max_concrete_mm.max(),4))
    print('Last recorded concrete/interface failures:',int(d.broken_concrete.iloc[-1]),int(d.broken_interface.iloc[-1]))