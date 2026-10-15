if RUN_STAGE19:
    stage19_data, stage19_summary = stage18_run(
        steps=STAGE19_STEPS, dt=p.dt, save_every=STAGE19_SAVE_EVERY,
        mode='conservative', limit_mm=50, verbose=True)
    import shutil
    shutil.copy2(STAGE18_DIR/'conservative_history.csv',STAGE19_DIR/'extended_impact_history.csv')
    shutil.copy2(STAGE18_DIR/'conservative_summary.json',STAGE19_DIR/'extended_impact_summary.json')
    print('Stage 19 completed:', stage19_summary)
else:
    print('Extended run disabled')