# Run only the short contact-reaching controlled demonstration by default.
RUN_STAGE18 = False
STAGE18_STEPS = 425
if RUN_STAGE18:
    stage18_data, stage18_summary = stage18_run(steps=STAGE18_STEPS,mode='conservative')
else:
    print('Stage 18 verification complete. Set RUN_STAGE18=True for the 425-step diagnostic.')