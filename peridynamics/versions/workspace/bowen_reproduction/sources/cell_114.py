RUN_STAGE17 = False  # Set True for contact-reaching diagnostic.
STAGE17_STEPS = 420
if RUN_STAGE17:
    history17,summary17=stage17_run(steps=STAGE17_STEPS,save_every=2)
else:
    print('Stage 17 is ready. Set RUN_STAGE17=True to run impact diagnostics.')