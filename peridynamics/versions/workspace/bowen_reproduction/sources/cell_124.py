# Configuration: full 5-mm concrete geometry, reference dt, unchanged model.
RUN_STAGE19 = True
STAGE19_STEPS = 550
STAGE19_SAVE_EVERY = 5
STAGE19_DIR = Path('bowen_stage19_outputs')
STAGE19_DIR.mkdir(exist_ok=True)
print('Time simulated (ms):', STAGE19_STEPS*p.dt*1000)
print('Constitutive/contact parameters: UNCHANGED from Stage 18')