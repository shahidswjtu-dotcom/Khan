# This is intentionally OFF. Stage 6 tests above run without a costly full dynamic simulation.
RUN_STAGE6_DIAGNOSTIC=False
STAGE6_DIAGNOSTIC_STEPS=first_step+100
STAGE6_SAVE_EVERY=10
if RUN_STAGE6_DIAGNOSTIC:
    dynamic_rows, dynamic_state=run_impact(
        num_steps=STAGE6_DIAGNOSTIC_STEPS,save_stride=STAGE6_SAVE_EVERY,verbose=True)
    if dynamic_rows:
        with (STAGE6_DIR/'diagnostic_history.csv').open('w',newline='') as f:
            writer=csv.DictWriter(f,fieldnames=list(dynamic_rows[0]))
            writer.writeheader();writer.writerows(dynamic_rows)
        max_force=max(r['impact_force_kN'] for r in dynamic_rows)
        max_defl=max(r['midspan_deflection_mm'] for r in dynamic_rows)
        print(f'DIAGNOSTIC ONLY: max stored contact force {max_force:.3f} kN')
        print(f'DIAGNOSTIC ONLY: max stored midspan deflection {max_defl:.3f} mm')
        if max(r['contact_points'] for r in dynamic_rows)==0:
            print('WARNING: no contact recorded; examine exact geometry and save frequency')
else:
    print('Diagnostic time integration not run. Set RUN_STAGE6_DIAGNOSTIC=True to run it.')
