"""Separate numerical verification from reproduction at printed precision."""
import argparse,json
from pathlib import Path
p=argparse.ArgumentParser()
p.add_argument('comparison')
p.add_argument('--out',default='outputs/reproduction_assessment.json')
a=p.parse_args()
data=json.loads(Path(a.comparison).read_text())
force=float(data['peak_force_kN']); displacement=float(data['max_observed_bottom_midspan_mm'])
force_matches=round(force,1)==359.8
displacement_matches=round(displacement,2)==22.70
assessment={
    'observed_force_kN':force,'target_force_kN':359.8,
    'force_matches_printed_precision':force_matches,
    'observed_bottom_midspan_mm':displacement,'target_max_deflection_mm':22.70,
    'observed_displacement_matches_printed_precision':displacement_matches,
    'original_interface_inputs_recovered':False,
    'original_displacement_measurement_matched':data.get('deflection_sensor_matched_to_thesis',False),
    'full_event_peak_established':data.get('full_event_peak_established',False),
    'published_crack_pattern_reproduced':False,
    'original_author_run_reproduction_demonstrated':False,
    'status':'NOT DEMONSTRATED',
    'reason':'The computed window and declared hypotheses do not establish the original Chapter 3 executable definition, original measurement, complete-event peaks or crack morphology. Passing component tests is numerical verification, not reproduction.'}
out=Path(a.out); out.parent.mkdir(parents=True,exist_ok=True)
out.write_text(json.dumps(assessment,indent=2)+'\n')
print(json.dumps(assessment,indent=2))
