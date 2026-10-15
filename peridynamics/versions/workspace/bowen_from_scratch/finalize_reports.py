"""Assemble complete fixed-window histories and source-labelled reports."""
from pathlib import Path
import csv,json,os,shutil,subprocess,sys,time
root=Path(__file__).parent
required=[root/'outputs'/name/'summary.json' for name in ('hooked_dt3e7','hooked_dt1p5e7','hooked_printed_history_complete')]
while not all(path.exists() for path in required): time.sleep(2)
combined=root/'outputs/hooked_printed_combined'; combined.mkdir(exist_ok=True)
segments=[root/'outputs/hooked_printed_history',root/'outputs/hooked_printed_history_complete']
rows=[]
for folder in segments:
    with (folder/'history.csv').open() as stream: rows.extend(csv.DictReader(stream))
with (combined/'history.csv').open('w',newline='') as stream:
    writer=csv.DictWriter(stream,fieldnames=rows[0].keys()); writer.writeheader(); writer.writerows(rows)
for name in ('manifest.json','checkpoint.npz','final_state.npz'):
    shutil.copyfile(segments[-1]/name,combined/name)
for folder in segments:
    for path in folder.glob('snapshot_*.npz'): shutil.copyfile(path,combined/path.name)
shutil.copyfile(segments[0]/'checkpoint.npz',combined/'snapshot_01500.npz')
shutil.copyfile(segments[-1]/'checkpoint.npz',combined/'snapshot_05000.npz')
summary=dict(final=rows[-1],combined_peak_force_kN=max(float(row['force_N']) for row in rows)/1000,
             combined_max_bottom_midspan_mm=max(float(row['bottom_midspan_down_m']) for row in rows)*1000,
             original_author_model_reproduced=False,history='instantaneous printed-force interpretation')
(combined/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
env=dict(os.environ,MPLCONFIGDIR='/tmp/bowen-mpl',XDG_CACHE_HOME='/tmp/bowen-cache')
commands=[['analyze.py'],['analyze_connectivity.py'],
          ['assess_reproduction.py','outputs/hooked_report/comparison.json','--out','outputs/hooked_report/reproduction_assessment.json'],
          ['analyze.py','--coarse','outputs/hooked_printed_combined','--fine','outputs/hooked_dt1p5e7','--report-dir','outputs/printed_report'],
          ['analyze_connectivity.py','--checkpoint','outputs/hooked_printed_combined/checkpoint.npz','--out','outputs/printed_report/connectivity.json'],
          ['assess_reproduction.py','outputs/printed_report/comparison.json','--out','outputs/printed_report/reproduction_assessment.json']]
for command in commands:
    subprocess.run([sys.executable,*command],cwd=root,env=env,check=True)
print('Reports finalized',flush=True)
