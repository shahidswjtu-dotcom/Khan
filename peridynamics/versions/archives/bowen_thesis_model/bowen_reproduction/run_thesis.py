"""Run the integrated equation-based reconstruction, without peak fitting."""
import argparse
from thesis_pd.integrator import run
p=argparse.ArgumentParser()
p.add_argument('--out',required=True);p.add_argument('--steps',type=int,default=5000)
p.add_argument('--dt',type=float,default=3e-7);p.add_argument('--resume')
a=p.parse_args();_,result=run(a.out,a.steps,'family_projection',a.dt,a.resume)
print(result)
