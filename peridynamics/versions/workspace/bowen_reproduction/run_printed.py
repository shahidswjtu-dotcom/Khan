"""Executable printed-force-history alternative; not certified Chapter3 identity."""
import argparse
from printed_pd.integrator import run
p=argparse.ArgumentParser();p.add_argument('--out',required=True);p.add_argument('--steps',type=int,default=550);p.add_argument('--dt',type=float,default=3e-7);p.add_argument('--resume');a=p.parse_args()
_,s=run(a.out,a.steps,'family_projection',a.dt,a.resume);print(s)
