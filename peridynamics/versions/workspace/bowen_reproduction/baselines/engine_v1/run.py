import argparse
from bowen_pd.integrator import run
p=argparse.ArgumentParser();p.add_argument("--out",required=True);p.add_argument("--steps",type=int,default=550);p.add_argument("--mode",default="legacy");p.add_argument("--dt",type=float,default=3e-7);p.add_argument("--resume");a=p.parse_args()
run(a.out,a.steps,a.mode,a.dt,a.resume)
