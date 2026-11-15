#!/usr/bin/env python3
"""MILP selection from a numerical algebraic intersection pool, with exact checks.

Selects a fixed number of additional points while preserving all seed points.
Line incidence constraints are generated numerically; consequently even solver
optimality would NOT certify the true optimum of the exact pool. Each returned
incumbent is recounted using exact number-field arithmetic before saving. The pool
and seed are explicitly bounded, not a global geometric search.
"""
import argparse,json,random,time
from pathlib import Path
import numpy as np,sympy as sp
from scipy.optimize import milp,Bounds,LinearConstraint
from scipy.sparse import coo_matrix
from algfield import parse_entry
from orchard import analyze,line_key,entry_text
from extend import to_float_pts,exact_point_from_lines
from beam_algebraic import floatscore
from beam_numeric_replay import candidates

def run(a):
 rng=random.Random(a.seed);F,P0=parse_entry(Path(a.entry).read_text());P0=list(P0);n=len(P0);r0=analyze(F,P0)
 if not r0['valid']:raise ValueError('invalid seed')
 sym=sp.Symbol('a');poly=sum(int(c)*sym**i for i,c in enumerate(F.p));root=max(complex(z).real for z in sp.nroots(poly,maxsteps=250) if abs(complex(z).imag)<1e-9);Q=to_float_pts(F,P0,root)
 pool=list(P0);approx=list(Q)
 if a.include:
  FF,known=parse_entry(Path(a.include).read_text())
  if FF.p!=F.p:raise ValueError('incompatible field')
  for q in known:
   if q not in pool:pool.append(q);approx.append(to_float_pts(F,[q],root)[0])
 sources=[a.entry]+(a.candidate_source or []);added=0;proposed_total=0
 for src in sources:
  FF,seedpts=parse_entry(Path(src).read_text())
  if FF.p!=F.p:raise ValueError('candidate source uses incompatible field')
  qpts=to_float_pts(F,seedpts,root);L=floatscore(qpts,a.tol)
  if L[5]:continue
  V=candidates(qpts,L,rng,a.pairs32,a.pairs22,a.tol);proposed_total+=len(V)
  V.sort(key=lambda z:z[0]+a.potential*z[1]+rng.random()*.15,reverse=True)
  count=0
  for gain,trips,(i,j),(k,l),v in V:
   q=exact_point_from_lines(F,line_key(F,seedpts[i],seedpts[j]),line_key(F,seedpts[k],seedpts[l]))
   if q is None or q in pool:continue
   qv=to_float_pts(F,[q],root)[0]
   if not np.all(np.isfinite(qv)) or np.linalg.norm(qv-v)>a.tol*max(1.,np.max(abs(qpts)))*10:continue
   pool.append(q);approx.append(qv);added+=1;count+=1
   if count>=a.extra:break
 points=np.array(approx);lines=floatscore(points,a.tol);rich=lines[4]+lines[5];m=len(pool);r=len(rich);rows=[];cols=[];vals=[]
 for i in range(m):rows.append(0);cols.append(i);vals.append(1.)
 for t,(_,_,idx) in enumerate(rich):
  row=1+2*t
  for i in idx:rows.extend((row,row+1));cols.extend((i,i));vals.extend((1.,-1.))
  rows.append(row+1);cols.append(m+t);vals.append(4.)
 C=coo_matrix((vals,(rows,cols)),shape=(1+2*r,m+r)).tocsc();lo=np.zeros(m+r);hi=np.ones(m+r);lo[:n]=1
 obj=np.r_[np.zeros(m),-np.ones(r)];lb=np.r_[0.,np.full(2*r,-np.inf)];ub=np.r_[0.,np.tile([4.,0.],r)]
 baseline={x['n']:x['score'] for x in json.loads(Path('records_int_now.json').read_text())['records']}
 for row in json.loads(Path('final_rows.json').read_text()):baseline[row['N']]=max(baseline.get(row['N'],0),row['score'])
 for nn in range(42,47):baseline[nn]=max(baseline.get(nn,0),3*nn-51)
 start=time.time();results=[];print(json.dumps({'source':a.entry,'pool':m,'candidate_intersections':proposed_total,'added':added,'numerical_rich_pool_lines':r,'seed':a.seed}),flush=True)
 for target in a.targets:
  lb[0]=ub[0]=float(target)
  sol=milp(obj,integrality=np.ones(m+r),bounds=Bounds(lo,hi),constraints=LinearConstraint(C,lb,ub),options={'time_limit':a.seconds,'mip_rel_gap':0.0})
  out={'target':target,'status':int(sol.status),'message':sol.message,'numerical_objective':None if sol.fun is None else -float(sol.fun),'exact_score':None,'valid':None,'baseline':baseline.get(target),'saved':None,'seconds_total':round(time.time()-start,2)}
  if sol.x is not None:
   chosen=[pool[i] for i in range(m) if sol.x[i]>.5]
   if len(chosen)!=target:raise ValueError('solver returned nonintegral selection')
   rr=analyze(F,chosen);out['valid']=rr['valid'];out['exact_score']=rr['score']
   if rr['valid'] and rr['score']>baseline.get(target,10**9):
    name=f'algebraic_milp_n{target}_s{rr["score"]}_seed{a.seed}.txt';Path('research/new_candidates').mkdir(exist_ok=True,parents=True)
    Path('research/new_candidates',name).write_text(entry_text(F,chosen,[int(v) for v in F.p])+'\n');out['saved']=name;baseline[target]=rr['score'];print('NEW',name,flush=True)
  results.append(out);print(json.dumps(out),flush=True);Path(a.result).write_text(json.dumps({'args':vars(a),'pool':m,'numerical_pool_lines':r,'results':results},indent=2)+'\n')
 print('DONE',flush=True)
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--entry',required=True);p.add_argument('--include');p.add_argument('--candidate-source',nargs='*');p.add_argument('--targets',type=int,nargs='+',required=True);p.add_argument('--extra',type=int,default=70);p.add_argument('--pairs32',type=int,default=1600);p.add_argument('--pairs22',type=int,default=1200);p.add_argument('--potential',type=float,default=.12);p.add_argument('--seconds',type=int,default=45);p.add_argument('--tol',type=float,default=1e-8);p.add_argument('--seed',type=int,default=1);p.add_argument('--result',required=True);run(p.parse_args())
