#!/usr/bin/env python3
"""Finite-grid strict-orchard MILP, independent of record-seeded beam search.

Select n lattice points from a rectangle and maximize the number of lines containing
exactly four of them, subject to no line containing five. Full pool line constraints
are enumerated using exact integer projective arithmetic. Time-limited results and
solver optimality statements apply ONLY to the chosen finite grid, not the real plane.
"""
import argparse, itertools, json, math, time
from pathlib import Path
import numpy as np
from scipy.optimize import milp,Bounds,LinearConstraint
from scipy.sparse import coo_matrix
from beam_orchard import cross,normalize,make_state,save_candidate

def solve(n,w,h,time_limit,result):
 t=time.time(); P=tuple(normalize((x,y,1)) for x in range(w) for y in range(h));m=len(P)
 line_masks={}
 for i,j in itertools.combinations(range(m),2):
  l=cross(P[i],P[j]);line_masks[l]=line_masks.get(l,0)|(1<<i)|(1<<j)
 lines=[v for v in line_masks.values() if v.bit_count()>=4];L=len(lines);nv=m+L
 rows=list(range(0)) ;cols=[];vals=[];lb=[float(n)];ub=[float(n)]
 for i in range(m):rows.append(0);cols.append(i);vals.append(1.)
 for j,mask in enumerate(lines):
  r=1+2*j
  while mask:
   v=mask&-mask;i=v.bit_length()-1;mask-=v
   rows.extend((r,r+1));cols.extend((i,i));vals.extend((1.,-1.))
  rows.append(r+1);cols.append(m+j);vals.append(4.)
  lb.extend((-np.inf,-np.inf));ub.extend((4.,0.))
 C=coo_matrix((vals,(rows,cols)),shape=(1+2*L,nv)).tocsc()
 objective=np.r_[np.zeros(m),-np.ones(L)]
 print(json.dumps({'n':n,'grid':[w,h],'points':m,'rich_pool_lines':L,'time_limit_s':time_limit}),flush=True)
 sol=milp(objective,integrality=np.ones(nv),bounds=Bounds(np.zeros(nv),np.ones(nv)),constraints=LinearConstraint(C,np.array(lb),np.array(ub)),options={'time_limit':time_limit,'mip_rel_gap':0.0})
 out={'n':n,'grid':[w,h],'status':int(sol.status),'message':sol.message,'objective':None if sol.fun is None else -float(sol.fun),'mip_dual_bound':getattr(sol,'mip_dual_bound',None),'seconds':round(time.time()-t,2),'score':None,'selected_points':None,'candidate':None}
 if sol.x is not None:
  pts=[P[i] for i in range(m) if sol.x[i]>.5]
  if len(pts)!=n:raise ValueError('nonintegral point selection')
  st=make_state(pts)
  if st is None:raise ValueError('five collinear')
  out['score']=st[2];out['selected_points']=[list(p) for p in st[0]]
  tab=json.loads(Path('records_int_now.json').read_text()); b=next(x['score'] for x in tab['records'] if x['n']==n)
  for row in json.loads(Path('final_rows.json').read_text()):
   if row['N']==n:b=max(b,row['score'])
  out['baseline']=b
  if st[2]>b:out['candidate']=save_candidate(st[0],st[2],n,{'method':'finite-grid MILP','grid':[w,h]},'research/candidates')
 Path(result).write_text(json.dumps(out,indent=2)+'\n');print('RESULT',json.dumps(out),flush=True)
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--n',type=int,required=True);ap.add_argument('--width',type=int,required=True);ap.add_argument('--height',type=int,required=True);ap.add_argument('--seconds',type=int,default=60);ap.add_argument('--result',required=True);a=ap.parse_args();solve(a.n,a.width,a.height,a.seconds,a.result)
