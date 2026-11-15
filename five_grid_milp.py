#!/usr/bin/env python3
"""Solve the at-least-five problem on a bounded 5 by 8 integer grid.

Binary selected-point and rich-line variables; line y requires >=5 selected
points. A zero MIP gap certifies optimality ONLY inside this finite grid.
"""
import json,math
from itertools import combinations
from pathlib import Path
import numpy as np
from scipy.optimize import milp, Bounds, LinearConstraint
from scipy.sparse import lil_matrix
from investigate_five import census
from verify_five_independent import census as independent_census
P=[(x,y,1) for x in range(8) for y in range(5)];_,lines=census(P)
rich=[m for m in lines.values() if m.bit_count()>=5];n=len(P);m=len(rich)
A=lil_matrix((m+1,n+m));A[0,:n]=1
for j,mask in enumerate(rich):
 for i in range(n):
  if mask>>i&1:A[1+j,i]=-1
 A[1+j,n+j]=5
c=np.r_[np.zeros(n),-np.ones(m)];lb=np.r_[25,np.full(m,-np.inf)];ub=np.r_[25,np.zeros(m)]
r=milp(c,integrality=np.ones(n+m),bounds=Bounds(np.zeros(n+m),np.ones(n+m)),constraints=LinearConstraint(A.tocsr(),lb,ub),options={'time_limit':120,'mip_rel_gap':0})
if r.x is None:raise RuntimeError(r.message)
Q=[P[i] for i in range(n) if r.x[i]>.5];hist,_=census(Q)
assert len(Q)==25
assert independent_census(Q)==hist
out={'grid':'{0,...,7} x {0,...,4}','n':25,'grid_rich_lines':m,'solver_status':int(r.status),'message':r.message,'reported_mip_gap':r.mip_gap,'points':Q,'histogram':hist,'exactly_five':hist.get(5,0),'at_least_five':sum(v for k,v in hist.items() if k>=5),'no_six':max(hist)<=5}
assert round(-r.fun)==out['at_least_five']
Path('research/five_grid_milp25.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
