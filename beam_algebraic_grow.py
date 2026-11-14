#!/usr/bin/env python3
"""Grow exact algebraic seeds with a score/potential beam (bounded heuristic).

Generates 3/3, 3/2, 2/2 line intersections numerically, creates selected points
exactly over the input field, and counts every retained state exactly. Tests sizes
n+1..n+steps against pinned strict baselines, saving improvements. Only affine
intersections are proposed. Does not establish optimality or literature novelty.
"""
import argparse,json,random,time
from pathlib import Path
import sympy as sp
import numpy as np
from algfield import parse_entry
from orchard import analyze, line_key, entry_text
from extend import exact_point_from_lines,to_float_pts
from beam_algebraic import floatscore,proposals

def grow(n,steps,width,branch,pairs32,pairs22,seed,tol,result,entry=None):
 rng=random.Random(seed);F,P0=parse_entry(Path(entry or f'entries/N{n}.txt').read_text());pts=tuple(P0)
 if len(pts)!=n:raise ValueError('entry size does not match --n')
 res=analyze(F,pts)
 if not res['valid']:raise ValueError('invalid seed')
 a=sp.Symbol('a');poly=sum(int(v)*a**i for i,v in enumerate(F.p));roots=sp.nroots(poly,maxsteps=250);real=[complex(x).real for x in roots if abs(complex(x).imag)<1e-9]
 if not real:raise ValueError('no real root')
 root=max(real);ints=json.loads(Path('records_int_now.json').read_text());base={x['n']:x['score'] for x in ints['records']}
 for row in json.loads(Path('final_rows.json').read_text()):base[row['N']]=max(base.get(row['N'],0),row['score'])
 layer=[(pts,res)];start=time.time();saved=[];reports=[];seen=set();counts={'proposed':0,'exact_checks':0,'valid_children':0}
 for level in range(steps):
  nxt=[]
  for points,r in layer:
   P=to_float_pts(F,points,root);ls=floatscore(P,tol)
   if ls[5]:continue
   opts=proposals(P,ls,rng,pairs32,pairs22,tol);counts['proposed']+=len(opts)
   opts.sort(key=lambda v:v[0]+0.12*min(v[1],12)+rng.random()*.12,reverse=True)
   for gain,trips,(i,j),(k,l),approx in opts[:branch]:
    q=exact_point_from_lines(F,line_key(F,points[i],points[j]),line_key(F,points[k],points[l]))
    if q is None or q in points:continue
    qp=to_float_pts(F,[q],root)[0]
    if not np.all(np.isfinite(qp)) or np.linalg.norm(qp-approx)>tol*max(1.,np.max(abs(P)))*10:continue
    child=points+(q,);sig=tuple(sorted(child))
    if sig in seen:continue
    seen.add(sig);r2=analyze(F,child);counts['exact_checks']+=1
    if not r2['valid']:continue
    counts['valid_children']+=1;score=r2['score'];nn=len(child);old=base.get(nn)
    if old is not None and score>old:
     name=f'algebraic_grow_n{nn}_s{score}_seed{seed}_{counts["exact_checks"]}'
     Path('research/candidates').mkdir(parents=True,exist_ok=True)
     Path('research/candidates',name+'.txt').write_text(entry_text(F,child,[int(v) for v in F.p])+'\n')
     saved.append(name);base[nn]=score;print('NEW',name,flush=True)
    quality=score+0.08*len(r2['sizes'].get(3,[]))+rng.random()*.03
    nxt.append((quality,child,r2))
  nxt.sort(key=lambda x:x[0],reverse=True);layer=[(p,r) for _,p,r in nxt[:width]]
  report={'source':n,'size':n+level+1,'beam_retained':len(layer),'max_score':max((r['score'] for p,r in layer),default=None),'baseline':base.get(n+level+1),'counts':counts.copy(),'seconds':round(time.time()-start,2)}
  reports.append(report);print(json.dumps(report),flush=True)
  if not layer:break
 out={'source':n,'entry':entry,'steps':steps,'width':width,'branch':branch,'pairs32':pairs32,'pairs22':pairs22,'seed':seed,'tol':tol,'reports':reports,'saved':saved,'seconds':round(time.time()-start,2)}
 Path(result).write_text(json.dumps(out,indent=2)+'\n');print('RESULT',json.dumps(out),flush=True)
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--n',type=int,required=True);p.add_argument('--entry',default=None);p.add_argument('--steps',type=int,default=3);p.add_argument('--width',type=int,default=3);p.add_argument('--branch',type=int,default=8);p.add_argument('--pairs32',type=int,default=500);p.add_argument('--pairs22',type=int,default=700);p.add_argument('--seed',type=int,default=1);p.add_argument('--tol',type=float,default=1e-8);p.add_argument('--result',required=True);a=p.parse_args();grow(a.n,a.steps,a.width,a.branch,a.pairs32,a.pairs22,a.seed,a.tol,a.result,a.entry)
