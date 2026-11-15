#!/usr/bin/env python3
"""Wide stochastic beam with exact replay of numerically promising lineages.

Each node stores floating coordinates plus the indices of two generating point-pair
lines for every addition. A fast vectorized incidence heuristic proposes points;
only apparent improvements are replayed from the exact algebraic seed and verified
with exact field arithmetic. Floating filters can miss or invent incidences. This is
an incomplete search, not a proof and not a finite-pool optimality certificate.
"""
import argparse,itertools,json,random,time
from pathlib import Path
import numpy as np
import sympy as sp
from algfield import parse_entry
from orchard import analyze,line_key,entry_text
from extend import exact_point_from_lines,to_float_pts
from beam_algebraic import floatscore

def candidates(P,L,rng,n32,n22,tol):
 two=L[2];three=L[3];four=L[4]+L[5];pairs=[]
 for A,B in itertools.combinations(three,2):
  if A[2].isdisjoint(B[2]):pairs.append((A,B))
 if two and three:
  for _ in range(n32):
   A=rng.choice(three);B=rng.choice(two)
   if A[2].isdisjoint(B[2]):pairs.append((A,B))
 if len(two)>1:
  for _ in range(n22):
   A,B=rng.sample(two,2)
   if A[2].isdisjoint(B[2]):pairs.append((A,B))
 if not pairs:return []
 a=np.array([q[0][0] for q in pairs]);b=np.array([q[1][0] for q in pairs]);det=a[:,0]*b[:,1]-a[:,1]*b[:,0]
 ok=np.abs(det)>1e-9;den=np.where(ok,det,1)
 xy=np.column_stack(((a[:,1]*b[:,2]-a[:,2]*b[:,1])/den,(a[:,2]*b[:,0]-a[:,0]*b[:,2])/den))
 scale=max(1.,float(np.max(np.abs(P))));ok &= np.max(abs(xy),axis=1)<1000*scale
 if four:
  abc=np.array([x[0] for x in four]);ok &= (np.min(np.abs(xy@abc[:,:2].T+abc[:,2]),axis=1)>=tol*scale)
 ok &= np.min(np.linalg.norm(xy[:,None,:]-P[None,:,:],axis=2),axis=1)>=tol*scale
 ids=np.flatnonzero(ok)
 if not len(ids):return []
 _,ii=np.unique(np.round(xy[ids]/(tol*scale)).astype(np.int64),axis=0,return_index=True);ids=ids[np.sort(ii)]
 if three:
  abc=np.array([x[0] for x in three]);gain=np.sum(np.abs(xy[ids]@abc[:,:2].T+abc[:,2])<tol*scale,axis=1)
 else:gain=np.zeros(len(ids),dtype=int)
 if two:
  abc=np.array([x[0] for x in two]);trips=np.sum(np.abs(xy[ids]@abc[:,:2].T+abc[:,2])<tol*scale,axis=1)
 else:trips=np.zeros(len(ids),dtype=int)
 return [(int(g),int(t),pairs[k][0][1],pairs[k][1][1],xy[k]) for k,g,t in zip(ids,gain,trips)]

def replay(F,P0,path):
 pts=list(P0)
 for (i,j),(k,l) in path:
  q=exact_point_from_lines(F,line_key(F,pts[i],pts[j]),line_key(F,pts[k],pts[l]))
  if q is None or q in pts:return None
  pts.append(q)
 res=analyze(F,pts)
 return (pts,res) if res['valid'] else None

def run(a):
 rng=random.Random(a.seed);F,P0=parse_entry(Path(a.entry).read_text());n=len(P0);res=analyze(F,P0)
 if not res['valid']:raise ValueError('invalid seed')
 sym=sp.Symbol('a');poly=sum(int(c)*sym**i for i,c in enumerate(F.p));roots=sp.nroots(poly,maxsteps=250);root=max(complex(z).real for z in roots if abs(complex(z).imag)<1e-9)
 P=to_float_pts(F,P0,root);base={x['n']:x['score'] for x in json.loads(Path('records_int_now.json').read_text())['records']}
 for x in json.loads(Path('final_rows.json').read_text()):base[x['N']]=max(base.get(x['N'],0),x['score'])
 for k in range(42,47):base[k]=max(base.get(k,0),3*k-51)
 start=time.time();counts={'states_expanded':0,'points_proposed':0,'exact_checks':0,'exact_invalid':0};reports=[];saved=[];layer=[(P,())]
 for step in range(1,a.steps+1):
  pool=[]
  for pts,path in layer:
   L=floatscore(pts,a.tol);counts['states_expanded']+=1
   if L[5]:continue
   cand=candidates(pts,L,rng,a.pairs32,a.pairs22,a.tol);counts['points_proposed']+=len(cand)
   for gain,trips,ij,kl,q in cand:
    score=len(L[4])+gain;potential=len(L[3])-gain+trips
    # The Gumbel/random jitter encourages genuinely different line topologies.
    quality=score+a.potential*potential+a.temperature*rng.gammavariate(1,1)
    pool.append((quality,score,pts,path,ij,kl,q))
  pool.sort(key=lambda x:x[0],reverse=True);layer=[];seen=set();checked=0
  for quality,score,pts,path,ij,kl,q in pool:
   # Use a coordinate-set signature to discard repeated intersections/path order.
   child=np.vstack([pts,q]);scale=max(1.,float(np.max(abs(child))))
   sig=tuple(sorted(map(tuple,np.round(child/(scale*a.tol)).astype(np.int64))))
   if sig in seen:continue
   seen.add(sig);newpath=path+((ij,kl),)
   if score>base.get(n+step,10**9) and checked<a.checks:
    checked+=1;counts['exact_checks']+=1;exact=replay(F,P0,newpath)
    if exact is None:counts['exact_invalid']+=1
    else:
     E,r=exact;actual=r['score'];threshold=base[n+step]
     print('EXACT',json.dumps({'n':n+step,'predicted':score,'score':actual,'baseline':threshold}),flush=True)
     if actual>threshold:
      name=f'numeric_replay_n{n+step}_s{actual}_seed{a.seed}_{counts["exact_checks"]}.txt'
      Path('research/new_candidates').mkdir(parents=True,exist_ok=True)
      Path('research/new_candidates',name).write_text(entry_text(F,E,[int(c) for c in F.p])+'\n');base[n+step]=actual;saved.append(name);print('NEW',name,flush=True)
   layer.append((child,newpath))
   if len(layer)>=a.width:break
  rec={'size':n+step,'beam':len(layer),'float_max':max((x[1] for x in pool),default=None),'baseline':base.get(n+step),'pool':len(pool),'counts':counts.copy(),'seconds':round(time.time()-start,2)};reports.append(rec);print(json.dumps(rec),flush=True)
  Path(a.result).write_text(json.dumps({'args':vars(a),'reports':reports,'saved':saved},indent=2)+'\n')
  if not layer:break
 print('DONE',flush=True)
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--entry',required=True);p.add_argument('--steps',type=int,default=8);p.add_argument('--width',type=int,default=40);p.add_argument('--pairs32',type=int,default=600);p.add_argument('--pairs22',type=int,default=500);p.add_argument('--potential',type=float,default=.08);p.add_argument('--temperature',type=float,default=.20);p.add_argument('--tol',type=float,default=1e-8);p.add_argument('--checks',type=int,default=4);p.add_argument('--seed',type=int,default=1);p.add_argument('--result',required=True);run(p.parse_args())
