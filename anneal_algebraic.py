#!/usr/bin/env python3
"""Stochastic fixed-size delete-and-replace walk in an algebraic incidence space.

Numerical geometry proposes and screens moves; every reported improvement is
replayed exactly from a trusted seed, then fully recounted in the number field.
Metropolis acceptance permits downhill moves, unlike greedy regrowth. Resets cap
lineage length. This is randomized and incomplete, not an optimality proof.
"""
import argparse,json,math,random,time
from pathlib import Path
import numpy as np,sympy as sp
from algfield import parse_entry
from orchard import analyze,line_key,entry_text
from extend import to_float_pts,exact_point_from_lines
from beam_algebraic import floatscore
from beam_numeric_replay import candidates

def replay(F,base,events):
 P=list(base)
 for rem,ij,kl in events:
  del P[rem]
  q=exact_point_from_lines(F,line_key(F,P[ij[0]],P[ij[1]]),line_key(F,P[kl[0]],P[kl[1]]))
  if q is None or q in P:return None
  P.append(q)
 r=analyze(F,P)
 return (P,r) if r['valid'] else None

def main(a):
 rng=random.Random(a.seed);F,seed=parse_entry(Path(a.entry).read_text());n=len(seed);r=analyze(F,seed)
 if not r['valid']:raise ValueError('invalid seed')
 x=sp.Symbol('a');p=sum(int(c)*x**i for i,c in enumerate(F.p));real=[complex(z).real for z in sp.nroots(p,maxsteps=250) if abs(complex(z).imag)<1e-9];root=max(real)
 base=max(r['score'],next((v['score'] for v in json.loads(Path('records_int_now.json').read_text())['records'] if v['n']==n),0))
 for v in json.loads(Path('final_rows.json').read_text()):
  if v['N']==n:base=max(base,v['score'])
 if 42<=n<=46:base=max(base,3*n-51)
 numerical_seed=to_float_pts(F,seed,root);P=numerical_seed.copy();events=[];best_float=r['score'];stats={'attempts':0,'accepted':0,'downhill_accepted':0,'proposed':0,'exact_checks':0,'exact_rejected':0,'resets':0};saved=[];start=time.time();prev=0
 for move in range(a.moves):
  if len(events)>=a.reset:
   P=numerical_seed.copy();events=[];stats['resets']+=1
  L=floatscore(P,a.tol)
  if L[5]:P=numerical_seed.copy();events=[];stats['resets']+=1;continue
  old=len(L[4]);degrees=[0]*n
  for _,_,ids in L[4]:
   for i in ids:degrees[i]+=1
  if rng.random()<.7:
   choices=sorted(range(n),key=lambda i:degrees[i]+rng.random()*3)[:max(5,n//3)];rem=rng.choice(choices)
  else:rem=rng.randrange(n)
  Q=np.delete(P,rem,axis=0);M=floatscore(Q,a.tol)
  V=candidates(Q,M,rng,a.pairs32,a.pairs22,a.tol);stats['attempts']+=1;stats['proposed']+=len(V)
  if not V:continue
  V.sort(key=lambda v:v[0]+a.potential*v[1]+rng.random()*.25,reverse=True)
  # Select across the best few options; geometrically distributed ranks encourage variety.
  rank=min(len(V)-1,int(rng.expovariate(1/max(1,a.top/4))));gain,trips,ij,kl,q=V[rank]
  score=len(M[4])+gain;delta=score-old
  temperature=a.temperature*(.5+.5*((move%a.cycle)/max(1,a.cycle-1)))
  if delta<0 and rng.random()>=math.exp(delta/temperature):continue
  P=np.vstack((Q,q));events.append((rem,ij,kl));stats['accepted']+=1;stats['downhill_accepted']+=delta<0;best_float=max(best_float,score)
  if score>base:
   stats['exact_checks']+=1;ans=replay(F,seed,events)
   if ans is None:stats['exact_rejected']+=1;print('EXACT_REJECT invalid',flush=True)
   else:
    E,rr=ans;actual=rr['score'];print('EXACT',json.dumps({'predicted':score,'actual':actual,'baseline':base,'events':len(events)}),flush=True)
    if actual>base:
     name=f'anneal_n{n}_s{actual}_seed{a.seed}_move{move+1}.txt';Path('research/new_candidates').mkdir(parents=True,exist_ok=True)
     Path('research/new_candidates',name).write_text(entry_text(F,E,[int(c) for c in F.p])+'\n');saved.append(name);base=actual;print('NEW',name,flush=True)
     seed=E;numerical_seed=to_float_pts(F,seed,root);P=numerical_seed.copy();events=[]
  if (move+1)%a.report==0:
   print(json.dumps({'move':move+1,'score_float':score,'best_float':best_float,'baseline':base,'stats':stats.copy(),'seconds':round(time.time()-start,2)}),flush=True)
 out={'args':vars(a),'n':n,'baseline_at_end':base,'best_float':best_float,'stats':stats,'saved':saved,'seconds':round(time.time()-start,2)}
 Path(a.result).write_text(json.dumps(out,indent=2)+'\n');print('RESULT',json.dumps(out),flush=True)
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--entry',required=True);p.add_argument('--moves',type=int,default=1000);p.add_argument('--pairs32',type=int,default=500);p.add_argument('--pairs22',type=int,default=700);p.add_argument('--temperature',type=float,default=1.0);p.add_argument('--cycle',type=int,default=100);p.add_argument('--potential',type=float,default=.08);p.add_argument('--top',type=int,default=25);p.add_argument('--reset',type=int,default=25);p.add_argument('--report',type=int,default=100);p.add_argument('--tol',type=float,default=1e-8);p.add_argument('--seed',type=int,default=1);p.add_argument('--result',required=True);main(p.parse_args())
