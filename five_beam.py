#!/usr/bin/env python3
"""Bounded exact projective deletion/regrowth beam for strict five-orchards.

All projective incidences use primitive integer cross products. All 4/4, 4/3,
3/3 intersections are proposed; 4/2, 3/2, 2/2 are sampled. Randomised beam,
coordinate bit cap, and finite rounds make this incomplete, not a proof.
"""
import argparse, itertools, json, random, time
from pathlib import Path
from beam_orchard import normalize, cross, incident, chart
from investigate_five import census

def state(P):
 P=tuple(sorted(P));h,lines=census(P)
 if max(h)>5:return None
 typ={i:[] for i in range(2,6)}
 for line,mask in lines.items():typ[mask.bit_count()].append((line,mask))
 return (P,typ,h.get(5,0))

def candidates(st,rng,samples,maxbits):
 P,typ,_=st;seen=set(P);pairs=[]
 for i,j in [(4,4),(4,3),(3,3)]:
  pairs.extend(itertools.combinations(typ[i],2) if i==j else itertools.product(typ[i],typ[j]))
 for i,j in [(4,2),(3,2),(2,2)]:
  if not typ[i] or not typ[j] or (i==j and len(typ[i])<2):continue
  for _ in range(samples):
   x,y=(rng.sample(typ[i],2) if i==j else (rng.choice(typ[i]),rng.choice(typ[j])))
   pairs.append((x,y))
 out=[]
 for (l,m),(k,n) in pairs:
  if m&n:continue
  p=cross(l,k)
  if p is None or p in seen or max(abs(v).bit_length() for v in p)>maxbits:continue
  seen.add(p)
  if any(incident(q,p) for q,_ in typ[5]):continue
  gain=sum(incident(q,p) for q,_ in typ[4]);pot=sum(incident(q,p) for q,_ in typ[3])
  out.append((p,gain,pot))
 return out

def main(a):
 rng=random.Random(a.seed);source=json.loads(Path('research/five_point_examples.json').read_text())[str(a.source)]['homogeneous_points'];source=state(tuple(map(tuple,source)));baseline=a.baseline
 start=time.time();total=0;checks=0;best=-1;winner=None;rows=[]
 for r in range(a.rounds):
  deg=[sum(bool(m>>i&1) for _,m in source[1][5]) for i in range(a.source)]
  inds=sorted(range(a.source),key=lambda i:deg[i]+rng.random()*3)[:max(a.delete*3,a.source//2)] if r%2==0 else list(range(a.source))
  removed=set(rng.sample(inds,a.delete));layer=[state([p for i,p in enumerate(source[0]) if i not in removed])]
  for step in range(a.delete+a.add):
   pool={}
   for st in layer:
    choices=candidates(st,rng,a.samples,a.maxbits);total+=len(choices)
    choices.sort(key=lambda z:z[1]+.14*z[2]+rng.random()*.2,reverse=True)
    for p,g,pot in choices[:a.branch]:
     child=state(st[0]+(p,));checks+=1
     if child is None:raise AssertionError('candidate introduced six-collinear line')
     assert child[2]==st[2]+g
     qual=child[2]+.08*len(child[1][4])+.02*len(child[1][3])+rng.random()*.1
     pool[child[0]]=(qual,child)
     if step==a.delete+a.add-1:
      if child[2]>best:best=child[2]
      if child[2]>baseline:
       baseline=child[2];winner={'points':child[0],'score':baseline,'histogram':{k:len(v) for k,v in child[1].items()},'round':r,'step':step}
       Path(a.result+'.candidate.json').write_text(json.dumps(winner,indent=2)+'\n')
       Q,cv=chart(child[0]);Path(a.result+'.candidate.txt').write_text('a = 0\n'+', '.join(f'({x}, {y})' for x,y in Q)+'\n')
       print('IMPROVEMENT',baseline,flush=True)
   layer=[v[1] for v in sorted(pool.values(),key=lambda z:z[0],reverse=True)[:a.width]]
   if not layer:break
  rows.append({'round':r+1,'best_completed':best,'proposed':total,'states_checked':checks})
  print(json.dumps(rows[-1]),flush=True)
 out={'args':vars(a),'seed_score':source[2],'best_completed':best,'final_baseline':baseline,'proposed':total,'states_checked':checks,'seconds':round(time.time()-start,2),'rounds':rows,'winner':winner}
 Path(a.result).write_text(json.dumps(out,indent=2)+'\n')
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--source',type=int,default=25);p.add_argument('--delete',type=int,default=3);p.add_argument('--add',type=int,default=0);p.add_argument('--rounds',type=int,default=12);p.add_argument('--width',type=int,default=8);p.add_argument('--branch',type=int,default=15);p.add_argument('--samples',type=int,default=100);p.add_argument('--maxbits',type=int,default=96);p.add_argument('--seed',type=int,default=17);p.add_argument('--baseline',type=int,default=18);p.add_argument('--result',required=True);main(p.parse_args())
