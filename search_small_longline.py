#!/usr/bin/env python3
"""Bounded exact grid search for the long-line case of T4(14).

Enumerates 9-point grid subsets, tests whether they have ten triple lines and no
four, then considers projective lines containing pairwise triple-line intersections
as potential supporting lines for five more points. This is NOT a complete search
of all real configurations. All operations use exact integer projective arithmetic.
"""
import argparse,itertools,json,time
from pathlib import Path
from beam_orchard import normalize,cross,incident,chart

def run(w,h,result):
 P=[normalize((x,y,1)) for x in range(w) for y in range(h)];D={}
 for i,j in itertools.combinations(range(len(P)),2):
  k=cross(P[i],P[j]);D[k]=D.get(k,0)|(1<<i)|(1<<j)
 rich=[(k,m) for k,m in D.items() if m.bit_count()>=3];start=time.time();seen=0;seeds=0;axes=0;best=0;out=None
 for ids in itertools.combinations(range(len(P)),9):
  mask=sum(1<<i for i in ids);triples=[];bad=False
  for k,m in rich:
   r=(m&mask).bit_count()
   if r>=4:bad=True;break
   if r==3:triples.append(k)
  seen+=1
  if bad or len(triples)!=10:continue
  seeds+=1;outside={P[i] for i in ids};hits={}
  for i,j in itertools.combinations(range(10),2):
   p=cross(triples[i],triples[j])
   if p not in outside:hits[p]=hits.get(p,0)|(1<<i)|(1<<j)
  points=list(hits);candidate_lines=set()
  for i,j in itertools.combinations(range(len(points)),2):
   ell=cross(points[i],points[j])
   if ell in candidate_lines or any(incident(ell,q) for q in outside):continue
   candidate_lines.add(ell);axes+=1;covered=0;on=[]
   for q in points:
    if incident(ell,q):covered|=hits[q];on.append(q)
   c=covered.bit_count();best=max(best,c)
   if c==10:
    chosen=None
    for size in range(1,min(5,len(on))+1):
     for subset in itertools.combinations(on,size):
      union=0
      for q in subset:union|=hits[q]
      if union.bit_count()==10:chosen=list(subset);break
     if chosen is not None:break
    if chosen is None:continue
    on=chosen
    # Extend line with exact projective points if necessary.
    for t in range(-100,101):
     if len(on)>=5:break
     a,b,c0=ell
     q=normalize((b*t-c0,-a*t, a)) if a else normalize((t,-c0,b))
     if q is not None and incident(ell,q) and q not in outside and q not in on:on.append(q)
    if len(on)!=5:continue
    allpts=list(outside)+on;aff,chartvec=chart(allpts)
    out={'grid':[w,h],'seed_indices':ids,'seed_projective':list(outside),'added_projective':on,'chart':chartvec,'integer_points':aff,'covered_lines':10};break
  if out:break
 data={'grid':[w,h],'tested_subsets':seen,'ten_triple_seeds':seeds,'axes_tested':axes,'best_covered':best,'found':out,'seconds':round(time.time()-start,2)}
 Path(result).write_text(json.dumps(data,indent=2)+'\n');print(json.dumps(data),flush=True)
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--width',type=int,required=True);p.add_argument('--height',type=int,required=True);p.add_argument('--result',required=True);a=p.parse_args();run(a.width,a.height,a.result)
