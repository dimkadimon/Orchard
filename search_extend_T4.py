"""Exact extension of the 60-point four-axis construction.
At each step, add an intersection of two current 3-point lines maximizing
its number of incident 3-point lines. Scores are certified by a full exact census.
"""
from fractions import Fraction as F
from itertools import combinations
from math import gcd, lcm
from functools import reduce
import json
P=[(F(1,180+x),F(0)) for x in range(13)] + [(F(0),F(1,60+y)) for y in range(13)] + [(F(1,240+z),F(1,240+z)) for z in range(4,21)] + [(F(1,120+w),-F(1,120+w)) for w in range(-8,9)]
def norm(v):
 d=lcm(*(x.denominator for x in v)); a=[int(x*d) for x in v]; g=reduce(gcd,a); a=[x//g for x in a]
 if next(x for x in a if x)<0:a=[-x for x in a]
 return tuple(a)
def census(P):
 lines={}
 for i,j in combinations(range(len(P)),2):
  x,y=P[i];u,v=P[j]; k=norm((y-v,u-x,x*v-u*y));lines.setdefault(k,set()).update((i,j))
 return lines
out=[]
for step in range(12):
 lines=census(P); tri=[k for k,v in lines.items() if len(v)==3]; cand={}
 for (a,b,c),(d,e,f) in combinations(tri,2):
  det=a*e-b*d
  if not det:continue
  p=(F(b*f-c*e,det),F(c*d-a*f,det))
  if p not in P: cand[p]=cand.get(p,0)+1
 if not cand:break
 p=max(cand,key=lambda p:sum(a*p[0]+b*p[1]+c==0 for a,b,c in tri))
 gain=sum(a*p[0]+b*p[1]+c==0 for a,b,c in tri)
 P.append(p); n=len(P); count=sum(len(v)>=4 for v in census(P).values()); base=(n*n+12)//28+4
 out.append({'n':n,'T4':count,'bound':base,'gain':gain,'point':list(map(str,p))})
 print(out[-1],flush=True)
with open('extend_T4_results.json','w') as f:json.dump(out,f,indent=2)
