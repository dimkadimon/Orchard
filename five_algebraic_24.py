#!/usr/bin/env python3
"""Reconstruct Du's posted 24-point/18-line five-orchard over Q(a).

Source: https://bbs.emath.ac.cn/thread-50370-3-1.html, post 23 (Oct. 2026).
We transform its projective coordinates to one finite affine chart over Q(a).
"""
from pathlib import Path
from sympy import symbols, invert, rem, expand
from verify_independent_algebraic import census

a=symbols('a'); f=a**3+a**2-1
P=[(-a*a+a+1,a*a,1),(0,a*a,1),(a*a+a,-a*a-a+2,1),(0,-a*a+1,1),
   (a*a+a,0,1),(-a*a+a+1,2*a*a-1,1),(a*a+a,-a*a+1,1),(1,a-1,0),
   (a+1,a*a,1),(1,0,0),(a*a+a,1-a,1),(a+1,0,1),(1,0,1),(1,a,1),
   (a*a,-a*a+1,1),(0,0,1),(1,a,0),(1,-a*a+1,1),(0,1,1),
   (1,a*a+a-1,1),(1,-1,0),(a,1-a,1),(0,1,0),(a,a*a,1)]

def mod(x):return expand(rem(x,f,a))
for b in range(-5,6):
 for c in range(-5,6):
  den=[mod(z+b*x+c*y) for x,y,z in P]
  if all(q!=0 for q in den):break
 else:continue
 break
Q=[(mod(x*invert(d,f,a)),mod(y*invert(d,f,a))) for (x,y,z),d in zip(P,den)]
path=Path('entries/five_Du24.txt');path.write_text('a^3+a^2-1 = 0\n'+', '.join(f'({x}, {y})' for x,y in Q).replace('**','^')+'\n')
h=census(str(path));assert h.get(5)==18 and max(h)==5
print('Du 24, no six:',h,'affine chart z+',b,'x+',c,'y')

from algfield import parse_entry
from orchard import analyze, entry_text
from itertools import combinations
import json
F,pts=parse_entry(path.read_text());res=analyze(F,pts)
h2={k:len(v) for k,v in res['sizes'].items()};assert h2==h
rich=[sum(1<<i for i in ids) for ids in res['sizes'][5]]
result={}
for n in range(18,24):
 best=-1;idsbest=None
 for ids in combinations(range(24),n):
  m=sum(1<<i for i in ids);s=sum((v&m)==v for v in rich)
  if s>best:best,idsbest=s,ids
 subset=[pts[i] for i in idsbest]
 out=Path(f'entries/five_Du24_subset_N{n}.txt');out.write_text(entry_text(F,subset,[-1,0,1,1])+'\n')
 hh=census(str(out));assert hh.get(5,0)==best and max(hh)<=5
 result[n]={'score':best,'indices_zero_based':idsbest,'histogram':hh,'entry':str(out)}
 print('exhaustive deletions',n,best)
Path('research/five_Du24_deletions.json').write_text(json.dumps({'source':'https://bbs.emath.ac.cn/thread-50370-3-1.html post 23','seed_histogram':h,'method':'exhaustive subsets of this one seed, not global','results':result},indent=2)+'\n')
