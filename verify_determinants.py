"""Independent incidence census by exact triple determinants (no line-key grouping).
Shares the exact field parser/arithmetic with the original verifier; it is not an
independent audit of that arithmetic. Run from the project root.
"""
from itertools import combinations
from collections import Counter
from math import comb
from algfield import parse_entry

EXPECTED={28:{2:126,4:42},36:{2:252,4:63},40:{2:240,3:40,4:70},48:{2:516,3:12,4:96}}

def census(n):
 F,P=parse_entry(open(f'entries/N{n}.txt').read())
 assert len(P)==n and len(set(P))==n
 pairs={}
 for i,j,k in combinations(range(n),3):
  x1,y1=P[i];x2,y2=P[j];x3,y3=P[k]
  z=F.sub(F.mul(F.sub(x2,x1),F.sub(y3,y1)),F.mul(F.sub(y2,y1),F.sub(x3,x1)))
  if F.is_zero(z):
   for p in ((i,j),(i,k),(j,k)):pairs.setdefault(p,set()).update((i,j,k))
 lines={frozenset(v) for v in pairs.values()}
 for line in lines:
  assert all(pairs.get(pair)==set(line) for pair in combinations(sorted(line),2))
 hist=Counter(map(len,lines));hist[2]=comb(n,2)-sum(comb(k,2)*v for k,v in hist.items())
 return dict(sorted(hist.items()))
for n,expected in EXPECTED.items():
 actual=census(n);print(n,actual,flush=True);assert actual==expected
