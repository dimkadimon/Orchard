#!/usr/bin/env python3
"""Independent exact determinant census for polynomial-coordinate certificates.

Uses SymPy only to parse polynomials and certify irreducibility and a real root;
collinearity is computed by a separate integer-polynomial quotient implementation,
without algfield.py, orchard.py or verify.py. Input must be a trusted project entry
file: SymPy's parser is not a sandbox. Does not verify cross-ratio honesty.
"""
import itertools,math,re,sys
from collections import Counter
from sympy import Symbol, Poly, QQ, ilcm
from sympy.parsing.sympy_parser import parse_expr, standard_transformations, implicit_multiplication_application

A=Symbol('a'); trans=standard_transformations+(implicit_multiplication_application,)
def expr(s):return parse_expr(s.replace('^','**'),local_dict={'a':A},transformations=trans)
def read_entry(path):
 text=open(path).read();head,body=text.split('\n',1);poly=Poly(expr(head.split('=')[0]),A,domain=QQ)
 if not poly.is_irreducible or not poly.count_roots(-float('inf'),float('inf')):raise ValueError('polynomial must be irreducible with a real root')
 poly=poly.monic();degree=poly.degree();coeff=[poly.nth(i) for i in range(degree+1)]
 if any(c.q!=1 for c in coeff):raise ValueError('expected monic integer polynomial')
 p=[int(c) for c in coeff]
 raw=re.findall(r'\(([^(),]+),\s*([^(),]+)\)',body)
 if not raw:raise ValueError('no points parsed')
 points=[];den=1
 for u,v in raw:
  coords=[]
  for s in (u,v):
   f=Poly(expr(s),A,domain=QQ).rem(poly); vals=[f.nth(i) for i in range(degree)]
   for x in vals:den=int(ilcm(den,x.q))
   coords.append(vals)
  points.append(coords)
 P=[tuple(tuple(int(z*den) for z in v) for v in pair) for pair in points]
 return p,P

def mul(u,v,p):
 d=len(u);r=[0]*(2*d-1)
 for i,x in enumerate(u):
  if x:
   for j,y in enumerate(v):r[i+j]+=x*y
 for k in range(2*d-2,d-1,-1):
  t=r[k]
  if t:
   for j in range(d):r[k-d+j]-=t*p[j]
 return r[:d]

def census(path):
 p,P=read_entry(path);n=len(P);d=len(p)-1
 if len(set(P))!=n:raise ValueError('duplicate points')
 pairs={};z=(0,)*d
 # Cache coordinate differences; arithmetic is integer-polynomial modulo p.
 dif={(i,j):(tuple(b-a for a,b in zip(P[i][0],P[j][0])),tuple(b-a for a,b in zip(P[i][1],P[j][1]))) for i,j in itertools.combinations(range(n),2)}
 for i,j,k in itertools.combinations(range(n),3):
  x,y=dif[i,j];u,v=dif[i,k]
  if mul(x,v,p)==mul(y,u,p):
   for pair in ((i,j),(i,k),(j,k)):pairs.setdefault(pair,set()).update((i,j,k))
 lines={frozenset(v) for v in pairs.values()}
 for line in lines:
  for pair in itertools.combinations(sorted(line),2):
   if pairs.get(pair)!=set(line):raise ValueError('inconsistent line incidence')
 hist=Counter(map(len,lines));hist[2]=math.comb(n,2)-sum(math.comb(k,2)*v for k,v in hist.items());hist=dict(sorted(hist.items()))
 assert sum(math.comb(k,2)*v for k,v in hist.items())==math.comb(n,2)
 print(f'{path}: N={n} field_degree={d} histogram={hist} no_five={max(hist)<5} no_six={max(hist)<6}',flush=True)
 return hist
if __name__=='__main__':
 if len(sys.argv)<2:raise SystemExit('usage: python verify_independent_algebraic.py TRUSTED_ENTRY ...')
 for path in sys.argv[1:]:census(path)
