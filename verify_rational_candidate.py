#!/usr/bin/env python3
"""Independent Fraction-based census for rational beam/MILP candidate files.

Candidate format: first line 'a = 0', followed by comma-separated integer coordinate
pairs. This intentionally does not use the projective beam or number-field arithmetic.
"""
import ast, sys
from fractions import Fraction
from itertools import combinations
from collections import Counter
from math import gcd, lcm
from functools import reduce

def census(path):
    text=open(path).read().splitlines()
    if text[0].strip() != 'a = 0':raise ValueError('expected rational a = 0 candidate')
    pts=ast.literal_eval('['+' '.join(text[1:])+']')
    pts=[(Fraction(x),Fraction(y)) for x,y in pts]
    if len(set(pts)) != len(pts):raise ValueError('duplicate points')
    lines={}
    for i,j in combinations(range(len(pts)),2):
        x,y=pts[i];u,v=pts[j]
        a,b,c=y-v,u-x,x*v-u*y
        d=lcm(a.denominator,b.denominator,c.denominator)
        A=[int(t*d) for t in (a,b,c)]; g=reduce(gcd,A);A=[k//g for k in A]
        if next(k for k in A if k)<0:A=[-k for k in A]
        lines.setdefault(tuple(A),set()).update((i,j))
    hist=dict(sorted(Counter(len(s) for s in lines.values()).items()))
    assert sum(k*(k-1)//2*v for k,v in hist.items())==len(pts)*(len(pts)-1)//2
    print(path, 'n=',len(pts),'hist=',hist,'strict=',max(hist,default=0)<5)
    if max(hist,default=0)>=5:raise ValueError('five or more collinear')
if __name__=='__main__':
    if len(sys.argv)<2:raise SystemExit('usage: python verify_rational_candidate.py FILE ...')
    for p in sys.argv[1:]:census(p)
