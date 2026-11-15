#!/usr/bin/env python3
"""Verify published 5-orchards and search all deletions of a published seed.

Projective coordinates transcribed from Zhao Hui Du's 2026 forum post:
https://bbs.emath.ac.cn/thread-50370-2-1.html (post #11, third solution).
No claim of novelty is made. All arithmetic uses integers/Fractions.
"""
from fractions import Fraction
from itertools import combinations
from collections import Counter
from math import gcd,lcm
from pathlib import Path
import json
from second_order_k import bound,melchior_k

def normalize(v):
    vals=[Fraction(x) for x in v];d=lcm(*(x.denominator for x in vals))
    t=[int(x*d) for x in vals];g=gcd(*t)
    if not g:raise ValueError('zero vector')
    t=tuple(x//g for x in t)
    return tuple(-x for x in t) if next(x for x in t if x)<0 else t

def cross(p,q):
    a,b,c=p;d,e,f=q
    return normalize((b*f-c*e,c*d-a*f,a*e-b*d))

def census(P):
    P=list(map(normalize,P));assert len(set(P))==len(P)
    lines={}
    for i,j in combinations(range(len(P)),2):
        l=cross(P[i],P[j]);lines[l]=lines.get(l,0)|(1<<i)|(1<<j)
    h=Counter(x.bit_count() for x in lines.values())
    assert sum(k*(k-1)//2*v for k,v in h.items())==len(P)*(len(P)-1)//2
    return dict(sorted(h.items())),lines

F=Fraction
P20=[(0,0,1),(1,-3,0),(2,-1,1),(1,2,0),(1,3,1),(2,9,1),
     (0,1,1),(1,0,0),(1,-2,0),(1,1,1),(1,0,1),(0,3,1),
     (2,-3,1),(1,-1,1),(2,1,1),(0,-3,1),(F(3,2),0,1),
     (F(1,2),0,1),(F(2,3),1,1),(0,1,0)]
def main():
    global P20
    P20=list(map(normalize,P20));h,lines=census(P20)
    assert h.get(5)==11 and max(h)==5
    rich=[m for m in lines.values() if m.bit_count()==5]
    results={}
    for n in range(9,21):
        best=-1;bestids=None
        for ids in combinations(range(20),n):
            mask=sum(1<<i for i in ids);score=sum((m&mask)==m for m in rich)
            if score>best:best,bestids=score,ids
        subset=[P20[i] for i in bestids];hist,_=census(subset)
        assert hist.get(5,0)==best and max(hist)<=5
        results[str(n)]={'score':best,'indices':bestids,'histogram':hist,'classical_upper':melchior_k(n,5),'second_order_upper':bound(n,5),'homogeneous_points':subset}
        print(f'n={n:2} strict lower={best:2} second-order upper={bound(n,5):2} hist={hist}')
    Path('research/five_point_seed_deletions.json').write_text(json.dumps({'source':'https://bbs.emath.ac.cn/thread-50370-2-1.html post 11 third solution','method':'exhaustive subsets of this one 20-point seed; not a global search','seed_histogram':h,'results':results},indent=2)+'\n')

if __name__ == "__main__":
    main()
