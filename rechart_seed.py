#!/usr/bin/env python3
"""Apply an exact invertible projective rechart to a trusted number-field entry.

Maps (x,y,1) to (x,y,1+u*x+v*y). This moves the old line at infinity into
an affine line, exposing intersections of formerly parallel lines. Does not itself
improve the orchard score. Choose rational u,v so no seed point goes to infinity.
"""
import argparse
from fractions import Fraction
from pathlib import Path
from algfield import parse_entry
from orchard import analyze,entry_text
p=argparse.ArgumentParser();p.add_argument('entry');p.add_argument('output');p.add_argument('--u',default='1/5');p.add_argument('--v',default='1/7');a=p.parse_args()
F,P=parse_entry(Path(a.entry).read_text());u=F.elem(Fraction(a.u));v=F.elem(Fraction(a.v));Q=[]
for x,y in P:
 d=F.add(F.ONE,F.add(F.mul(u,x),F.mul(v,y)))
 if F.is_zero(d):raise ValueError('a seed point maps to infinity; use another chart')
 Q.append((F.div(x,d),F.div(y,d)))
r=analyze(F,P);s=analyze(F,Q)
assert r['valid'] and s['valid'] and r['score']==s['score']
Path(a.output).parent.mkdir(parents=True,exist_ok=True);Path(a.output).write_text(entry_text(F,Q,[int(c) for c in F.p])+'\n')
print(f'{a.entry}: n={len(P)} score={r["score"]} -> {a.output}; u={a.u}, v={a.v}')
