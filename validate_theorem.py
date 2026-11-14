"""Sanity check: the second-order inequality sum_p w(r_p) <= C(L,2)-3 must hold
for EVERY real configuration. Test on all 50 integer board records + our 50 algebraic."""
import json, itertools, math
from math import comb
from collections import defaultdict
import numpy as np
from algfield import parse_entry
from orchard import float_lines
from make_figures import numeric_points

def w(j):
    if j<=2: return 0
    if j==3: return 3
    return comb(j,2)+(j-3)

def check(P,tag):
    n=len(P); P=np.asarray(P,float)
    scale=max(1.0,np.abs(P).max())
    lines={}
    for i in range(n):
        for j in range(i+1,n):
            d=P[j]-P[i]; Ln=np.hypot(*d)
            if Ln<1e-12: continue
            cr=np.abs((P[:,0]-P[i,0])*d[1]-(P[:,1]-P[i,1])*d[0])/Ln
            on=tuple(np.where(cr<1e-7*scale)[0])
            if len(on)>=4: lines[on]=len(on)
    L=len(lines)
    if L<3: return None
    r=np.zeros(n,int)
    for q in lines:
        for i in q: r[i]+=1
    # pencil?
    concurrent = len(set(lines.keys()))>0 and all(len(set(q)&set(q2))==1 for q in lines for q2 in lines) and len(set.intersection(*[set(q) for q in lines]))==1
    lhs=int(sum(w(int(x)) for x in r)); rhs=comb(L,2)-3
    ok = lhs<=rhs or concurrent
    return (tag,n,L,lhs,rhs,ok,bool(concurrent),int(max(r)),(n-1)//3)

rows=[]
ints=json.load(open('integer_solutions.json'))
for k,v in ints.items():
    pts=v['points'] if isinstance(v,dict) else v
    P=[(float(a),float(b)) for a,b in pts]
    rows.append(check(P,f'int N={k}'))
for N in range(11,61):
    try:
        F,pts=parse_entry(open(f'entries/N{N}.txt').read())
        best=None
        for rt,P in numeric_points(F,pts):
            if len(set(map(tuple,np.round(P,9))))<len(P): continue
            sz=float_lines(P,tol=1e-7)
            if any(s>=5 for s in sz): continue
            sc=len(sz.get(4,[]))
            if best is None or sc>best[0]: best=(sc,P)
        if best: rows.append(check(best[1],f'alg N={N}'))
    except Exception as e:
        print('skip',N,e)
bad=[r for r in rows if r and not r[5]]
print(f"checked {len([r for r in rows if r])} configurations")
print("VIOLATIONS:",len(bad))
for b in bad: print("  ",b)
print("\nsample (tag,n,L,sum_w,C(L,2)-3,ok,pencil,maxdeg,cap):")
for r in rows[:6]: print("  ",r)
for r in rows:
    if r and r[1] in (11,13,28,36,48,60): print("  ",r)
json.dump([r for r in rows if r],open('validate_theorem_results.json','w'),indent=1)
