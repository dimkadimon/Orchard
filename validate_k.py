"""Validate the general-k second-order inequality on real configurations (grids,
pencils, random sets) for k=4..8. Must never be violated."""
import numpy as np, itertools, math, random
from math import comb
from collections import defaultdict
def w(j):
    if j<=2: return 0
    if j==3: return 3
    return comb(j,2)+(j-3)
def census(P):
    n=len(P); lines=defaultdict(set)
    for i,j in itertools.combinations(range(n),2):
        x1,y1=P[i];x2,y2=P[j]
        A=y1-y2;B=x2-x1;C=x1*y2-x2*y1
        g=math.gcd(math.gcd(abs(A),abs(B)),abs(C)) or 1
        A//=g;B//=g;C//=g
        if A<0 or (A==0 and B<0):A,B,C=-A,-B,-C
        lines[(A,B,C)].update((i,j))
    return lines
def test(P,k,tag):
    n=len(P); lines=census(P)
    rich=[v for v in lines.values() if len(v)>=k]
    L=len(rich)
    if L<3: return None
    r=[0]*n
    for v in rich:
        for i in v: r[i]+=1
    inter=set.intersection(*[set(v) for v in rich])
    pencil=len(inter)==1
    lhs=sum(w(x) for x in r); rhs=comb(L,2)-3
    ok=(lhs<=rhs) or pencil
    return (tag,n,k,L,lhs,rhs,ok,pencil)
sets=[]
for a in range(2,9):
    for b in range(2,9):
        if 4<=a*b<=90: sets.append((f'grid{a}x{b}',[(x,y) for x in range(a) for y in range(b)]))
# pencils: m lines through origin with t points each
for m in range(3,7):
    for t in range(3,7):
        P=[(0,0)]
        for i in range(m):
            ang=i*math.pi/m
            for s in range(1,t+1):
                P.append((round(s*math.cos(ang)*1e6),round(s*math.sin(ang)*1e6)))
        if len(P)<=90: sets.append((f'pencil{m}x{t}',P))
rng=random.Random(1)
for trial in range(40):
    n=rng.randint(8,40)
    P=list({(rng.randint(-8,8),rng.randint(-8,8)) for _ in range(n*3)})[:n]
    sets.append((f'rand{trial}',P))
bad=[];cnt=0
for tag,P in sets:
    for k in range(4,9):
        r=test(P,k,tag)
        if r:
            cnt+=1
            if not r[6]: bad.append(r)
print("checks:",cnt,"VIOLATIONS:",len(bad))
for b in bad[:10]: print("  ",b)
