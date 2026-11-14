"""Optimize the four-axis (Xiang Li-type) construction for T_4(n):
lines with AT LEAST 4 points.  Points lie on four concurrent axes; a transversal
contains one point from each axis iff x+y in C and x-y in D.  Long axes count.
We search interval sets A,B,C,D and shifts.  Exact rational scoring for best.
"""
import json, time, itertools, sys
from fractions import Fraction
T0=time.time(); BUDGET=float(sys.argv[1]) if len(sys.argv)>1 else 900.0
OUT='atleast4_li_results.json'
REC_INT={r['n']:r['score'] for r in json.load(open('records_int_now.json'))['records']}
REC_ALG={e['N']:e['score'] for e in json.load(open('final_rows.json'))}

def count_rect(a,b,c,d):
    # grid of (s=x+y, t=x-y) for x in [0,a), y in [0,b)
    S=a+b-1; T=a+b-1  # t offset b-1 -> 0..a+b-2
    grid=[[0]*T for _ in range(S)]
    for x in range(a):
        for y in range(b):
            grid[x+y][x-y+b-1]=1
    pref=[[0]*(T+1) for _ in range(S+1)]
    for i in range(S):
        row=grid[i]; pr=pref[i]; pn=pref[i+1]
        s=0
        for j in range(T):
            s+=row[j]; pn[j+1]=pr[j+1]+s
    best=0; bu=bv=0
    if c>S or d>T: return 0,0,0
    for u in range(S-c+1):
        p1=pref[u]; p2=pref[u+c]
        for v in range(T-d+1):
            val=p2[v+d]-p2[v]-p1[v+d]+p1[v]
            if val>best:
                best=val; bu=u; bv=v
    return best,bu,bv-(b-1)

def li_default(n):
    s=3*(n+1)//7
    a=s//2; b=(s+1)//2
    c=(n-s)//2; d=(n-s+1)//2
    u=(s-c-1)//2; v=(s-d)//2  # ceil((s-d-1)/2)
    # D start = -(b-1)+v
    return a,b,c,d,u,-(b-1)+v

def score_sets(n,a,b,c,d,u,v):
    A=[3*n+x for x in range(a)]
    B=[n+x for x in range(b)]
    C=[4*n+u+x for x in range(c)]
    D=[2*n+v+x for x in range(d)]
    if min(A+B+C+D)<=0: return -1
    pts=[(Fraction(1,al),Fraction(0)) for al in A]
    pts+=[(Fraction(0),Fraction(1,be)) for be in B]
    pts+=[(Fraction(1,ga),Fraction(1,ga)) for ga in C]
    pts+=[(Fraction(1,de),Fraction(-1,de)) for de in D]
    # exact line count
    lines={}
    for i,j in itertools.combinations(range(n),2):
        x1,y1=pts[i]; x2,y2=pts[j]
        Acoef=y1-y2; Bcoef=x2-x1; Ccoef=x1*y2-x2*y1
        # normalize by first nonzero
        for c0 in (Acoef,Bcoef,Ccoef):
            if c0!=0:
                Acoef/=c0; Bcoef/=c0; Ccoef/=c0; break
        lines.setdefault((Acoef,Bcoef,Ccoef),set()).update((i,j))
    t4=sum(1 for vset in lines.values() if len(vset)==4)
    T4=sum(1 for vset in lines.values() if len(vset)>=4)
    return T4

results={}
for n in range(11,61):
    if time.time()-T0>BUDGET: break
    best=0; bestpars=None
    # default Li
    a,b,c,d,u,v=li_default(n)
    if min(a,b,c,d)>=1:
        cnt,_,_=count_rect(a,b,c,d)
        T=cnt+sum(1 for z in (a,b,c,d) if z>=4)
        if T>best: best=T; bestpars=(a,b,c,d,u,v,'default')
    # search compositions
    for a in range(1,n-2):
        if time.time()-T0>BUDGET: break
        for b in range(1,n-a-1):
            rem=n-a-b
            if rem<2: continue
            # quick upper bound
            if a*b + 4 <= best: continue
            for c in range(1,rem):
                d=rem-c
                if d<1: continue
                if min(a*b,c*d)+4 <= best: continue
                cnt,u,v=count_rect(a,b,c,d)
                T=cnt+sum(1 for z in (a,b,c,d) if z>=4)
                if T>best:
                    best=T; bestpars=(a,b,c,d,u,v,'search')
    # verify exact
    if bestpars:
        a,b,c,d,u,v,tag=bestpars
        exact=score_sets(n,a,b,c,d,u,v)
        if exact!=best:
            print('mismatch',n,best,exact,bestpars,flush=True)
            best=exact
    prev=max(REC_INT.get(n,0),REC_ALG.get(n,0))
    results[n]={'T4':best,'pars':bestpars,'prev_best':prev,'improves':best>prev}
    print(f'n={n} T4={best} prev={prev} pars={bestpars} t={time.time()-T0:.0f}s',flush=True)
    json.dump(results,open(OUT,'w'),indent=1)
json.dump(results,open(OUT,'w'),indent=1)
print('DONE',flush=True)
