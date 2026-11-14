"""Optimise 4-axis orchard constructions over the choice of axis slopes.
Points on axis i at (1/(c_i+t))*(p_i,q_i), t in interval T_i.
A transversal uX+vY=1 meets axis i at t = u p_i + v q_i - c_i.
With axes (1,0),(0,1),(p3,q3),(p4,q4) and x=t1,y=t2:
   0<=x<=a, 0<=y<=b, 0<=p3 x+q3 y-o3<=c, 0<=p4 x+q4 y-o4<=d.
Asymptotic transversal count = area of this convex polygon; maximise area/n^2
subject to a+b+c+d=n.  Xiang Li = (1,1),(1,-1).
"""
import itertools, json, math, random

def clip(poly, A, B, C):
    # keep A*x+B*y <= C
    out=[]
    n=len(poly)
    for i in range(n):
        P=poly[i]; Q=poly[(i+1)%n]
        fp=A*P[0]+B*P[1]-C; fq=A*Q[0]+B*Q[1]-C
        if fp<=1e-12: out.append(P)
        if (fp<0 and fq>0) or (fp>0 and fq<0):
            t=fp/(fp-fq)
            out.append((P[0]+t*(Q[0]-P[0]), P[1]+t*(Q[1]-P[1])))
    return out

def area(poly):
    if len(poly)<3: return 0.0
    s=0.0
    for i in range(len(poly)):
        x1,y1=poly[i]; x2,y2=poly[(i+1)%len(poly)]
        s+=x1*y2-x2*y1
    return abs(s)/2

def region_area(a,b,c,d,o3,o4,p3,q3,p4,q4):
    poly=[(0,0),(a,0),(a,b),(0,b)]
    for A,B,C in [(1,0,a),(-1,0,0),(0,1,b),(0,-1,0),
                  (p3,q3,o3+c),(-p3,-q3,-o3),
                  (p4,q4,o4+d),(-p4,-q4,-o4)]:
        poly=clip(poly,A,B,C)
        if len(poly)<3: return 0.0
    return area(poly)

def optimise(p3,q3,p4,q4,n=1.0,tries=4000,seed=0):
    rng=random.Random(seed)
    best=0; bestp=None
    # variables: a,b,c (d = n-a-b-c), o3,o4
    for it in range(tries):
        if it<tries//2:
            w=[rng.random() for _ in range(3)]; s=sum(w)+1e-9
            a,b,c=[n*x/(s+1) for x in w]
            d=n-a-b-c
            o3=rng.uniform(-max(p3,q3,0)*a-max(p3,q3,0)*b, (abs(p3)+abs(q3))*max(a,b))
            o4=rng.uniform(-max(p4,q4,0)*a-max(p4,q4,0)*b, (abs(p4)+abs(q4))*max(a,b))
        else:
            if bestp is None: break
            a,b,c,d,o3,o4=bestp
            sc=0.08*n*max(0.05,1-it/tries)
            a=max(1e-6,a+rng.gauss(0,sc)); b=max(1e-6,b+rng.gauss(0,sc))
            c=max(1e-6,c+rng.gauss(0,sc))
            if a+b+c>=n: 
                f=n/(a+b+c+1e-9)*0.999; a*=f;b*=f;c*=f
            d=n-a-b-c
            o3+=rng.gauss(0,sc); o4+=rng.gauss(0,sc)
        if d<=0: continue
        A=region_area(a,b,c,d,o3,o4,p3,q3,p4,q4)
        if A>best: best=A; bestp=(a,b,c,d,o3,o4)
    return best,bestp

if __name__=='__main__':
    cands=[]
    for p3 in range(-3,4):
        for q3 in range(-3,4):
            if (p3,q3)==(0,0): continue
            for p4 in range(-3,4):
                for q4 in range(-3,4):
                    if (p4,q4)==(0,0): continue
                    if p3*q4-p4*q3==0: continue
                    # normalise: axis direction sign irrelevant
                    if (p3,q3)<(0,0): p3,q3=-p3,-q3
                    if (p4,q4)<(0,0): p4,q4=-p4,-q4
                    cands.append((p3,q3,p4,q4))
    cands=sorted(set(cands))
    print("candidate slope pairs:",len(cands))
    res=[]
    for (p3,q3,p4,q4) in cands:
        A,pars=optimise(p3,q3,p4,q4,tries=1200,seed=1)
        res.append((A,(p3,q3,p4,q4),pars))
    res.sort(key=lambda r:-r[0])
    print("TOP 15 by area/n^2  (Xiang Li (1,1),(1,-1) should be ~1/28=0.035714)")
    for A,sl,pars in res[:15]:
        print("%.6f  1/%.3f  slopes %s  pars %s"%(A,1/A if A>0 else 0,sl,[round(x,4) for x in pars] if pars else None))
    for A,sl,pars in res:
        if sl==(1,1,1,-1) or sl==(1,-1,1,1):
            print("XiangLi ref:",round(A,6),sl)
    json.dump([{'area':A,'slopes':sl,'pars':pars} for A,sl,pars in res],open('slope_opt_results.json','w'),indent=1)
