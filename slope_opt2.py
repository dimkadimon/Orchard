"""4-axis orchard: four DISTINCT concurrent axes (1,0),(0,1),(p3,q3),(p4,q4).
Transversal parameters (x,y) in Z^2; t3=p3 x+q3 y-o3 in [0,c], t4=p4 x+q4 y-o4 in [0,d].
Asymptotic count = area of polygon. Requires p3,q3,p4,q4 all nonzero and directions distinct.
"""
import json, math, random
from math import gcd
from slope_opt import region_area

def optimise(p3,q3,p4,q4,tries=6000,seed=7):
    rng=random.Random(seed); best=0.0; bestp=None
    D=abs(p3*q4-p4*q3)
    for it in range(tries):
        if bestp is None or it%400==0 or rng.random()<0.15:
            w=[rng.random()+0.05 for _ in range(4)]; s=sum(w)
            a,b,c,d=[x/s for x in w]
            o3=rng.uniform(-(abs(p3)+abs(q3)), abs(p3)+abs(q3))
            o4=rng.uniform(-(abs(p4)+abs(q4)), abs(p4)+abs(q4))
        else:
            a,b,c,d,o3,o4=bestp
            sc=0.06*max(0.02,(1-it/tries))**1.5
            a=max(1e-5,a+rng.gauss(0,sc)); b=max(1e-5,b+rng.gauss(0,sc))
            c=max(1e-5,c+rng.gauss(0,sc)); d=max(1e-5,d+rng.gauss(0,sc))
            f=1.0/(a+b+c+d); a*=f;b*=f;c*=f;d*=f
            o3+=rng.gauss(0,sc*2); o4+=rng.gauss(0,sc*2)
        A=region_area(a,b,c,d,o3,o4,p3,q3,p4,q4)
        if A>best: best=A; bestp=(a,b,c,d,o3,o4)
    # local polish
    for rnd in range(6):
        sc=0.01/(2**rnd)
        for it in range(1500):
            a,b,c,d,o3,o4=bestp
            aa=max(1e-6,a+rng.gauss(0,sc));bb=max(1e-6,b+rng.gauss(0,sc))
            cc=max(1e-6,c+rng.gauss(0,sc));dd=max(1e-6,d+rng.gauss(0,sc))
            f=1.0/(aa+bb+cc+dd);aa*=f;bb*=f;cc*=f;dd*=f
            A=region_area(aa,bb,cc,dd,o3+rng.gauss(0,sc),o4+rng.gauss(0,sc),p3,q3,p4,q4)
            if A>best: best=A;bestp=(aa,bb,cc,dd,o3,o4)
            # also try shifting offsets only
            A2=region_area(a,b,c,d,o3+rng.gauss(0,sc),o4+rng.gauss(0,sc),p3,q3,p4,q4)
            if A2>best: best=A2;bestp=(a,b,c,d,o3,o4)
    return best,bestp

if __name__=='__main__':
    cands=set()
    for p3 in range(-4,5):
        for q3 in range(-4,5):
            if p3==0 or q3==0: continue
            if gcd(abs(p3),abs(q3))!=1: continue
            if (p3,q3)<(0,0): p3n,q3n=-p3,-q3
            else: p3n,q3n=p3,q3
            for p4 in range(-4,5):
                for q4 in range(-4,5):
                    if p4==0 or q4==0: continue
                    if gcd(abs(p4),abs(q4))!=1: continue
                    p4n,q4n=(-p4,-q4) if (p4,q4)<(0,0) else (p4,q4)
                    if (p3n,q3n)==(p4n,q4n): continue
                    if p3n*q4n-p4n*q3n==0: continue
                    key=tuple(sorted([(p3n,q3n),(p4n,q4n)]))
                    cands.add((key[0][0],key[0][1],key[1][0],key[1][1]))
    cands=sorted(cands)
    print("distinct-direction candidates:",len(cands))
    res=[]
    for sl in cands:
        A,pars=optimise(*sl,tries=2500)
        res.append((A,sl,pars,abs(sl[0]*sl[3]-sl[2]*sl[1])))
    res.sort(key=lambda r:-r[0])
    print("\nTOP 20 (Xiang Li = 1/28 = 0.0357143)")
    for A,sl,pars,D in res[:20]:
        print("%.7f = 1/%6.3f  det=%d  slopes (1,0),(0,1),%s,%s" % (A,1/A,D,(sl[0],sl[1]),(sl[2],sl[3])))
    print("\nXiang Li entry:")
    for A,sl,pars,D in res:
        if set([(sl[0],sl[1]),(sl[2],sl[3])])=={(1,1),(1,-1)}:
            print("  %.7f  pars=%s"%(A,[round(x,5) for x in pars]))
    json.dump([{'area':A,'slopes':sl,'pars':pars,'det':D} for A,sl,pars,D in res],open('slope_opt2_results.json','w'),indent=1)
