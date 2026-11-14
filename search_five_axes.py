"""Search concurrent-axis interval constructions for T4(60).
All scoring is exact integer arithmetic. Deterministic seeded simulated annealing;
not an exhaustive search and no optimality claim is made.
"""
import random, math, json, time, sys
from collections import Counter
from math import gcd

# Point on axis (p,q) is (p/D,q/D), D=p*180+q*60+t.
# We choose directions with positive base denominator and intervals avoiding zero.
def census(axes,lengths,offsets, details=False):
    lines={}
    for i in range(len(axes)):
      p,q=axes[i]; off=offsets[i]; li=lengths[i]; base=p*180+q*60
      for j in range(i+1,len(axes)):
       s,t=axes[j]; bj=s*180+t*60; ofj=offsets[j]; lj=lengths[j]
       C=p*t-q*s
       if C==0: continue
       for a in range(off,off+li):
        D=base+a
        if D<=0:return (-10000,None) if details else -10000
        for b in range(ofj,ofj+lj):
         E=bj+b
         if E<=0:return (-10000,None) if details else -10000
         A=q*E-D*t; B=D*s-p*E
         g=gcd(gcd(abs(A),abs(B)),abs(C)); A//=g;B//=g; CC=C//g
         if A<0 or (A==0 and B<0):A=-A;B=-B;CC=-CC
         key=(A,B,CC); lines[key]=lines.get(key,0)|(1<<i)|(1<<j)
    counts=Counter(v.bit_count() for v in lines.values())
    score=counts[4]+counts[5]+sum(x>=4 for x in lengths)
    return (score,dict(counts)) if details else score

def main():
    n=60; rng=random.Random(314159); axes=[(1,0),(0,1),(1,1),(1,-1),(1,2)]
    # offsets of first 2 can be fixed by shifts of line parameters (translation).
    best=(-1,None); start=time.time(); restarts=int(sys.argv[1]) if len(sys.argv)>1 else 20
    for restart in range(restarts):
      lengths=[13,13,17,13,4] if restart%3==0 else [rng.randint(7,16) for _ in axes]
      while sum(lengths)!=n:
       i=rng.randrange(5); d=1 if sum(lengths)<n else -1
       if lengths[i]+d>=4:lengths[i]+=d
      offsets=[0,0,4,-6,10] if restart%3==0 else [0,0,rng.randint(-10,10),rng.randint(-10,10),rng.randint(-10,30)]
      val=census(axes,lengths,offsets)
      for step in range(2600):
       lens=lengths.copy(); offs=offsets.copy()
       if rng.random()<.55:
        i,j=rng.sample(range(5),2)
        if lens[i]<=4:continue
        lens[i]-=1;lens[j]+=1
       else:
        i=rng.randrange(2,5);offs[i]+=rng.choice([-2,-1,1,2])
       val2=census(axes,lens,offs)
       temp=max(.2,2.0*(1-(step%650)/650))
       if val2>=val or rng.random()<math.exp((val2-val)/temp): lengths,offsets,val=lens,offs,val2
       if val>best[0]:
        best=(val,(lengths.copy(),offsets.copy(),census(axes,lengths,offsets,True)[1]))
        print('best',best,'restart',restart,'step',step,'secs',round(time.time()-start,1),flush=True)
    data={'axes':axes,'n':n,'best':best,'seed':314159,'restarts':restarts,'elapsed':time.time()-start}
    with open('five_axes_search_results.json','w') as f:json.dump(data,f,indent=2)
    print('FINAL',data,flush=True)
if __name__=='__main__':main()
