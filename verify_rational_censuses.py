"""Independent exact integer/rational line-census checks for the revised paper."""
from fractions import Fraction as F
from itertools import combinations
from collections import Counter
from functools import reduce
from math import gcd,lcm

def line(p,q):
 x,y=p; u,v=q
 a,b,c=y-v,u-x,x*v-u*y
 denom=lcm(a.denominator,b.denominator,c.denominator)
 k=[int(t*denom) for t in (a,b,c)]
 g=reduce(gcd,k)
 k=[t//g for t in k]
 if next(t for t in k if t)<0:k=[-t for t in k]
 return tuple(k)

def census(P):
 assert len(set(P))==len(P)
 D={}
 for i,j in combinations(range(len(P)),2):D.setdefault(line(P[i],P[j]),set()).update((i,j))
 hist=dict(sorted(Counter(map(len,D.values())).items()))
 assert sum(m*(m-1)//2*n for m,n in hist.items())==len(P)*(len(P)-1)//2
 return hist

P12=[(0,0),(42,84),(56,112),(70,70),(84,28),(84,63),(84,84),(84,168),(96,72),(105,105),(112,84),(126,42)]
P12=[(F(x),F(y)) for x,y in P12]
P60=[(F(1,180+x),F(0)) for x in range(13)]+[(F(0),F(1,60+y)) for y in range(13)]+[(F(1,240+z),F(1,240+z)) for z in range(4,21)]+[(F(1,120+w),-F(1,120+w)) for w in range(-8,9)]
a=census(P12);b=census(P60)
print('n=12:',a);print('n=60:',b)
assert a.get(4)==7 and max(a)==4
assert b=={2:400,3:56,4:129,13:2,17:2}
