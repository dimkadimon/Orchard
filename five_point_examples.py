#!/usr/bin/env python3
"""Exact projective certificates of five-point constructions, no six collinear."""
from fractions import Fraction as F
from itertools import combinations,permutations
from pathlib import Path
import json
from investigate_five import normalize,census
from verify_five_independent import census as determinant_census

# Intersections of y=i*x+i^2, i=0,...,5; a direct generic six-line construction.
P15=[(-(i+j),-i*j,1) for i,j in combinations(range(6),2)]
# Du, post 9, https://bbs.emath.ac.cn/thread-50370-1-1.html
P18=[(-68,11,-257),(37,-10,205),(-269,62,-1241),(8,-3,45),(373,-101,2137),
     (600,-147,2733),(-29,23,-191),(14,-5,83),(222,-66,1356),(-28,7,-133),
     (167,-44,923),(9,-3,51),(-240,66,-1374),(112,-22,490),(-217,49,-973),
     (-5,2,-29),(147,-42,903),(-80,40,-440)]
# Du, post 10, same thread. Earlier claims in that thread were corrected.
P19=[(1,0,1),(F(-1,4),0,1),(F(-1,2),1,0),(0,4,1),(F(-1,2),-1,1),
     (F(1,2),-1,1),(F(3,4),1,1),(F(3,2),1,1),(F(1,2),3,1),(F(1,2),1,0),
     (1,2,1),(0,-2,1),(0,1,1),(F(1,4),1,0),(1,0,0),(0,0,1),
     (F(1,2),0,1),(F(1,2),1,1),(0,1,0)]
# Elkies' 25-point example in barycentric/projective coordinates, as recorded
# by Pegg: https://blog.wolfram.com/2018/02/02/cultivating-new-solutions-for-the-orchard-planting-problem/
P25=sorted({p for t in [(1,1,1),(0,0,1),(1,2,3),(1,1,2),(0,1,1),(1,2,2),(0,1,2)] for p in permutations(t)})
P24=[p for i,p in enumerate(P25) if i!=3]
examples={}
for n,P,score,source in [(15,P15,6,'six parabolic tangents'),
  (18,P18,9,'https://bbs.emath.ac.cn/thread-50370-1-1.html post 9'),
  (19,P19,10,'https://bbs.emath.ac.cn/thread-50370-1-1.html post 10'),
  (24,P24,15,'one-point deletion from Elkies/Pegg 25-point example'),
  (25,P25,18,'https://blog.wolfram.com/2018/02/02/cultivating-new-solutions-for-the-orchard-planting-problem/')]:
 P=list(map(normalize,P));h,_=census(P);h2=determinant_census(P)
 assert len(P)==n and h==h2 and h.get(5)==score and max(h)==5
 examples[str(n)]={'score':score,'source':source,'histogram':h,'homogeneous_points':P}
 print(f'n={n} strict score={score}, histogram={h}; two independent exact censuses agree')
Path('research/five_point_examples.json').write_text(json.dumps(examples,indent=2)+'\n')
