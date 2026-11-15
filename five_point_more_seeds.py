#!/usr/bin/env python3
"""Exact census and exhaustive selected-size subsets of a published 23-point seed."""
from fractions import Fraction as F
from itertools import combinations
from pathlib import Path
import json
from investigate_five import normalize,census

# Zhao Hui Du, https://bbs.emath.ac.cn/thread-50370-3-1.html, post 22, rational #0.
P=[(0,0,1),(0,1,1),(1,0,0),(F(3,5),F(1,3),1),(1,0,1),(1,F(-5,6),0),
   (F(1,5),F(1,2),1),(F(4,5),F(1,3),1),(1,F(-5,8),0),(F(2,5),F(5,12),1),
   (F(8,5),0,1),(F(2,5),F(3,8),1),(F(4,5),F(1,4),1),(F(3,5),F(1,4),1),
   (0,F(1,2),1),(1,F(-5,12),0),(F(2,5),F(1,3),1),(F(4,5),F(1,6),1),
   (F(2,5),F(1,2),1),(F(4,5),0,1),(F(8,15),F(1,3),1),(0,1,0),(0,F(2,3),1)]
P=list(map(normalize,P));h,lines=census(P);assert h.get(5)==15 and max(h)==5
rich=[m for m in lines.values() if m.bit_count()==5]
out={}
for n in (15,18,19,20,21,22,23):
 best=-1;ids_best=None
 for ids in combinations(range(len(P)),n):
  mask=sum(1<<i for i in ids);score=sum(m&mask==m for m in rich)
  if score>best:best,ids_best=score,ids
 subset=[P[i] for i in ids_best]; hist,_=census(subset)
 assert hist.get(5,0)==best and max(hist)<=5
 out[str(n)]={'score':best,'indices':ids_best,'histogram':hist,'homogeneous_points':subset}
 print(f'23-point seed, n={n}: {best} strict lines; histogram {hist}',flush=True)
Path('research/five_point_23_seed_deletions.json').write_text(json.dumps({'source':'https://bbs.emath.ac.cn/thread-50370-3-1.html post 22 rational #0','method':'exhaustive subsets of this one seed at listed sizes','seed_histogram':h,'results':out},indent=2)+'\n')
