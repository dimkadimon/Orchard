#!/usr/bin/env python3
"""Exhaustively add intersections of 4-point and spanned lines to Du's #0 seed."""
import itertools,json
from pathlib import Path
from investigate_five import census,cross
from beam_orchard import chart
from verify_independent_algebraic import census as independent
P=list(map(tuple,json.loads(Path('research/five_point_23_seed_deletions.json').read_text())['results']['23']['homogeneous_points']))
h,lines=census(P);four=[(l,m) for l,m in lines.items() if m.bit_count()==4];candidates=set()
for l,m in four:
 for k,n in lines.items():
  if m&n:continue
  q=cross(l,k)
  if q is not None and q not in P:candidates.add(q)
valid=[]
for q in candidates:
 hh,_=census(P+[q])
 if max(hh)<=5:valid.append((hh.get(5,0),q,hh))
valid.sort(reverse=True);score,q,hh=valid[0];assert score>=16
Q,cv=chart(P+[q]);path=Path('entries/five_rational_N24.txt');path.write_text('a = 0\n'+', '.join(f'({x}, {y})' for x,y in Q)+'\n')
assert independent(str(path))==hh
out={'seed':'Du 2026 forum post 22 rational #0','exhaustive_scope':'intersections of seed four-point lines with all disjoint seed spanned lines','candidate_count':len(candidates),'valid_count':len(valid),'best_score':score,'histogram':hh,'added_point':q,'homogeneous_points':P+[q],'entry':str(path)}
Path('research/five_rational_N24.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
