"""Exact deletion search from certified configurations; exhaustive only for the
specified source sets and deletion sizes. Lower-bound comparison uses the pinned
record files; the script emits candidate coordinates only on improvement.
"""
from itertools import combinations
from algfield import parse_entry
from orchard import analyze,entry_text
import json,time,sys
records={x['n']:x['score'] for x in json.load(open('records_int_now.json'))['records']}
for x in json.load(open('final_rows.json')):records[x['N']]=max(records.get(x['N'],0),x['score'])
MAXD=int(sys.argv[1]) if len(sys.argv)>1 else 5
results=[]; start=time.time()
for n in [28,36,40,48]:
 F,pts=parse_entry(open(f'entries/N{n}.txt').read());res=analyze(F,pts)
 fours=res['sizes'][4]; masks=[0]*n
 for i,line in enumerate(fours):
  for v in line:masks[v]|=1<<i
 for d in range(1,MAXD+1):
  best=-1; arg=None; count=0
  for remove in combinations(range(n),d):
   mask=0
   for i in remove:mask|=masks[i]
   score=len(fours)-mask.bit_count();count+=1
   if score>best:best=score;arg=remove
  prev=records.get(n-d,0); r={'source':n,'removed':d,'size':n-d,'score':best,'prior':prev,'improved':best>prev,'sets':count,'indices':arg};results.append(r)
  print(r,'seconds',round(time.time()-start,1),flush=True)
  if best>prev:
   remaining=[p for i,p in enumerate(pts) if i not in arg]
   poly=[int(x) for x in F.p];fn=f'entries/improved_N{n-d}_from_{n}.txt'
   open(fn,'w').write(entry_text(F,remaining,poly)+'\n')
   records[n-d]=best
with open('deletion_search_results.json','w') as f:json.dump(results,f,indent=2)
