#!/usr/bin/env python3
"""Exhaustive deletion neighbourhood of a strict algebraic configuration.

Scores are exact using the seed's exact four-line incidence lists. Exhaustive only
for the chosen source and deletion cardinalities; does not imply global optimality.
"""
import argparse,itertools,json,time
from pathlib import Path
from algfield import parse_entry
from orchard import analyze,entry_text

def main():
 p=argparse.ArgumentParser();p.add_argument('--entry',required=True);p.add_argument('--max-delete',type=int,default=7);p.add_argument('--result',required=True);a=p.parse_args()
 F,pts=parse_entry(Path(a.entry).read_text());n=len(pts);res=analyze(F,pts)
 if not res['valid']:raise ValueError('source has five collinear')
 lines=res['sizes'][4];masks=[0]*n
 for j,v in enumerate(lines):
  for i in v:masks[i]|=1<<j
 base={x['n']:x['score'] for x in json.loads(Path('records_int_now.json').read_text())['records']}
 for row in json.loads(Path('final_rows.json').read_text()):base[row['N']]=max(base.get(row['N'],0),row['score'])
 for size in range(42,47):base[size]=max(base.get(size,0),3*size-51)
 records=[];start=time.time()
 for d in range(1,min(a.max_delete,n)+1):
  best=-1;arg=None;number=0
  for removal in itertools.combinations(range(n),d):
   mask=0
   for idx in removal:mask|=masks[idx]
   score=len(lines)-mask.bit_count();number+=1
   if score>best:best,arg=score,removal
  target=n-d;baseline=base.get(target);out={'source':a.entry,'n':n,'deleted':d,'target':target,'sets':number,'best':best,'baseline':baseline,'removed':arg,'seconds':round(time.time()-start,2),'saved':None}
  if baseline is not None and best>baseline:
   Path('research/new_candidates').mkdir(parents=True,exist_ok=True)
   name=f'delete_n{target}_s{best}_from{n}.txt';removed=set(arg)
   Path('research/new_candidates',name).write_text(entry_text(F,[v for i,v in enumerate(pts) if i not in removed],[int(v) for v in F.p])+'\n');out['saved']=name;print('NEW',name,flush=True)
  records.append(out);print(json.dumps(out),flush=True);Path(a.result).write_text(json.dumps(records,indent=2)+'\n')
 print('DONE',flush=True)
if __name__=='__main__':main()
