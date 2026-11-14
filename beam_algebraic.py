#!/usr/bin/env python3
"""Algebraic-field beam search using floating-point proposal generation.

Start from exact algebraic configurations, delete several points, and beam-regrow at
intersections of 2- and 3-point lines. Numerical incidence grouping is used only as a
heuristic. Each surviving proposed point is constructed exactly in Q(a), and any
candidate exceeding the record is re-censused exactly before it is written. This
bounded search is incomplete and cannot prove an upper bound.
"""
import argparse,itertools,json,random,time
from pathlib import Path
import numpy as np
import sympy as sp
from algfield import parse_entry
from orchard import analyze,line_key,entry_text
from extend import exact_point_from_lines,to_float_pts

def floatscore(P,tol):
    n=len(P);scale=max(1.,float(np.max(np.abs(P)))); groups={}
    for i in range(n):
      d=P[i+1:]-P[i];norm=np.linalg.norm(d,axis=1)
      for jj in range(len(d)):
       j=i+1+jj
       if norm[jj]<tol*scale:continue
       v=d[jj]; dist=np.abs((P[:,0]-P[i,0])*v[1]-(P[:,1]-P[i,1])*v[0])/norm[jj]
       ids=tuple(np.flatnonzero(dist<tol*scale))
       if len(ids)>=2 and ids not in groups:
        a,b=v[1]/norm[jj],-v[0]/norm[jj];c=-a*P[i,0]-b*P[i,1]
        groups[ids]=((a,b,c),(i,j),set(ids))
    lines={k:[] for k in (2,3,4,5)}
    for ids,(abc,ij,idx) in groups.items():lines[min(len(ids),5)].append((abc,ij,idx))
    return lines

def proposals(P,lines,rng,pairs32,pairs22,tol):
    two,three,four=lines[2],lines[3],lines[4]+lines[5]
    pairs=[]
    for A,B in itertools.combinations(three,2):
     if A[2].isdisjoint(B[2]):pairs.append((A,B))
    if two and three:
     for _ in range(pairs32):
      A=rng.choice(three);B=rng.choice(two)
      if A[2].isdisjoint(B[2]):pairs.append((A,B))
    if len(two)>1:
     for _ in range(pairs22):
      A,B=rng.sample(two,2)
      if A[2].isdisjoint(B[2]):pairs.append((A,B))
    scale=max(1.,float(np.max(np.abs(P)))); dedupe=set();out=[]
    for A,B in pairs:
      (a,b,c)=A[0];(d,e,f)=B[0]; det=a*e-b*d
      if abs(det)<1e-9:continue
      x=(b*f-c*e)/det;y=(c*d-a*f)/det
      if max(abs(x),abs(y))>1e4*scale:continue
      p=np.array([x,y]);key=(round(x/(scale*tol)),round(y/(scale*tol)))
      if key in dedupe:continue
      dedupe.add(key)
      if np.min(np.linalg.norm(P-p,axis=1))<tol*scale:continue
      if any(abs(a*x+b*y+c)<tol*scale for (a,b,c),_,_ in four):continue
      gain=sum(abs(a*x+b*y+c)<tol*scale for (a,b,c),_,_ in three)
      trips=sum(abs(a*x+b*y+c)<tol*scale for (a,b,c),_,_ in two)
      out.append((gain,trips,A[1],B[1],p))
    return out

def run(n,depth,width,rounds,seednum,pairs32,pairs22,tol,result):
    rng=random.Random(seednum); F,P0=parse_entry(Path(f'entries/N{n}.txt').read_text()); exact=analyze(F,P0)
    if not exact['valid']:raise ValueError('invalid starting set')
    a=sp.Symbol('a');p=sum(int(v)*a**i for i,v in enumerate(F.p));roots=sp.nroots(p,maxsteps=250);real=[complex(x).real for x in roots if abs(complex(x).imag)<1e-9];root=max(real)
    tab=json.loads(Path('records_int_now.json').read_text());base=max(exact['score'],next(x['score'] for x in tab['records'] if x['n']==n))
    for row in json.loads(Path('final_rows.json').read_text()):
      if row['N']==n:base=max(base,row['score'])
    checked=0;generated=0;best=exact['score'];start=time.time();saved=[]
    degrees=[sum(i in line for line in exact['sizes'][4]) for i in range(n)]
    for runno in range(rounds):
      if runno%2:
       rem=set(rng.sample(range(n),depth))
      else:
       pool=sorted(range(n),key=lambda i:degrees[i]+rng.random()*3)[:max(depth*4,min(n,16))];rem=set(rng.sample(pool,depth))
      layer=[tuple(p for i,p in enumerate(P0) if i not in rem)];level_results=[]
      for level in range(depth):
       nxt=[];seen=set()
       for points in layer:
        P=to_float_pts(F,points,root);ls=floatscore(P,tol)
        if ls[5]:continue
        cand=proposals(P,ls,rng,pairs32,pairs22,tol);generated+=len(cand)
        cand.sort(key=lambda v:v[0]+.1*min(v[1],10)+rng.random()*.1,reverse=True)
        for gain,trips,(i,j),(k,l),approx in cand[:max(width*2,8)]:
         L1=line_key(F,points[i],points[j]);L2=line_key(F,points[k],points[l]);q=exact_point_from_lines(F,L1,L2)
         if q is None or q in points:continue
         qp=to_float_pts(F,[q],root)[0]
         if not np.all(np.isfinite(qp)) or np.linalg.norm(qp-approx)>tol*max(1.,np.max(abs(P)))*10:continue
         child=points+(q,);sig=tuple(sorted(child))
         if sig in seen:continue
         seen.add(sig)
         if len(child)==n:
          checked+=1;res=analyze(F,child)
          if not res['valid']:continue
          score=res['score'];best=max(best,score)
          if score>base:
           Path('research/candidates').mkdir(exist_ok=True,parents=True)
           name=f'algebraic_n{n}_s{score}_{seednum}_{runno}_{checked}'
           Path(f'research/candidates/{name}.txt').write_text(entry_text(F,child,[int(v) for v in F.p])+'\n')
           saved.append(name);base=score;print('NEW',name,'score',score,flush=True)
          quality=score
         else:
          quality=len(ls[4])+gain+.07*(len(ls[3])+trips)+rng.random()*.02
         nxt.append((quality,child))
       nxt.sort(key=lambda t:t[0],reverse=True);layer=[x[1] for x in nxt[:width]];level_results.append([len(nxt),len(layer)])
       if not layer:break
      print(json.dumps({'n':n,'round':runno+1,'depth':depth,'source':exact['score'],'baseline':base,'best_exact':best,'generated':generated,'exact_final_checks':checked,'beam':level_results,'seconds':round(time.time()-start,2)}),flush=True)
    out={'n':n,'depth':depth,'width':width,'rounds':rounds,'source':exact['score'],'baseline_at_end':base,'best_exact':best,'generated':generated,'exact_final_checks':checked,'saved':saved,'seed':seednum,'pairs32':pairs32,'pairs22':pairs22,'tol':tol,'seconds':round(time.time()-start,2)}
    Path(result).write_text(json.dumps(out,indent=2)+'\n');print('RESULT',json.dumps(out),flush=True)
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--n',type=int,required=True);ap.add_argument('--depth',type=int,default=3);ap.add_argument('--width',type=int,default=3);ap.add_argument('--rounds',type=int,default=5);ap.add_argument('--seed',type=int,default=123);ap.add_argument('--pairs32',type=int,default=300);ap.add_argument('--pairs22',type=int,default=200);ap.add_argument('--tol',type=float,default=1e-8);ap.add_argument('--result',required=True);a=ap.parse_args();run(a.n,a.depth,a.width,a.rounds,a.seed,a.pairs32,a.pairs22,a.tol,a.result)
