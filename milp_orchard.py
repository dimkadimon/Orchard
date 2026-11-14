#!/usr/bin/env python3
"""Restricted candidate-pool MILP for the no-five four-point orchard problem.

The candidate pool consists of a known integer construction plus exact projective
intersections drawn from several partially deleted seeds. Every line containing four
or more pool points is enumerated exactly. Binary x_i selects points and binary z_l
selects four-point lines; the model enforces |S|=n, <=4 points per line, and z_l=1
only if four pool points on the line are selected. Maximizing sum z_l therefore gives
the true strict score in the *finite chosen pool*. A time-limited solver may not prove
optimality even within that pool. Any new incumbent is verified by a fresh exact census.
"""
import argparse, itertools, json, math, random, time
from pathlib import Path
import numpy as np
from scipy.optimize import milp, Bounds, LinearConstraint
from scipy.sparse import coo_matrix
from beam_orchard import load, make_state, candidate_points, cross, save_candidate


def pool_lines(P):
    lines={}
    for i,j in itertools.combinations(range(len(P)),2):
        l=cross(P[i],P[j]);lines[l]=lines.get(l,0)|(1<<i)|(1<<j)
    return [(l,m) for l,m in lines.items() if m.bit_count()>=4]

def search(n, extra, delete, trials, seednum, seconds, outdir, result):
    rng=random.Random(seednum); start=time.time(); source=load(n,'integer');base=source[0]
    possible={}
    for trial in range(trials):
        rem=set(rng.sample(range(n),delete));st=make_state([p for i,p in enumerate(base) if i not in rem])
        for p,(gain,trips) in candidate_points(st,rng,700,700,64).items():
            if p in base:continue
            quality=gain+0.13*min(trips,10)+0.2*rng.random()
            possible[p]=max(possible.get(p,0),quality)
    cand=sorted(possible,key=lambda p:possible[p],reverse=True)[:extra]
    P=tuple(base)+tuple(cand);m=len(P); lines=pool_lines(P);L=len(lines); nv=m+L
    rows=[];cols=[];vals=[]; lb=[];ub=[]
    # n points exactly
    for i in range(m):rows.append(0);cols.append(i);vals.append(1.)
    lb.append(float(n));ub.append(float(n))
    for t,(_,mask) in enumerate(lines):
        # cap: sum x <= 4; implication: 4z - sum x <= 0
        row1=1+2*t;row2=row1+1
        while mask:
            low=mask&-mask;i=low.bit_length()-1;mask-=low
            rows.extend([row1,row2]);cols.extend([i,i]);vals.extend([1.,-1.])
        rows.append(row2);cols.append(m+t);vals.append(4.)
        lb.extend([-np.inf,-np.inf]);ub.extend([4.,0.])
    mat=coo_matrix((vals,(rows,cols)),shape=(1+2*L,nv)).tocsc()
    c=np.concatenate([np.zeros(m),-np.ones(L)])
    print(json.dumps({'n':n,'candidate_pool':m,'new_points':len(cand),'possible':len(possible),'rich_pool_lines':L,'variables':nv,'time_limit_s':seconds,'seed':seednum}),flush=True)
    sol=milp(c,integrality=np.ones(nv),bounds=Bounds(np.zeros(nv),np.ones(nv)),constraints=LinearConstraint(mat,np.array(lb),np.array(ub)),options={'time_limit':seconds,'mip_rel_gap':0.0,'presolve':True})
    status={'n':n,'source':source[2],'pool':m,'pool_lines':L,'seed':seednum,'extra':extra,'deletions':delete,'trials':trials,'time_limit':seconds,'status':int(sol.status),'message':sol.message,'solver_fun':float(sol.fun) if sol.fun is not None else None,'mip_dual_bound':float(sol.mip_dual_bound) if getattr(sol,'mip_dual_bound',None) is not None else None,'seconds':round(time.time()-start,2),'selected_score':None,'candidate':None}
    if sol.x is not None:
        chosen=[P[i] for i in range(m) if sol.x[i]>0.5]
        if len(chosen)!=n:raise RuntimeError('solver returned nonintegral point selection')
        st=make_state(chosen)
        if st is None:raise RuntimeError('solver returned a forbidden five-point line')
        status['selected_score']=st[2]
        tab=json.loads(Path('records_int_now.json').read_text());b=max(source[2],next(x['score'] for x in tab['records'] if x['n']==n))
        for row in json.loads(Path('final_rows.json').read_text()):
            if row['N']==n:b=max(b,row['score'])
        status['baseline']=b
        if st[2]>b:status['candidate']=save_candidate(st[0],st[2],n,{'method':'MILP','seed':seednum,'pool':m},outdir)
    Path(result).write_text(json.dumps(status,indent=2)+'\n');print('RESULT',json.dumps(status),flush=True)

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--n',type=int,required=True);ap.add_argument('--extra',type=int,default=28);ap.add_argument('--delete',type=int,default=4);ap.add_argument('--trials',type=int,default=5);ap.add_argument('--seed',type=int,default=1234);ap.add_argument('--seconds',type=int,default=45);ap.add_argument('--outdir',default='research/candidates');ap.add_argument('--result',required=True);a=ap.parse_args();search(a.n,a.extra,a.delete,a.trials,a.seed,a.seconds,a.outdir,a.result)
