#!/usr/bin/env python3
"""Exact projective delete-and-regrow beam search for no-five four-point orchards.

States are sets of primitive integer homogeneous points, including points at infinity.
Each generation adds an intersection of lines spanned by the current points. Candidate
2/2 intersections permit temporarily score-neutral steps; 3/2 and 3/3 intersections
can produce new quadruple lines. Branching is bounded and randomized, so this is NOT
an exhaustive search. All candidate incidence tests and line counts use integer
arithmetic. Record candidates are output in a rational affine chart and can be
independently checked by verify_rational_candidate.py.
"""
import argparse, collections, itertools, json, math, os, random, time
from fractions import Fraction
from pathlib import Path


def normalize(v):
    g=math.gcd(math.gcd(abs(v[0]),abs(v[1])),abs(v[2]))
    if not g: return None
    v=tuple(int(x//g) for x in v)
    return tuple(-x for x in v) if next(x for x in v if x)<0 else v

def cross(p,q):
    a,b,c=p; d,e,f=q
    return normalize((b*f-c*e,c*d-a*f,a*e-b*d))

def incident(l,p):return sum(a*b for a,b in zip(l,p))==0

def make_state(points):
    P=tuple(sorted(points)); assert len(P)==len(set(P))
    lines={}
    for i,j in itertools.combinations(range(len(P)),2):
        l=cross(P[i],P[j]); m=lines.get(l,0);lines[l]=m|(1<<i)|(1<<j)
    if any(v.bit_count()>4 for v in lines.values()):return None
    typ={j:[] for j in (2,3,4)}
    for l,v in lines.items():typ[v.bit_count()].append((l,v))
    return P,typ,len(typ[4])

def candidate_points(state, rng, pair32=800, pair22=300, maxbits=96):
    P,typ,_=state; two=typ[2];three=typ[3];four=typ[4]; seen=set(P); candidates={}
    # 3/3 intersections are always considered, 3/2 and 2/2 are bounded samples.
    ps=[]
    for A,B in itertools.combinations(three,2):
        if not A[1]&B[1]:ps.append((A,B))
    if three and two:
        for _ in range(pair32):
            A=rng.choice(three);B=rng.choice(two)
            if not A[1]&B[1]:ps.append((A,B))
    if len(two)>1:
        for _ in range(pair22):
            A,B=rng.sample(two,2)
            if not A[1]&B[1]:ps.append((A,B))
    for A,B in ps:
        p=cross(A[0],B[0])
        if p is None or p in seen or max(abs(x).bit_length() for x in p)>maxbits:continue
        seen.add(p)
        if any(incident(l,p) for l,_ in four):continue
        gain=sum(incident(l,p) for l,_ in three)
        # New triple lines can unlock future gains. Actual gain never assumes them.
        triples=sum(incident(l,p) for l,_ in two)
        candidates[p]=(gain,triples)
    return candidates

def chart(P):
    # Invertible projective map (X,Y,Z)->(X,Y,aX+bY+Z). Finite set permits a,b.
    for a in range(-20,21):
      for b in range(-20,21):
        den=[a*x+b*y+z for x,y,z in P]
        if all(den):
            Q=[(Fraction(x,d),Fraction(y,d)) for (x,y,z),d in zip(P,den)]
            scale=math.lcm(*(x.denominator for p in Q for x in p))
            return [(int(x*scale),int(y*scale)) for x,y in Q], (a,b,1)
    raise ValueError('failed to choose affine chart')

def exact_census(P):
    s=make_state(P)
    if s is None:raise ValueError('invalid no-five set')
    return {j:len(s[1][j]) for j in (2,3,4)}

def save_candidate(P,score,n,metadata,outdir):
    q,chartvec=chart(P);hist=exact_census(P)
    assert hist[4]==score and len(q)==n
    key=f'n{n}_s{score}_{time.time_ns()}'
    os.makedirs(outdir,exist_ok=True)
    txt='a = 0\n'+', '.join(f'({x}, {y})' for x,y in q)+'\n'
    Path(outdir,f'{key}.txt').write_text(txt)
    Path(outdir,f'{key}.json').write_text(json.dumps({'n':n,'score':score,'histogram':hist,'chart':chartvec,'homogeneous':[list(x) for x in P],'search':metadata},indent=2)+'\n')
    return key

def load(n,source):
    if source=='integer':
        data=json.loads(Path('integer_solutions.json').read_text())[str(n)]
        P=[normalize((int(x),int(y),1)) for x,y in data['points']]
    else:
        from algfield import parse_entry
        F, pts=parse_entry(Path(f'entries/N{n}.txt').read_text())
        if F.d!=1:raise ValueError('entry is not rational; use integer source')
        P=[]
        for x,y in pts:
            a,b=x[0],y[0];d=math.lcm(a.denominator,b.denominator)
            P.append(normalize((int(a*d),int(b*d),d)))
    st=make_state(P)
    if st is None:raise ValueError('source has five collinear')
    return st

def beam(seed,n,depth,width,rounds,rng,limits,baseline,outdir):
    best=seed[2]; evaluated=0; generated=0; max_final=None; start=time.time(); saved=[]
    for run in range(rounds):
        # Vary the deletion difficulty: not only the smallest-loss deletion sets.
        degree=[0]*n
        for l,mask in seed[1][4]:
            for i in range(n):degree[i]+=bool(mask>>i&1)
        if run%2==0:
            pool=sorted(range(n),key=lambda i:degree[i]+rng.random()*3.0)[:max(depth*4,min(n,16))]
            removed=set(rng.sample(pool,depth))
        else: removed=set(rng.sample(range(n),depth))
        root=make_state([p for i,p in enumerate(seed[0]) if i not in removed]); assert root
        layer=[root]; depth_report=[]
        for level in range(depth):
            proposals={}
            for st in layer:
                candidates=candidate_points(st,rng,*limits)
                generated+=len(candidates)
                ranked=sorted(candidates.items(),key=lambda kv:(kv[1][0]+0.11*min(kv[1][1],10)+rng.random()*0.12),reverse=True)
                # Favor some unusual low-immediate-gain alternatives, not just greedy adds.
                ranked=ranked[:max(width*5,35)]
                for p,(gain,triples) in ranked:
                    child=make_state(st[0]+(p,));evaluated+=1
                    if child is None:continue
                    assert child[2]==st[2]+gain
                    if child[0] in proposals:continue
                    # Score plus potential; preserve varying topologies with exact state deduplication.
                    potential=len(child[1][3]); quality=child[2]+0.06*potential+rng.random()*0.08
                    proposals[child[0]]=(quality,child)
                    if len(child[0])==n:
                        max_final=child[2] if max_final is None else max(max_final,child[2])
                        if child[2]>best:best=child[2]
                        if child[2]>baseline:
                            meta={'round':run,'depth':depth,'width':width,'source_score':seed[2],'baseline':baseline}
                            key=save_candidate(child[0],child[2],n,meta,outdir);saved.append(key);baseline=child[2]
                            print('NEW',key,'score',child[2],'elapsed',round(time.time()-start,1),flush=True)
            layer=[x[1] for x in sorted(proposals.values(),key=lambda x:x[0],reverse=True)[:width]]
            depth_report.append([len(proposals),len(layer)])
            if not layer:break
        print(json.dumps({'n':n,'round':run+1,'depth':depth,'source':seed[2],'baseline':baseline,'best_completed':max_final,'generated':generated,'states_evaluated':evaluated,'beam':depth_report,'seconds':round(time.time()-start,2)}),flush=True)
    return {'n':n,'depth':depth,'width':width,'rounds':rounds,'source':seed[2],'baseline_at_end':baseline,'best_completed':max_final,'generated':generated,'states_evaluated':evaluated,'saved':saved,'seconds':round(time.time()-start,2)}

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--n',type=int,required=True);ap.add_argument('--depth',type=int,default=3);ap.add_argument('--width',type=int,default=8);ap.add_argument('--rounds',type=int,default=12);ap.add_argument('--seed',type=int,default=12345);ap.add_argument('--source',choices=['integer','entry'],default='integer');ap.add_argument('--pairs32',type=int,default=800);ap.add_argument('--pairs22',type=int,default=300);ap.add_argument('--maxbits',type=int,default=96);ap.add_argument('--outdir',default='research/candidates');ap.add_argument('--result',default='');a=ap.parse_args()
    tab=json.loads(Path('records_int_now.json').read_text()); baselines={x['n']:x['score'] for x in tab['records']}
    for row in json.loads(Path('final_rows.json').read_text()):baselines[row['N']]=max(baselines.get(row['N'],0),row['score'])
    s=load(a.n,a.source);baseline=max(s[2],baselines.get(a.n,0))
    result=beam(s,a.n,a.depth,a.width,a.rounds,random.Random(a.seed),(a.pairs32,a.pairs22,a.maxbits),baseline,a.outdir)
    result['seed']=a.seed;result['source_kind']=a.source;result['pairs32']=a.pairs32;result['pairs22']=a.pairs22;result['maxbits']=a.maxbits
    if a.result:Path(a.result).write_text(json.dumps(result,indent=2)+'\n')
    print('RESULT',json.dumps(result),flush=True)
if __name__=='__main__':main()
