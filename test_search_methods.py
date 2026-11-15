#!/usr/bin/env python3
"""Deterministic smoke tests of proposal, replay and projective rechart paths."""
import random
from pathlib import Path
import numpy as np,sympy as sp
from algfield import parse_entry
from orchard import analyze
from extend import to_float_pts
from beam_algebraic import floatscore
from beam_numeric_replay import candidates,replay
from anneal_algebraic import replay as replay_moves

def main():
 F,P=parse_entry(Path('entries/N40.txt').read_text());x=sp.Symbol('a');poly=sum(int(c)*x**i for i,c in enumerate(F.p))
 root=max(complex(z).real for z in sp.nroots(poly,maxsteps=250) if abs(complex(z).imag)<1e-9);Q=to_float_pts(F,P,root);L=floatscore(Q,1e-8)
 assert len(L[4])==70 and not L[5]
 V=candidates(Q,L,random.Random(8),500,500,1e-8); assert V
 gain,trips,ij,kl,q=max(V,key=lambda v:v[0]);result=replay(F,P,((ij,kl),));assert result is not None
 E,res=result;assert res['score']==70+gain and len(E)==41
 rem=0;part=np.delete(Q,rem,axis=0);ls=floatscore(part,1e-8);V2=candidates(part,ls,random.Random(9),500,500,1e-8);assert V2
 g,t,ij,kl,q=max(V2,key=lambda v:v[0]);moved=replay_moves(F,P,[(rem,ij,kl)]);assert moved is not None
 E,r=moved;assert len(E)==40 and r['score']==len(ls[4])+g
 R,T=parse_entry(Path('research/recharts/n40_u1-5_v1-7.txt').read_text());rr=analyze(R,T);assert rr['valid'] and rr['score']==70
 print('PASS: numeric proposal -> exact growth replay, exact replacement replay, and projective rechart')
if __name__=='__main__':main()
