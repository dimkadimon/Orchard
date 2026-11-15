#!/usr/bin/env python3
"""Independent determinant census of the recorded projective five-point seeds."""
from itertools import combinations
from collections import Counter
from math import comb
import json

def det(a,b,c):
 return a[0]*(b[1]*c[2]-b[2]*c[1])-a[1]*(b[0]*c[2]-b[2]*c[0])+a[2]*(b[0]*c[1]-b[1]*c[0])

def census(P):
 # Each pair's complete point membership is calculated without line normalization.
 # Distinct lines are deduplicated by their point sets.
 n=len(P);assert len(set(map(tuple,P)))==n
 lines={frozenset(k for k in range(n) if det(P[i],P[j],P[k])==0) for i,j in combinations(range(n),2)}
 h=dict(sorted(Counter(map(len,lines)).items()))
 assert sum(comb(k,2)*v for k,v in h.items())==comb(n,2)
 return h

def verify_file(path):
    D=json.load(open(path))
    for n,row in D['results'].items():
        h=census(row['homogeneous_points'])
        assert h=={int(k):v for k,v in row['histogram'].items()}
        assert max(h)<=5 and h.get(5,0)==row['score']
        print(path, 'n='+n, 'histogram='+str(h), 'no-six=True')

if __name__ == '__main__':
    verify_file('research/five_point_seed_deletions.json')
    verify_file('research/five_point_23_seed_deletions.json')
