"""Second-order Melchior bound for the no-five 4-orchard problem.

Let S be n points, no 5 collinear, L = t_4 rich (4-point) lines.
Let A = arrangement of the L rich lines (projective).
 s_j = #{p in S : exactly j rich lines through p}
 e_j = #{vertices of A NOT in S, of multiplicity j}
 D   = # pairs of rich lines meeting outside S = sum_j C(j,2) e_j >= e_2
Melchior on A (L>=3, not all concurrent):
 s_2 + e_2 >= 3 + sum_{j>=4} (j-3)(s_j+e_j)
 => s_2 + D >= 3 + sum_{j>=4}(j-3) s_j
Substituting D = C(L,2) - sum_j C(j,2) s_j:
 sum_{j>=3} w_j s_j <= C(L,2) - 3,  w_j = C(j,2)+(j-3) for j>=4, w_3=3, w_<=2=0.
Also D>=0: sum_j C(j,2) s_j <= C(L,2).
Both left sides are CONVEX in j (w: 0,0,0,3,7,12,18,25,... diffs 0,0,3,4,5,6,7).
With sum s_j = n and sum j s_j = 4L fixed, the BALANCED degree sequence
majorisation-minimises BOTH simultaneously => exact feasibility test.
"""
from math import comb
import json

def w(j):
    if j<=2: return 0
    if j==3: return 3
    return comb(j,2)+(j-3)

def feasible(n,L):
    if L<3: return True,"trivial"
    R=(n-1)//3
    if L> n*R//4: return False,"degree cap"
    if n < 1+3*L:
        pass
    else:
        # all-concurrent pencil is possible; Melchior on A not applicable
        concurrent_ok=True
    q,rem=divmod(4*L,n)
    if q+ (1 if rem else 0) > R: return False,"degree cap (balanced)"
    s={q:n-rem}
    if rem: s[q+1]=rem
    D=comb(L,2)-sum(c*comb(j,2) for j,c in s.items())
    if D<0: return False,"pair count"
    lhs=sum(c*w(j) for j,c in s.items())
    if lhs>comb(L,2)-3:
        if n>=1+3*L: return True,"pencil escape"
        return False,"2nd-order Melchior"
    return True,"ok"

def melchior(n): return (comb(n,2)-3)//7

if __name__=='__main__':
    known={11:6,12:7,13:9,14:10,15:12,16:15,17:16,18:18,19:20,20:23,21:28}
    out={}
    print(" n  Melchior  2nd-order  gain  reason")
    for n in range(11,201):
        M=melchior(n)
        best=M
        for L in range(M,0,-1):
            ok,why=feasible(n,L)
            if ok:
                best=L; break
        gain=M-best
        out[n]={'melchior':M,'second_order':best,'gain':gain}
        if n<=60 or gain>0:
            print("%3d %8d %10d %5d  %s"%(n,M,best,gain,why if gain else ''))
    json.dump(out,open('second_order_results.json','w'),indent=1)
