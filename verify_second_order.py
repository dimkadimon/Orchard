"""Brute-force verification: enumerate ALL degree sequences (not just balanced)
to confirm the majorisation claim, for n=11..18."""
from math import comb
from itertools import product
def w(j):
    if j<=2: return 0
    if j==3: return 3
    return comb(j,2)+(j-3)

def gen(n,R,M):
    # degree sequences as counts s_0..s_R with sum n, sum j*s_j = M
    def rec(j,rem_n,rem_m,acc):
        if j==R:
            if rem_m==R*rem_n:
                yield acc+[rem_n]
            return
        maxc=min(rem_n, rem_m//j if j>0 else rem_n)
        for c in range(maxc+1):
            yield from rec(j+1,rem_n-c,rem_m-j*c,acc+[c])
    yield from rec(0,n,M,[])

def feas_seq(s,L):
    D=comb(L,2)-sum(c*comb(j,2) for j,c in enumerate(s))
    if D<0: return False
    if sum(c*w(j) for j,c in enumerate(s))>comb(L,2)-3: return False
    return True

print("n  L  balanced-feasible  any-sequence-feasible  agree")
allok=True
for n in range(11,19):
    R=(n-1)//3
    M=(comb(n,2)-3)//7
    for L in range(max(3,M-3),M+1):
        if 1+3*L<=n: continue   # pencil escape, skip
        seqs=list(gen(n,R,4*L))
        anyf=any(feas_seq(s,L) for s in seqs)
        q,rem=divmod(4*L,n)
        bal=[0]*(R+1)
        if q>R or (rem and q+1>R):
            balf=False
        else:
            bal[q]+=n-rem
            if rem: bal[q+1]+=rem
            balf=feas_seq(bal,L)
        ok=(anyf==balf); allok&=ok
        print("%2d %3d   %-5s   %-5s  (%d seqs)  %s"%(n,L,balf,anyf,len(seqs),"OK" if ok else "MISMATCH"))
print("ALL AGREE" if allok else "DISAGREEMENT FOUND")
