"""General-k second-order test. Melchior's threshold is always 4, so w is the SAME
function for every k; only sum r_p >= kL and the cap R=floor((n-1)/(k-1)) change."""
from math import comb
def w(j):
    if j<=2: return 0
    if j==3: return 3
    return comb(j,2)+(j-3)
def melchior_k(n,k):
    # C(n,2) >= 3 + sum_{j>=k}[C(j,2)+(j-3)] t_j ; min bracket at j=k
    return (comb(n,2)-3)//(comb(k,2)+k-3)
def feasible(n,L,k):
    if L<3: return True
    if 1+(k-1)*L<=n: return True      # pencil escape
    R=(n-1)//(k-1)
    q,rem=divmod(k*L,n)
    if q+(1 if rem else 0)>R: return False
    s=[(q,n-rem)]+([(q+1,rem)] if rem else [])
    if sum(c*comb(j,2) for j,c in s)>comb(L,2): return False
    if sum(c*w(j) for j,c in s)>comb(L,2)-3: return False
    return True
def bound(n,k):
    M=melchior_k(n,k)
    for L in range(M,0,-1):
        if feasible(n,L,k): return L
    return 0
if __name__=='__main__':
    for k in (4,5,6,7):
        gains=[]
        for n in range(k+1,121):
            M=melchior_k(n,k); B=bound(n,k)
            if B<M: gains.append((n,M,B))
        print(f"k={k}: improvements over Melchior for n<=%d: %s"%(120,gains))
