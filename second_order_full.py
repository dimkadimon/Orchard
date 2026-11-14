"""Second-order Melchior bound, valid for BOTH t_4^(5)(n) and T_4(n).
Also derive exact values by combining with known constructions."""
from math import comb
import json
def w(j):
    if j<=2: return 0
    if j==3: return 3
    return comb(j,2)+(j-3)
def feasible(n,L):
    if L<3: return True
    if 1+3*L<=n: return True          # pencil: Melchior on A not applicable
    R=(n-1)//3
    q,rem=divmod(4*L,n)
    if q+(1 if rem else 0)>R: return False
    s=[(q,n-rem)]+([(q+1,rem)] if rem else [])
    if sum(c*comb(j,2) for j,c in s)>comb(L,2): return False
    if sum(c*w(j) for j,c in s)>comb(L,2)-3: return False
    return True
def bound(n):
    M=(comb(n,2)-3)//7
    for L in range(M,0,-1):
        if feasible(n,L): return L
    return 0
lower_no5={11:6,12:7,13:9,14:10,15:12,16:15,17:16,18:18,19:20,20:23}
if __name__=='__main__':
    rows=[]
    print(" n | Melchior | 2nd-order | known t4^(5) | T4 lower | T4 exact?")
    for n in range(11,61):
        M=(comb(n,2)-3)//7; B=bound(n)
        lo=lower_no5.get(n)
        ex = (lo is not None and B==lo)
        rows.append({'n':n,'melchior':M,'second_order':B,'known_no5':lo,'exact':ex})
        if n<=22:
            print("%3d | %8d | %9d | %12s | %8s | %s"%(n,M,B,lo,lo,("YES = %d"%B) if ex else ""))
    json.dump(rows,open('second_order_full.json','w'),indent=1)
    print("\nExact values obtained (2nd-order bound == known construction):")
    print([r['n'] for r in rows if r['exact']])
    print("\nImprovements over Melchior:", [(r['n'],r['melchior'],r['second_order']) for r in rows if r['second_order']<r['melchior']])
