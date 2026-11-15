#!/usr/bin/env python3
"""Integer checks for all upper-bound comparisons displayed in the paper."""
from math import comb
from second_order_k import bound, melchior_k, w

old=[3,4,6,7,9,10,12,14]
new=[2,3,5,6,8,9,11,13]
for n,c,b in zip(range(8,16),old,new):
    assert (melchior_k(n,4),bound(n,4))==(c,b)

for n,k,limit,classical in [(24,5,19,22),(25,5,23,24),(32,6,19,27)]:
    assert melchior_k(n,k)==classical and bound(n,k)==limit
    print(f'n={n}, k={k}: classical={classical}, combined={limit}')
    for L in range(limit+1,classical+1):
        assert n<1+(k-1)*L  # nonconcurrency hypothesis
        q,r=divmod(k*L,n)
        W=(n-r)*w(q)+r*w(q+1)
        C=comb(L,2)-3
        assert W>C
        print(f'  L={L}: balanced weight {W} > {C}')
print('All claimed numerical upper-bound comparisons verified.')
