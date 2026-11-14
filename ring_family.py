"""Numerical exploration of a four-ring cosine-ratio candidate family.

Radii are rounded and restricted to a finite numerical range; collinearities are
counted using floating-point tolerance. Neither enumeration completeness nor
family optimality is proved by this script. In particular, the old parity-only
chord criterion and the assertion that every pair shares a nondiametral chord
are false in general. See paper/orchard.tex and docs/ring_structure.md.
"""
import numpy as np, json, math, itertools, sys, time

def census_fast(P, tol=1e-9):
    """Robust line census by angular grouping at each point.

    At point i, group the directions to all other points modulo pi. A group of size g
    means a line with g+1 points. Count it only when i is its smallest index, so every
    line is counted once. Tolerance-based grouping avoids the rounding-boundary splits
    that a hash-of-rounded-coordinates approach suffers from.
    """
    n = len(P)
    n4 = n5 = 0
    for i in range(n):
        d = P[i+1:] - P[i]
        if len(d) == 0: break
        ang = np.arctan2(d[:, 1], d[:, 0]) % np.pi
        order = np.argsort(ang)
        a = ang[order]; idx = np.arange(i+1, n)[order]
        # split into groups of near-equal angle
        brk = np.where(np.diff(a) > tol)[0] + 1
        groups = np.split(np.arange(len(a)), brk)
        # merge wrap-around group (angle near pi with angle near 0)
        if len(groups) > 1 and (np.pi - a[-1]) + a[0] <= tol:
            groups[0] = np.concatenate([groups[-1], groups[0]])
            groups = groups[:-1]
        for g in groups:
            r = len(g) + 1
            if r < 4: continue
            if idx[g].min() < i: continue   # not canonical
            if r == 4: n4 += 1
            else: n5 += 1
    return n4, n5

def points_of(rs, clss, m):
    P = []
    for r, c in zip(rs, clss):
        a = 2*np.pi*np.arange(m)/m + c*np.pi/m
        P.append(np.stack([r*np.cos(a), r*np.sin(a)], axis=1))
    return np.concatenate(P, axis=0)

def enumerate_family(m, k, target, verbose=True, keep=8):
    cs = [abs(math.cos(math.pi*a/m)) for a in range(1, m//2+1)]
    cs = [c for c in cs if c > 1e-12]
    ratios = sorted({round(a/b, 10) for a in cs for b in cs})
    ratios = [r for r in ratios if 0.02 < r < 50]
    cand = sorted(set([1.0] + ratios))
    best = -1; bestcfg = None; tested = 0; hits = []
    for combo in itertools.combinations_with_replacement(cand, k-1):
        rs = (1.0,) + combo
        if len({round(x, 9) for x in rs}) < k: continue
        for cl in itertools.product((0, 1), repeat=k-1):
            clss = (0,) + cl
            P = points_of(rs, clss, m)
            n4, n5 = census_fast(P)
            tested += 1
            if n5: continue
            if n4 > best:
                best = n4; bestcfg = (rs, clss)
            if n4 > target:
                hits.append({'n': k*m, 'm': m, 'n4': n4, 'radii': list(rs), 'classes': list(clss)})
                if verbose: print(f"    IMPROVE n={k*m} n4={n4} > {target} radii={[round(x,8) for x in rs]} cls={clss}", flush=True)
    return {'m': m, 'k': k, 'n': k*m, 'tested': tested, 'max_n4': best,
            'radii': list(bestcfg[0]) if bestcfg else None,
            'classes': list(bestcfg[1]) if bestcfg else None,
            'target': target, 'hits': hits}

if __name__ == '__main__':
    from ring_collision import full_count
    # cross-check the fast census against the slow one on the four known records
    known = [(7, [1.0,4.04891734,0.356895868,0.692021472], [0,0,1,1], 42),
             (9, [1.0,0.8152074690959046,0.5320888862379561,0.18479253090409542], [0,1,0,1], 63),
             (10,[1.0,0.85065080835204,0.6180339887498949,0.32491969623290634], [0,1,0,1], 70),
             (12,[1.0,0.8965754721680536,0.7320508075688773,0.5176380902050416], [0,1,0,1], 96)]
    for m, rs, cl, exp in known:
        P = points_of(rs, cl, m)
        f = census_fast(P); s = full_count(list(zip(rs, [c*np.pi/m for c in cl])), m)
        print(f"check m={m}: fast={f} slow={s[:2]} expected={exp}", flush=True)
        assert f[0] == exp == s[0], (m, f, s, exp)
    print("census cross-check OK\n", flush=True)

    records = {r['N']: r['score'] for r in json.load(open('final_rows.json'))}
    ints = {r['n']: r['score'] for r in json.load(open('records_int_now.json'))['records']}
    cases = [(7,4),(9,4),(10,4),(12,4),(8,4),(11,4),(6,4),(5,4),(13,4),(14,4)]
    if len(sys.argv) > 1 and sys.argv[1] == 'all':
        cases += [(m,k) for k in (3,5,6) for m in range(5,21) if 16 <= k*m <= 60]
    res = []
    for m, k in cases:
        n = k*m
        target = max(records.get(n,0), ints.get(n,0))
        t = time.time()
        r = enumerate_family(m, k, target)
        res.append(r)
        flag = "MATCHES TARGET (NUMERICAL)" if r['max_n4'] == target else ("BEATS" if r['max_n4'] > target else "below")
        print(f"m={m:2d} k={k} n={n:2d}: tested {r['tested']:6d}  max n4={r['max_n4']:3d}  target={target:3d}  {flag}  ({time.time()-t:.1f}s)", flush=True)
    json.dump(res, open('ring_family_complete.json','w'), indent=1)
    print("\nTOTAL IMPROVEMENTS:", sum(len(r['hits']) for r in res), flush=True)
