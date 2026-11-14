"""NEW APPROACH 2: distance-collision search over a real cyclotomic field.

Observation (verified on all four record configurations): almost every 4-point line
meets TWO rings TWICE each -- type (2,2). For a ring of radius r with m equally spaced
points, a chord subtending span a lies on the line
        x . (cos t, sin t) = r*cos(pi*a/m),   t = phi + pi*(2*j+a)/m.
So a (2,2) line exists exactly when
        r_i * cos(pi*a/m) = r_j * cos(pi*b/m)          <-- LINEAR in the radii
with matching normal direction. Hence the useful radii RATIOS are ratios of cosines,
a finite explicit set. This turns an intractable continuous search into a finite one.

We therefore:
  1. build the candidate ratio set R = { cos(pi a/m) / cos(pi b/m) },
  2. precompute, for each ratio and each offset-class difference, how many (2,2)
     lines a pair of rings contributes,
  3. exhaustively maximise the total over k rings,
  4. reject any candidate with a line of >= 5 points, and take a full census.
"""
import numpy as np, json, math, itertools, sys, time
from collections import defaultdict

TOL = 1e-9

def ring_chords(m):
    """All chords of the unit ring with m points, offset 0.
    Returns arrays (c, t) with c = cos(pi*a/m) signed, t = normal angle."""
    ang = 2*np.pi*np.arange(m)/m
    cs, ts = [], []
    for j in range(m):
        for k in range(j+1, m):
            mid = (ang[j]+ang[k])/2
            half = (ang[k]-ang[j])/2
            cs.append(math.cos(half)); ts.append(mid % (2*math.pi))
    return np.array(cs), np.array(ts)

def canon(d, t):
    """Canonical line key: x.(cos t,sin t)=d with d>=0."""
    d = np.array(d, float); t = np.array(t, float)
    flip = d < 0
    d = np.where(flip, -d, d)
    t = np.where(flip, t+np.pi, t) % (2*np.pi)
    # for d==0 the normal is defined mod pi
    t = np.where(d < 1e-14, t % np.pi, t)
    return d, t

def keys_for(r, delta, m, base_c, base_t):
    d, t = canon(r*base_c, base_t + delta)
    return set(zip(np.round(d, 8), np.round(t, 8)))

def pair_table(m, ratios):
    """table[(rho, dcls)] = number of (2,2) lines between ring i (r=1,cls=0)
    and ring j (r=rho, cls=dcls)."""
    bc, bt = ring_chords(m)
    A = keys_for(1.0, 0.0, m, bc, bt)
    tab = {}
    for rho in ratios:
        for dc in (0, 1):
            B = keys_for(rho, dc*np.pi/m, m, bc, bt)
            tab[(rho, dc)] = len(A & B)
    return tab

def full_count(ms_rings, m):
    """Exact-ish full census of the actual point set. Returns (n4, n5, n)."""
    pts = []
    for (r, phi) in ms_rings:
        j = np.arange(m)
        a = 2*np.pi*j/m + phi
        pts.append(np.stack([r*np.cos(a), r*np.sin(a)], axis=1))
    P = np.concatenate(pts, axis=0)
    n = len(P)
    scale = max(1.0, np.abs(P).max())
    seen = {}
    for i in range(n):
        d = P - P[i]
        L = np.hypot(d[:,0], d[:,1])
        for j in range(i+1, n):
            if L[j] < 1e-12: continue
            cr = np.abs(d[:,0]*d[j,1] - d[:,1]*d[j,0])/L[j]
            on = tuple(np.where(cr < 1e-7*scale)[0])
            if len(on) >= 4: seen[on] = len(on)
    n4 = sum(1 for v in seen.values() if v == 4)
    n5 = sum(1 for v in seen.values() if v >= 5)
    return n4, n5, n

def search(m, k, records, ints, topn=40, verbose=True):
    n = k*m
    if n > 60: return []
    half = m//2
    cs = [abs(math.cos(math.pi*a/m)) for a in range(1, half+1)]
    cs = [c for c in cs if c > 1e-12]
    ratios = set()
    for a in cs:
        for b in cs:
            ratios.add(round(a/b, 10))
    ratios = sorted(r for r in ratios if 0.05 < r < 20)
    tab = pair_table(m, ratios)
    # only ratios that actually produce collisions are useful
    useful = [r for r in ratios if tab[(r,0)] > 0 or tab[(r,1)] > 0]
    if verbose:
        print(f"  m={m} k={k} n={n}: {len(ratios)} ratios, {len(useful)} useful", flush=True)
    best = []
    target = max(records.get(n,0), ints.get(n,0))
    # radii: r_0 = 1; others = 1 or a useful ratio (relative to r_0)
    cand = [1.0] + useful
    for combo in itertools.product(cand, repeat=k-1):
        rs = (1.0,) + combo
        if len(set(round(x,9) for x in rs)) < k: continue
        for cls in itertools.product((0,1), repeat=k-1):
            clss = (0,) + cls
            tot = 0
            for i in range(k):
                for j in range(i+1, k):
                    rho = round(rs[j]/rs[i], 10)
                    dc = (clss[j]-clss[i]) % 2
                    tot += tab.get((rho,dc), 0)
            if tot > target:
                best.append((tot, rs, clss))
    best.sort(key=lambda x: -x[0])
    out = []
    for tot, rs, clss in best[:topn]:
        rings = [(rs[i], clss[i]*np.pi/m) for i in range(k)]
        n4, n5, _ = full_count(rings, m)
        if n5 == 0 and n4 > target:
            out.append({'m':m,'k':k,'n':n,'pairbound':tot,'n4':n4,
                        'radii':list(rs),'classes':list(clss),'record':target})
            if verbose:
                print(f"    CANDIDATE n={n} n4={n4} (record {target}) radii={[round(x,8) for x in rs]} classes={clss}", flush=True)
    return out

if __name__ == '__main__':
    records = {r['N']: r['score'] for r in json.load(open('final_rows.json'))}
    ints = {r['n']: r['score'] for r in json.load(open('records_int_now.json'))['records']}
    allout = []
    t0 = time.time()
    budget = float(sys.argv[1]) if len(sys.argv)>1 else 600.0
    cases = [(m,k) for k in (3,4,5,6) for m in range(5,21) if 16 <= k*m <= 60]
    # prioritise 4 rings (the known mechanism), largest n first
    cases.sort(key=lambda mk: (mk[1]!=4, -mk[0]*mk[1]))
    for m,k in cases:
        if time.time()-t0 > budget:
            print("budget reached", flush=True); break
        try:
            allout += search(m,k,records,ints)
        except Exception as e:
            print("  err",m,k,e, flush=True)
    json.dump(allout, open('ring_collision_results.json','w'), indent=1)
    print("IMPROVEMENTS:", [(o['n'],o['n4'],o['record']) for o in allout], flush=True)
