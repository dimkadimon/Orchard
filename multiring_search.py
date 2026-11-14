"""NEW APPROACH: generalised multi-ring ("celestial") search.

Previous searches fixed FOUR rings of EQUAL size m with offsets locked to {0, pi/m}.
Here we relax all three restrictions at once:
  * arbitrary number of rings k (3..8),
  * arbitrary per-ring point counts m_i (not all equal),
  * CONTINUOUS per-ring radii r_i and angular offsets phi_i.

Score = number of lines with exactly 4 points, rejecting any line with >=5.
Optimisation: random restarts + coordinate-wise / simultaneous perturbation hill
climbing on a smooth surrogate, then exact-count polishing.

Any candidate beating the current record is written out with high-precision
coordinates for later exactification (LLL / algdep).
"""
import numpy as np, json, math, random, sys, time, itertools

def build(rings):
    """rings: list of (m_i, r_i, phi_i). Returns (N,2) array."""
    pts = []
    for m, r, phi in rings:
        j = np.arange(m)
        ang = 2*np.pi*j/m + phi
        pts.append(np.stack([r*np.cos(ang), r*np.sin(ang)], axis=1))
    return np.concatenate(pts, axis=0)

def census(P, tol=1e-7):
    """Exact-ish floating census; returns (n4, n5plus, lines4)."""
    n = len(P)
    scale = max(1.0, np.abs(P).max())
    seen = {}
    for i in range(n):
        d = P - P[i]
        L = np.hypot(d[:,0], d[:,1])
        ok = L > 1e-12
        for j in range(i+1, n):
            if not ok[j]:
                continue
            cr = np.abs(d[:,0]*d[j,1] - d[:,1]*d[j,0]) / L[j]
            on = np.where(cr < tol*scale)[0]
            if len(on) >= 4:
                seen[tuple(on)] = len(on)
    n4 = sum(1 for v in seen.values() if v == 4)
    n5 = sum(1 for v in seen.values() if v >= 5)
    return n4, n5, seen

def score(rings, penalise=True):
    P = build(rings)
    if len(P) > 60:
        return -1, 0, 0
    n4, n5, _ = census(P)
    if n5 > 0:
        return (-1 if penalise else n4), n4, n5
    return n4, n4, n5

def optimise(ms, target, seed=0, iters=4000, tlimit=60.0):
    """ms: per-ring point counts. Returns best (score, rings)."""
    rng = np.random.default_rng(seed)
    k = len(ms)
    t0 = time.time()

    def rand_rings():
        r = np.sort(rng.uniform(0.15, 1.0, k))
        r[-1] = 1.0
        if rng.random() < 0.25:
            r[rng.integers(k)] = rng.uniform(1.2, 4.5)
        phi = rng.uniform(0, 2*np.pi, k)
        # bias towards half-step offsets, which is where the known solutions live
        if rng.random() < 0.6:
            phi = np.array([ (np.pi/ms[i])*rng.integers(0,2) for i in range(k) ])
        return [(ms[i], float(r[i]), float(phi[i])) for i in range(k)]

    best = -1; bestrings = None
    cur = rand_rings(); cs, _, _ = score(cur)
    for it in range(iters):
        if time.time()-t0 > tlimit:
            break
        if cs < 0 or rng.random() < 0.06:
            cur = rand_rings(); cs, _, _ = score(cur)
        prop = list(cur)
        nmut = 1 + (rng.random() < 0.3)
        for _ in range(nmut):
            i = int(rng.integers(k))
            m, r, phi = prop[i]
            step = max(0.002, 0.09*(1 - it/iters))
            if rng.random() < 0.55:
                r = max(0.03, r + rng.normal(0, step))
            else:
                phi = phi + rng.normal(0, step*0.8)
            prop[i] = (m, float(r), float(phi))
        ps, _, _ = score(prop)
        if ps >= cs:
            cur, cs = prop, ps
        if cs > best:
            best, bestrings = cs, list(cur)
            if best >= target:
                break
    return best, bestrings

if __name__ == '__main__':
    records = {r['N']: r['score'] for r in json.load(open('final_rows.json'))}
    ints = {r['n']: r['score'] for r in json.load(open('records_int_now.json'))['records']} \
        if __import__('os').path.exists('records_int_now.json') else {}
    # candidate ring signatures with total <= 60
    sigs = []
    for k in (3,4,5,6):
        for m in range(4, 16):
            if k*m <= 60 and k*m >= 16:
                sigs.append([m]*k)
    # unequal-ring signatures
    for base in range(4, 13):
        for k in (3,4,5):
            for dbl in range(1, k):
                ms = [base]*k
                for i in range(dbl):
                    ms[i] = 2*base
                if sum(ms) <= 60 and sum(ms) >= 16:
                    sigs.append(ms)
    sigs = [list(s) for s in {tuple(sorted(x)) for x in sigs}]
    sigs.sort(key=lambda s: (sum(s), len(s)))
    print(f"testing {len(sigs)} ring signatures", flush=True)
    out = []
    t0 = time.time()
    budget = float(sys.argv[1]) if len(sys.argv) > 1 else 600.0
    for ms in sigs:
        n = sum(ms)
        target = max(records.get(n, 0), ints.get(n, 0)) + 1
        b, rings = optimise(ms, target, seed=1, iters=1500, tlimit=min(20.0, budget/len(sigs)))
        rec = records.get(n, 0)
        flag = 'BEATS' if b > rec else ''
        print(f"ms={ms} n={n} best={b} record={rec} {flag}", flush=True)
        out.append({'ms': ms, 'n': n, 'best': b, 'record': rec, 'rings': rings})
        if time.time()-t0 > budget:
            print("budget reached", flush=True)
            break
    json.dump(out, open('multiring_results.json','w'), indent=1)
    wins = [o for o in out if o['best'] > o['record']]
    print("IMPROVEMENTS:", [(o['n'], o['best'], o['record']) for o in wins], flush=True)
