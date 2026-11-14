"""Optimality of the four records WITHIN the cosine-ratio family, plus mixed-modulus rings.

Part A (complete/exhaustive).
  If every pair of rings contributes at least one (2,2) line then r_i/r_0 lies in the
  finite set R_m = { cos(pi a/m)/cos(pi b/m) }. Normalising r_0 = 1, an exhaustive
  enumeration of R_m^{k-1} x {0,1}^{k-1} therefore covers EVERY such configuration.
  We report the exact maximum, and full-census the top candidates.

Part B (new direction).
  Rings of DIFFERENT moduli m_i, sharing a common offset grid. Chord normals of a
  modulus-m ring live on the grid of step pi/m, so rings with commensurate moduli can
  still collide. This strictly enlarges the search space beyond Part A.
"""
import numpy as np, json, math, itertools, sys, time
from collections import defaultdict
from ring_collision import ring_chords, canon, keys_for, full_count

def partA(m, k, records, ints, topn=12):
    n = k*m
    cs = [abs(math.cos(math.pi*a/m)) for a in range(1, m//2+1)]
    cs = [c for c in cs if c > 1e-12]
    ratios = sorted({round(a/b, 10) for a in cs for b in cs} - {1.0})
    ratios = [r for r in ratios if 0.02 < r < 50]
    bc, bt = ring_chords(m)
    A = keys_for(1.0, 0.0, m, bc, bt)
    tab = {}
    for rho in ratios + [1.0]:
        for dc in (0, 1):
            tab[(rho, dc)] = len(A & keys_for(rho, dc*np.pi/m, m, bc, bt))
    cand = [1.0] + ratios
    scored = []
    for combo in itertools.product(cand, repeat=k-1):
        rs = (1.0,) + combo
        if len({round(x, 9) for x in rs}) < k: continue
        for cls in itertools.product((0, 1), repeat=k-1):
            clss = (0,) + cls
            tot = 0
            for i in range(k):
                for j in range(i+1, k):
                    tot += tab.get((round(rs[j]/rs[i], 10), (clss[j]-clss[i]) % 2), 0)
            scored.append((tot, rs, clss))
    scored.sort(key=lambda x: -x[0])
    best = scored[0][0] if scored else 0
    out = []
    seen = set()
    for tot, rs, clss in scored:
        if len(out) >= topn: break
        sig = tuple(sorted(round(x, 8) for x in rs))
        if sig in seen: continue
        seen.add(sig)
        rings = [(rs[i], clss[i]*np.pi/m) for i in range(k)]
        n4, n5, _ = full_count(rings, m)
        out.append({'m': m, 'k': k, 'n': n, 'pairbound': tot, 'n4': n4, 'n5': n5,
                    'radii': [round(x, 10) for x in rs], 'classes': list(clss)})
    target = max(records.get(n, 0), ints.get(n, 0))
    valid = [o for o in out if o['n5'] == 0]
    bestvalid = max((o['n4'] for o in valid), default=0)
    return {'m': m, 'k': k, 'n': n, 'max_pairbound': best, 'best_valid_n4': bestvalid,
            'target': target, 'top': out}

# ---------------- Part B: mixed moduli ----------------
def chord_keys(m, r, phi):
    bc, bt = ring_chords(m)
    return keys_for(r, phi, m, bc, bt)

def mixed_eval(spec, grid):
    """spec = [(m, r, phi_index)]; returns (pair_total, n4, n5, n)."""
    sets = []
    npts = 0
    for (m, r, pi) in spec:
        sets.append(chord_keys(m, r, pi*np.pi/grid))
        npts += m
    tot = 0
    for i, j in itertools.combinations(range(len(spec)), 2):
        tot += len(sets[i] & sets[j])
    return tot, npts

def mixed_points(spec, grid):
    P = []
    for (m, r, pi) in spec:
        a = 2*np.pi*np.arange(m)/m + pi*np.pi/grid
        P.append(np.stack([r*np.cos(a), r*np.sin(a)], axis=1))
    return np.concatenate(P, axis=0)

def census_pts(P):
    n = len(P); scale = max(1.0, np.abs(P).max()); seen = {}
    for i in range(n):
        d = P - P[i]; L = np.hypot(d[:, 0], d[:, 1])
        for j in range(i+1, n):
            if L[j] < 1e-12: continue
            cr = np.abs(d[:, 0]*d[j, 1] - d[:, 1]*d[j, 0])/L[j]
            on = tuple(np.where(cr < 1e-7*scale)[0])
            if len(on) >= 4: seen[on] = len(on)
    return sum(1 for v in seen.values() if v == 4), sum(1 for v in seen.values() if v >= 5)

def partB(moduli, records, ints, budget=120.0, verbose=True, deadline=None):
    """Greedy + exhaustive-over-radii for a fixed multiset of moduli."""
    n = sum(moduli)
    if n > 60: return []
    grid = math.lcm(*moduli)
    if grid > 60:
        if verbose: print(f"  mixed {moduli}: grid {grid} > 60, skipped", flush=True)
        return []
    target = max(records.get(n, 0), ints.get(n, 0))
    cs = sorted({round(abs(math.cos(math.pi*a/m)), 12) for m in set(moduli)
                 for a in range(1, m//2+1)} - {0.0})
    cs = [c for c in cs if c > 1e-9]
    ratios = sorted({round(a/b, 9) for a in cs for b in cs})
    ratios = [r for r in ratios if 0.03 < r < 30]
    if verbose:
        print(f"  mixed {moduli} n={n} grid={grid}: {len(ratios)} radii ratios, target {target}", flush=True)
    # precompute chord sets per (m, r, phi) lazily
    cache = {}
    def S(m, r, pi):
        key = (m, round(r, 9), pi % (2*grid))
        if key not in cache: cache[key] = chord_keys(m, r, (pi % (2*grid))*np.pi/grid)
        return cache[key]
    out = []
    t0 = time.time()
    k = len(moduli)
    # Ring i is invariant under rotation by 2*pi/m_i, so its offset matters only
    # modulo 2*pi/m_i. Fixing phi_0 = 0 (global rotation), ring i needs at most
    # 2*grid/m_i distinct offsets on the pi/grid lattice -- NOT 2*grid of them.
    phis = [list(range(max(1, 2*grid//m))) for m in moduli]
    rad = [1.0] + ratios
    if k <= 3:
        best = []
        ncomb = (len(rad)**(k-1)) * math.prod(len(p) for p in phis[1:])
        if ncomb > 4_000_000:
            if verbose: print(f"    {ncomb} combos > cap, skipped", flush=True)
            return []
        for rs in itertools.product(rad, repeat=k-1):
            if deadline and time.time() > deadline:
                if verbose: print("    deadline", flush=True)
                break
            rr = (1.0,) + rs
            if len({round(x, 8) for x in rr}) < k: continue
            for ps in itertools.product(*phis[1:]):
                pp = (0,) + ps
                spec = [(moduli[i], rr[i], pp[i]) for i in range(k)]
                tot, _ = mixed_eval(spec, grid)
                if tot > target: best.append((tot, spec))
        best.sort(key=lambda x: -x[0])
        if verbose and best:
            print(f"    best (2,2) total = {best[0][0]} (target {target})", flush=True)
        for tot, spec in best[:20]:
            n4, n5 = census_pts(mixed_points(spec, grid))
            if n5 == 0 and n4 > target:
                out.append({'moduli': moduli, 'n': n, 'n4': n4, 'spec': [[m, r, p] for m, r, p in spec], 'target': target})
                if verbose: print(f"    MIXED CANDIDATE n={n} n4={n4} > {target} : {spec}", flush=True)
    return out

if __name__ == '__main__':
    records = {r['N']: r['score'] for r in json.load(open('final_rows.json'))}
    ints = {r['n']: r['score'] for r in json.load(open('records_int_now.json'))['records']}
    mode = sys.argv[1] if len(sys.argv) > 1 else 'A'
    res = []
    if mode == 'A':
        for m, k in [(7,4),(9,4),(10,4),(12,4),(8,4),(11,4),(13,4),(14,4),(15,4),(6,4),(5,4)]:
            t = time.time()
            r = partA(m, k, records, ints)
            res.append(r)
            print(f"m={m} k={k} n={r['n']}: max(2,2)={r['max_pairbound']}  best valid n4={r['best_valid_n4']}  target={r['target']}  ({time.time()-t:.1f}s)", flush=True)
        json.dump(res, open('ring_optimality.json', 'w'), indent=1)
    else:
        t0 = time.time()
        budget = float(sys.argv[2]) if len(sys.argv) > 2 else 600
        mods = []
        for a in range(5, 21):
            for b in range(a, 21):
                if 24 <= a+b <= 60: mods.append([a, b])
        for a in range(5, 16):
            for b in range(a, 16):
                for c in range(b, 16):
                    if 30 <= a+b+c <= 60: mods.append([a, b, c])
        mods = [mm for mm in mods if math.lcm(*mm) <= 60]
        mods.sort(key=lambda L: (-sum(L), math.lcm(*L)))
        for mm in mods:
            if time.time()-t0 > budget: print("budget reached", flush=True); break
            try: res += partB(mm, records, ints, deadline=t0+budget)
            except Exception as e: print("  err", mm, e, flush=True)
        json.dump(res, open('ring_mixed_results.json', 'w'), indent=1)
    print("IMPROVEMENTS:", [(o.get('n'), o.get('n4')) for o in res if o.get('n4', 0) > o.get('target', 10**9)], flush=True)
