"""Grow orchards exactly: add points at intersections of two 3-point lines.
Also try projective re-views (bringing directions / points at infinity into
the affine plane) before extending, and delete+regrow for other sizes.
All arithmetic stays in the same number field.
"""
import numpy as np
from itertools import combinations
from fractions import Fraction

from orchard import analyze, honesty, line_key


def to_float_pts(F, pts, root):
    return np.array([[F.to_float(x, root), F.to_float(y, root)] for x, y in pts])


def lines_of(F, pts):
    """exact lines with >=2 points: list (frozenset, (A,B,C))."""
    res = analyze(F, pts)
    if not res["valid"]:
        return None, res
    out = [(v, k) for v, k in res["lines"].items()]
    return out, res


def float_candidates(P, min_gain=1):
    """Float pass: intersections of pairs of 3-point lines; gain = #3-lines through.
    Returns list of (point, gain) sorted by gain desc, deduped."""
    n = len(P)
    scale = max(1.0, np.abs(P).max())
    tol = 1e-7 * scale
    # collect 3-point lines (as index tuples + line coeffs)
    lines3 = []
    lines4 = []
    seen = set()
    for i in range(n):
        for j in range(i + 1, n):
            d = P[j] - P[i]
            L = np.hypot(*d)
            if L < 1e-12:
                continue
            cr = np.abs((P[:, 0] - P[i, 0]) * d[1] - (P[:, 1] - P[i, 1]) * d[0]) / L
            on = tuple(np.where(cr < tol)[0])
            if on in seen:
                continue
            seen.add(on)
            # line coeffs A x + B y + C = 0, normalized
            A, B = d[1] / L, -d[0] / L
            C = -(A * P[i, 0] + B * P[i, 1])
            if len(on) == 3:
                lines3.append((on, (A, B, C)))
            elif len(on) >= 4:
                lines4.append((on, (A, B, C)))
    cands = []
    for (o1, l1), (o2, l2) in combinations(lines3, 2):
        if len(set(o1) & set(o2)) >= 1:
            continue  # meet at an existing point
        A1, B1, C1 = l1
        A2, B2, C2 = l2
        det = A1 * B2 - A2 * B1
        if abs(det) < 1e-12:
            continue
        x = (B1 * C2 - B2 * C1) / det
        y = (C1 * A2 - C2 * A1) / det
        p = np.array([x, y])
        if np.abs(p).max() > 1e6:
            continue
        # gain: number of 3-lines through p; invalid if on a 4-line
        gain = 0
        bad = False
        for on, (A, B, C) in lines3:
            if abs(A * x + B * y + C) < tol:
                gain += 1
        for on, (A, B, C) in lines4:
            if abs(A * x + B * y + C) < tol:
                bad = True
                break
        if bad or gain < min_gain:
            continue
        if np.min(np.linalg.norm(P - p, axis=1)) < 1e-8:
            continue
        cands.append((p, gain))
    cands.sort(key=lambda t: -t[1])
    # dedupe
    out = []
    for p, g in cands:
        if all(np.linalg.norm(p - q) > 1e-7 for q, _ in out):
            out.append((p, g))
    return out


def exact_point_from_lines(F, l1, l2):
    A1, B1, C1 = l1
    A2, B2, C2 = l2
    den = F.sub(F.mul(A1, B2), F.mul(A2, B1))
    if F.is_zero(den):
        return None
    x = F.div(F.sub(F.mul(B1, C2), F.mul(B2, C1)), den)
    y = F.div(F.sub(F.mul(C1, A2), F.mul(C2, A1)), den)
    return (x, y)


def grow(F, pts, root, max_add=10, verbose=False):
    """Greedy exact growth. Returns (pts, score, added)."""
    added = 0
    cur = analyze(F, pts)
    if not cur["valid"]:
        return pts, -1, 0
    while added < max_add:
        P = to_float_pts(F, pts, root)
        cands = float_candidates(P, min_gain=2)
        if not cands:
            cands = float_candidates(P, min_gain=1)
        if not cands:
            break
        lines, res = lines_of(F, pts)
        l3 = [(v, k) for v, k in lines if len(v) == 3]
        l4 = [(v, k) for v, k in lines if len(v) >= 4]
        best = None
        for p, g in cands[:40]:
            # find exact pairs of 3-lines through p (float match then exact)
            through = []
            for v, k in l3:
                vv = sorted(v)
                # distance test in float
                A, B, C = [F.to_float(c, root) for c in k]
                if abs(A * p[0] + B * p[1] + C) < 1e-6 * max(1, abs(A), abs(B)):
                    through.append((v, k))
            if len(through) < 2:
                continue
            q = None
            for (v1, k1), (v2, k2) in combinations(through, 2):
                q = exact_point_from_lines(F, k1, k2)
                if q is not None:
                    break
            if q is None:
                continue
            if q in pts:
                continue
            # exact gain / validity
            gain = 0
            bad = False
            for v, k in lines:
                A, B, C = k
                val = F.add(F.add(F.mul(A, q[0]), F.mul(B, q[1])), C)
                if F.is_zero(val):
                    if len(v) == 3:
                        gain += 1
                    elif len(v) >= 4:
                        bad = True
                        break
            if bad or gain < 1:
                continue
            if best is None or gain > best[1]:
                best = (q, gain)
        if best is None:
            break
        pts = pts + [best[0]]
        added += 1
        if verbose:
            print(f"  +point gain {best[1]} -> N={len(pts)} score={cur['score'] + best[1]}", flush=True)
        cur = analyze(F, pts)
        if not cur["valid"]:
            raise RuntimeError("grew invalid?!")
    return pts, cur["score"], added


def projective_view(F, pts, t, swap=False):
    """(x:y:1) -> (x:y:1+t x); optional coordinate swap first. Returns None if a
    point goes to infinity."""
    out = []
    for x, y in pts:
        if swap:
            x, y = y, x
        den = F.add(F.ONE, F.mul(t, x))
        if F.is_zero(den):
            return None
        out.append((F.div(x, den), F.div(y, den)))
    return out


def grow_multiview(F, pts, root, ts=None, max_add=10, verbose=False):
    """Try growth in several projective views; return best (pts, score, view)."""
    if ts is None:
        ts = [F.elem(1), F.elem(-1), F.elem(2), F.elem(Fraction(1, 2)),
              F.A, F.neg(F.A), F.add(F.A, F.ONE)]
    best = (pts, analyze(F, pts)["score"], None)
    for swap in (False, True):
        for t in ts:
            v = projective_view(F, pts, t, swap)
            if v is None:
                continue
            try:
                g, sc, ad = grow(F, v, root, max_add=max_add, verbose=verbose)
            except Exception:
                continue
            if sc > best[1]:
                best = (g, sc, (swap, t))
    return best
