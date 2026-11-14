"""Exact orchard scoring over a number field, plus float search utilities."""
from fractions import Fraction
from itertools import combinations
import numpy as np


# ---------------- exact ----------------

def line_key(F, p, q):
    """Normalized (A,B,C) with A x + B y + C = 0 through p, q (exact)."""
    x1, y1 = p
    x2, y2 = q
    A = F.sub(y1, y2)
    B = F.sub(x2, x1)
    C = F.sub(F.mul(x1, y2), F.mul(x2, y1))
    for c in (A, B, C):
        if not F.is_zero(c):
            iv = F.inv(c)
            return (F.mul(A, iv), F.mul(B, iv), F.mul(C, iv))
    raise ValueError("identical points")


def analyze(F, pts, check=True):
    """Return dict: lines (frozenset -> size), score, valid, reasons."""
    n = len(pts)
    for i, j in combinations(range(n), 2):
        if pts[i] == pts[j]:
            return {"valid": False, "reason": f"duplicate points {i},{j}", "score": 0}
    pairlines = {}
    for i, j in combinations(range(n), 2):
        k = line_key(F, pts[i], pts[j])
        pairlines.setdefault(k, []).append((i, j))
    lines = {}
    for k, prs in pairlines.items():
        m = len(prs)
        # m = C(s,2): s = (1+sqrt(1+8m))/2
        s = int((1 + (1 + 8 * m) ** 0.5) / 2 + 1e-9)
        assert s * (s - 1) // 2 == m
        verts = frozenset(x for pr in prs for x in pr)
        assert len(verts) == s
        lines[verts] = k
    sizes = {}
    for v in lines:
        sizes.setdefault(len(v), []).append(v)
    bad = [v for v in lines if len(v) >= 5]
    score = len(sizes.get(4, []))
    return {
        "valid": not bad,
        "reason": ("lines with >=5 points: " + str([sorted(v) for v in bad])) if bad else "",
        "score": score,
        "lines": lines,
        "sizes": sizes,
    }


def cross_ratio(F, pts, quad):
    """Cross-ratio of 4 collinear points (projective invariant)."""
    quad = sorted(quad)
    p1, p2, p3, p4 = [pts[i] for i in quad]
    # parameter: use x unless all equal, then y
    t = [p[0] for p in (p1, p2, p3, p4)]
    if len(set(t)) < 4:
        t = [p[1] for p in (p1, p2, p3, p4)]
    t1, t2, t3, t4 = t
    num = F.mul(F.sub(t3, t1), F.sub(t4, t2))
    den = F.mul(F.sub(t3, t2), F.sub(t4, t1))
    return F.div(num, den)


def honesty(F, pts, res):
    certs = []
    for quad in res["sizes"].get(4, []):
        cr = cross_ratio(F, pts, quad)
        if not F.is_rational(cr):
            certs.append((sorted(quad), cr))
    all_rational = all(F.is_rational(c) for p in pts for c in p)
    return {"all_coords_rational": all_rational, "irrational_cross_ratios": certs}


def entry_text(F, pts, poly_coeffs):
    x = 'a'
    terms = []
    d = len(poly_coeffs) - 1
    for i in range(d, -1, -1):
        c = poly_coeffs[i]
        if c == 0:
            continue
        sign = '-' if c < 0 else '+'
        ac = abs(c)
        if i == 0:
            t = str(ac)
        else:
            xs = x if i == 1 else f"{x}^{i}"
            t = xs if ac == 1 else f"{ac}{xs}"
        if not terms:
            terms.append(('-' if c < 0 else '') + t)
        else:
            terms.append(f" {sign} {t}")
    poly = ''.join(terms) + " = 0"
    body = ", ".join(f"({F.fmt(p[0])}, {F.fmt(p[1])})" for p in pts)
    return poly + "\n" + body


# ---------------- float (search) ----------------

def float_lines(P, tol=1e-7):
    """P: (n,2) float array. Returns dict size -> list of tuples of indices,
    counting lines with >= 3 points. tol relative to coordinate scale."""
    P = np.asarray(P, float)
    n = len(P)
    scale = max(1.0, np.abs(P).max())
    found = {}
    for i in range(n):
        for j in range(i + 1, n):
            d = P[j] - P[i]
            L = np.hypot(*d)
            if L < 1e-12:
                continue
            cr = np.abs((P[:, 0] - P[i, 0]) * d[1] - (P[:, 1] - P[i, 1]) * d[0]) / L
            on = np.where(cr < tol * scale)[0]
            if len(on) >= 3:
                key = tuple(on)
                found[key] = len(on)
    sizes = {}
    for k, s in found.items():
        sizes.setdefault(s, []).append(k)
    return sizes


def float_score(P, tol=1e-7):
    sizes = float_lines(P, tol)
    if any(s >= 5 for s in sizes):
        return -1, sizes
    return len(sizes.get(4, [])), sizes
