"""Turn a float symmetric ansatz into an exact orchard over Q(cos(pi/(2m))).

All angles are multiples of pi/(2m); cos(k pi/(2m)) = T_k(c) with c = cos(pi/(2m)).
Radii are matched to exact cosine ratios, or recognized by PSLQ in Q(c).
"""
import numpy as np
import sympy as sp
import mpmath as mp
from fractions import Fraction

from algfield import Field
from orchard import analyze, honesty, entry_text

mp.mp.dps = 60


def base_field(m):
    c = sp.cos(sp.pi / (2 * m))
    mpoly = sp.minpoly(c, x=sp.Symbol('t'))
    coeffs = [int(x) for x in sp.Poly(mpoly, sp.Symbol('t')).all_coeffs()][::-1]
    lc = coeffs[-1]
    coeffs = [Fraction(x, lc) for x in coeffs]
    # Field requires monic integer coefficients; scale variable instead if lc != 1
    if lc != 1:
        # use generator 2c instead: minpoly of 2c is monic
        c2 = 2 * c
        mpoly = sp.minpoly(c2, x=sp.Symbol('t'))
        coeffs = [int(x) for x in sp.Poly(mpoly, sp.Symbol('t')).all_coeffs()][::-1]
        assert coeffs[-1] == 1
        gen = c2
        scale = Fraction(1, 2)  # c = gen * scale
    else:
        gen = c
        scale = Fraction(1)
    F = Field(coeffs)
    return F, gen, scale, coeffs


class Exactifier:
    def __init__(self, m):
        self.m = m
        F, gen, scale, coeffs = base_field(m)
        self.F = F
        self.coeffs = coeffs
        self.scale = scale  # c = scale * A
        self.c = F.scale(F.A, scale)
        self.cnum = mp.cos(mp.pi / (2 * m))
        # Chebyshev table T_k(c), k = 0 .. 4m+2
        self.T = [F.ONE, self.c]
        for k in range(2, 4 * m + 3):
            self.T.append(F.sub(F.scale(F.mul(self.c, self.T[-1]), 2), self.T[-2]))

    def point(self, o_phase, i, r):
        """exact point: r * (cos((4i+2ph) pi/(2m)), sin(...)); r in F."""
        k = 4 * i + 2 * o_phase
        x = self.F.mul(r, self.T[k])
        y = self.F.mul(r, self.T[abs(self.m - k)])
        return (x, y)

    def ratio_match(self, r_target, phases, a, b, tol=None):
        """match r_b/r_a to +-cos(pi j/m)/cos(pi d/m) with direction condition."""
        m = self.m
        if tol is None:
            tol = mp.mpf('1e-9')
        for j in range(1, m):
            if 2 * j == m:
                continue
            for l in range(m):
                for h in range(m):
                    if l == h:
                        continue
                    if (j - (phases[b] - phases[a]) - l - h) % m != 0:
                        continue
                    t = (phases[a] + j) - (phases[b] + l + h)
                    sgn = -1 if (t // m) % 2 else 1
                    d = h - l
                    cd = mp.cos(mp.pi * d / m)
                    if abs(cd) < mp.mpf('1e-50'):
                        continue
                    rv = sgn * mp.cos(mp.pi * j / m) / cd
                    if rv <= 0:
                        continue
                    if abs(rv - r_target) < tol:
                        # exact: cos(pi j/m) = T_{2j}(c), cos(pi d/m) = T_{2|d|}(c)
                        num = self.T[2 * j]
                        den = self.T[2 * abs(d)]
                        rex = self.F.div(num, den)
                        if sgn < 0:
                            rex = self.F.neg(rex)
                        # verify numerically to 40 digits
                        if abs(self.evalf(rex) - r_target) < tol:
                            return rex
        return None

    def evalf(self, u):
        # evaluate field element at the generator A = c/scale numerically
        cv = self.cnum * mp.mpf(self.scale.denominator) / self.scale.numerator
        r = mp.mpf(0)
        for co in reversed(u):
            r = r * cv + mp.mpf(co.numerator) / co.denominator
        return r

    def pslq_recognize(self, r_target):
        """recognize r in Q(c) via PSLQ over powers of c."""
        D = self.F.d
        cv = self.cnum
        basis = [r_target] + [cv ** i for i in range(D)]
        rel = mp.pslq(basis, tol=mp.mpf('1e-45'), maxcoeff=10 ** 8)
        if rel is None or rel[0] == 0:
            return None
        # rel[0]*r + sum rel[i+1] c^i = 0
        rex = self.F.ZERO
        for i in range(D):
            if rel[i + 1]:
                rex = self.F.add(rex, self.F.scale(self.T[i], Fraction(-rel[i + 1], rel[0])))
        if abs(self.evalf(rex) - r_target) < mp.mpf('1e-35'):
            return rex
        return None

    def exactify(self, phases, radii, center, refined=None):
        """radii: floats (r0=1). refined: optional dict orbit->60-digit mp value."""
        F = self.F
        rex = [F.ONE]
        for o in range(1, len(radii)):
            target = mp.mpf(refined[o]) if refined and o in refined else mp.mpf(radii[o])
            r = None
            for a in range(o):
                # ratio r_o/r_a; r_a exact already
                rv = target / self.evalf(rex[a])
                r = self.ratio_match(rv, phases, a, o)
                if r is not None:
                    r = F.mul(r, rex[a])
                    break
            if r is None:
                r = self.pslq_recognize(target)
            if r is None:
                return None, f"radius {o} = {radii[o]} not recognized"
            rex.append(r)
        pts = []
        for o, ph in enumerate(phases):
            for i in range(self.m):
                pts.append(self.point(ph, i, rex[o]))
        if center:
            pts.append((F.ZERO, F.ZERO))
        return pts, rex


def refine_radii(m, phases, tup, x0):
    """Refine a k=3 (2,1,1) solution to 60 digits with mpmath."""
    mp.mp.dps = 60
    ang = [[2 * mp.pi * i / m + ph * mp.pi / m for i in range(m)] for ph in phases]
    cs = [[mp.cos(a) for a in row] for row in ang]
    sn = [[mp.sin(a) for a in row] for row in ang]
    (o1, i1), (o2, i2), (o3, i3), (o4, i4) = tup

    def pts(rv):
        r = [mp.mpf(1), rv[0], rv[1]]
        return [(r[o] * cs[o][i], r[o] * sn[o][i]) for o, i in [(o1, i1), (o2, i2), (o3, i3), (o4, i4)]]

    def f(r1, r2):
        p1, p2, p3, p4 = pts([r1, r2])
        d1 = (p2[0] - p1[0]) * (p3[1] - p1[1]) - (p2[1] - p1[1]) * (p3[0] - p1[0])
        d2 = (p2[0] - p1[0]) * (p4[1] - p1[1]) - (p2[1] - p1[1]) * (p4[0] - p1[0])
        return (d1, d2)

    try:
        sol = mp.findroot(f, (mp.mpf(x0[0]), mp.mpf(x0[1])), tol=mp.mpf('1e-55'), maxsteps=100)
        return {1: sol[0], 2: sol[1]}
    except Exception as e:
        return None


def process_result(r, verbose=True):
    m, k, phases, center, radii = r["m"], r["k"], tuple(r["phases"]), r["center"], r["radii"]
    ex = Exactifier(m)
    refined = None
    if r.get("tuple"):
        refined = refine_radii(m, phases, [tuple(t) for t in r["tuple"]], radii[1:])
    pts, info = ex.exactify(phases, radii, center, refined)
    if pts is None:
        if verbose:
            print(f"N={r['N']} m={m}: FAILED {info}")
        return None
    res = analyze(ex.F, pts)
    h = honesty(ex.F, pts, res)
    if verbose:
        print(f"N={len(pts)} m={m} k={k} ph={phases} c={center}: float {r['score']} -> exact {res['score']} valid={res['valid']} irr_cr={len(h['irrational_cross_ratios'])} field deg={ex.F.d}")
    if not res["valid"] or res["score"] < r["score"]:
        if verbose:
            print("   mismatch!", res.get("reason", ""))
        return None
    txt = entry_text(ex.F, pts, ex.coeffs)
    return {"N": len(pts), "score": res["score"], "text": txt, "meta": {**r, "field_deg": ex.F.d, "radii_exact": [ex.F.fmt(u) for u in info]}}
