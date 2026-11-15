"""Exact arithmetic in Q[a]/p(a) for the algebraic orchard problem.

Elements are tuples of Fractions (coefficients of 1, a, a^2, ...).
p must be monic with integer coefficients; irreducibility is checked.
"""
from fractions import Fraction
import sympy as sp


class Field:
    def __init__(self, poly):
        # poly: list of integer coefficients, low degree first, monic
        self.p = [Fraction(c) for c in poly]
        assert self.p[-1] == 1
        self.d = len(poly) - 1
        if self.d == 1:
            self.p = [Fraction(0), Fraction(1)]  # a = 0, rational field
        x = sp.Symbol('x')
        expr = sum(int(c) * x**i for i, c in enumerate(poly))
        if self.d >= 1 and not sp.Poly(expr, x).is_irreducible:
            raise ValueError(f"polynomial {poly} is not irreducible")
        # reduction table: a^k for k = d .. 2d-2 as coefficient vectors
        self.red = []
        for k in range(self.d, 2 * self.d - 1):
            v = [Fraction(0)] * (k + 1)
            v[k] = Fraction(1)
            for j in range(k, self.d - 1, -1):
                c = v[j]
                if c:
                    v[j] = Fraction(0)
                    for i in range(self.d):
                        v[j - self.d + i] -= c * self.p[i]
            self.red.append(v[:self.d])
        self.ZERO = (Fraction(0),) * self.d
        self.ONE = (Fraction(1),) + (Fraction(0),) * (self.d - 1)
        self.A = (Fraction(0), Fraction(1)) + (Fraction(0),) * (self.d - 2) if self.d >= 2 else self.ZERO

    def elem(self, c):
        if isinstance(c, tuple):
            v = list(c) + [Fraction(0)] * (self.d - len(c))
            return tuple(Fraction(x) for x in v[:self.d])
        if isinstance(c, (int, Fraction)):
            return (Fraction(c),) + (Fraction(0),) * (self.d - 1)
        raise TypeError(c)

    def is_zero(self, u):
        return not any(u)

    def is_rational(self, u):
        return not any(u[1:])

    def add(self, u, v):
        return tuple(x + y for x, y in zip(u, v))

    def neg(self, u):
        return tuple(-x for x in u)

    def sub(self, u, v):
        return tuple(x - y for x, y in zip(u, v))

    def mul(self, u, v):
        d = self.d
        w = [Fraction(0)] * (2 * d - 1)
        for i, ui in enumerate(u):
            if ui:
                for j, vj in enumerate(v):
                    if vj:
                        w[i + j] += ui * vj
        for k in range(2 * d - 2, d - 1, -1):
            c = w[k]
            if c:
                w[k] = Fraction(0)
                rk = self.red[k - d]
                for i in range(d):
                    w[i] += c * rk[i]
        return tuple(w[:d])

    def scale(self, u, c):
        c = Fraction(c)
        return tuple(x * c for x in u)

    def inv(self, u):
        # solve v*u = 1 mod p by linear algebra over Q (d small)
        d = self.d
        if self.is_zero(u):
            raise ZeroDivisionError
        if d == 1:
            return (1 / u[0],)
        # matrix M: column j = u * a^j mod p
        M = []
        v = list(u)
        for j in range(d):
            M.append(list(v))
            v = list(self.mul(tuple(v), self.A))
        # Gaussian elimination on [M | e1]
        rows = [[M[j][i] for j in range(d)] + [Fraction(1 if i == 0 else 0)] for i in range(d)]
        for col in range(d):
            piv = next((r for r in range(col, d) if rows[r][col] != 0), None)
            if piv is None:
                raise ZeroDivisionError("non-invertible element")
            rows[col], rows[piv] = rows[piv], rows[col]
            pv = rows[col][col]
            rows[col] = [x / pv for x in rows[col]]
            for r in range(d):
                if r != col and rows[r][col] != 0:
                    f = rows[r][col]
                    rows[r] = [a - f * b for a, b in zip(rows[r], rows[col])]
        return tuple(rows[i][d] for i in range(d))

    def div(self, u, v):
        return self.mul(u, self.inv(v))

    def pow(self, u, n):
        if n < 0:
            return self.pow(self.inv(u), -n)
        r = self.ONE
        while n:
            if n & 1:
                r = self.mul(r, u)
            u = self.mul(u, u)
            n >>= 1
        return r

    def to_float(self, u, root):
        # Horner at a real root
        r = 0.0
        for c in reversed(u):
            r = r * root + float(c)
        return r

    def fmt(self, u, var='a'):
        """Pretty polynomial string with fraction coefficients."""
        terms = []
        for i in range(self.d - 1, -1, -1):
            c = u[i]
            if c == 0:
                continue
            sign = '-' if c < 0 else '+'
            ac = abs(c)
            if i == 0:
                cs = str(ac.numerator) if ac.denominator == 1 else f"{ac.numerator}/{ac.denominator}"
                t = cs
            else:
                xs = var if i == 1 else f"{var}^{i}"
                if ac.numerator == 1:
                    t = xs if ac.denominator == 1 else f"{xs}/{ac.denominator}"
                elif ac.denominator == 1:
                    t = f"{ac.numerator}{xs}"
                else:
                    t = f"{ac.numerator}{xs}/{ac.denominator}"
            terms.append((sign, t))
        if not terms:
            return "0"
        out = ""
        for j, (sign, t) in enumerate(terms):
            if j == 0:
                out += ("-" if sign == '-' else "") + t
            else:
                out += f" {sign} {t}"
        return out


def field_from_string(s):
    """Parse 'a^3 + a^2 - 2a - 1 = 0' into a Field."""
    s = s.strip()
    if '=' in s:
        s = s.split('=')[0]
    x = sp.Symbol('a')
    from sympy.parsing.sympy_parser import parse_expr as sp_parse, standard_transformations, implicit_multiplication_application, convert_xor
    expr = sp_parse(s, local_dict={'a': x},
                    transformations=standard_transformations + (implicit_multiplication_application, convert_xor))
    poly = sp.Poly(sp.expand(expr), x)
    coeffs = [int(c) for c in poly.all_coeffs()][::-1]
    assert poly.LC() == 1, "polynomial must be monic"
    return Field(coeffs)


# ---- expression parser: rational expressions in a, evaluated in a Field ----
import re

TOKEN = re.compile(r'\s*(\d+|a|\^|\+|-|\*|/|\(|\))')


def parse_expr(s, F):
    toks = []
    pos = 0
    while pos < len(s):
        m = TOKEN.match(s, pos)
        if not m:
            if s[pos:].strip() == '':
                break
            raise ValueError(f"bad token at {s[pos:]}")
        toks.append(m.group(1))
        pos = m.end()
    it = iter(toks)
    tok = next(it, None)

    def peek():
        return tok

    def eat(t=None):
        nonlocal tok
        if t is not None and tok != t:
            raise ValueError(f"expected {t} got {tok}")
        cur = tok
        tok = next(it, None)
        return cur

    def atom():
        t = peek()
        if t == '(':
            eat('(')
            v = expr()
            eat(')')
            return v
        if t == 'a':
            eat('a')
            return F.A
        if t is not None and re.fullmatch(r'\d+', t):
            eat()
            return F.elem(int(t))
        raise ValueError(f"unexpected token {t}")

    def power():
        v = atom()
        while peek() == '^':
            eat('^')
            e = eat()
            v = F.pow(v, int(e))
        return v

    def unary():
        if peek() == '-':
            eat('-')
            return F.neg(unary())
        if peek() == '+':
            eat('+')
            return unary()
        return power()

    def term():
        v = unary()
        while True:
            t = peek()
            if t == '*':
                eat('*')
                v = F.mul(v, unary())
            elif t == '/':
                eat('/')
                v = F.div(v, unary())
            elif t in ('a', '(') or (t is not None and re.fullmatch(r'\d+', t)):
                v = F.mul(v, unary())  # implicit multiplication
            else:
                return v

    def expr():
        v = term()
        while peek() in ('+', '-'):
            op = eat()
            w = term()
            v = F.add(v, w) if op == '+' else F.sub(v, w)
        return v

    v = expr()
    if tok is not None:
        raise ValueError(f"trailing tokens: {tok}")
    return v


def parse_entry(text):
    """Parse a full submission: polynomial line, then points."""
    lines = [l for l in text.strip().splitlines() if l.strip()]
    F = field_from_string(lines[0])
    body = ' '.join(lines[1:])
    pts = []
    for m in re.finditer(r'\(([^()]*(?:\([^()]*\)[^()]*)*)\)', body):
        inner = m.group(1)
        # split on the top-level comma
        depth = 0
        for i, ch in enumerate(inner):
            if ch == '(':
                depth += 1
            elif ch == ')':
                depth -= 1
            elif ch == ',' and depth == 0:
                x = parse_expr(inner[:i], F)
                y = parse_expr(inner[i + 1:], F)
                pts.append((x, y))
                break
        else:
            raise ValueError(f"point without comma: ({inner})")
    return F, pts
