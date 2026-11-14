#!/usr/bin/env python3
"""Independent exact verifier for algebraic orchards (Infinite Search Space,
Orchard Planting - algebraic variant).

Usage:
    python3 verify.py entries/N60.txt          # verify one entry
    python3 verify.py entries/*.txt            # verify many
    cat entry.txt | python3 verify.py -        # read stdin

Checks, with EXACT arithmetic in Q[a]/p(a) (no rounding):
  * p is monic, integer, irreducible (so every root gives the same lines)
  * no two trees coincide
  * score = number of lines through exactly four trees
  * invalid if any line passes through >= 5 trees
  * honesty: not all coordinates rational, and at least one 4-tree line has an
    irrational cross-ratio (a projective invariant) -> provably NOT a stretched
    / tilted whole-number orchard.
Requires: sympy  (pip install sympy)
"""
import sys
from algfield import parse_entry
from orchard import analyze, honesty


def verify(text):
    F, pts = parse_entry(text)
    res = analyze(F, pts)
    out = {"N": len(pts), "field_degree": F.d, "valid": res["valid"]}
    if not res["valid"]:
        out["reason"] = res["reason"]
        out["score"] = 0
        return out
    h = honesty(F, pts, res)
    out["score"] = res["score"]
    out["all_coords_rational"] = h["all_coords_rational"]
    out["lines_with_irrational_cross_ratio"] = len(h["irrational_cross_ratios"])
    out["honest_algebraic"] = (not h["all_coords_rational"]) and bool(h["irrational_cross_ratios"])
    out["line_size_histogram"] = {k: len(v) for k, v in sorted(res["sizes"].items())}
    return out


def main():
    args = sys.argv[1:]
    if not args:
        print(__doc__)
        return
    texts = []
    if args == ["-"]:
        texts = [("<stdin>", sys.stdin.read())]
    else:
        for f in args:
            texts.append((f, open(f).read()))
    ok = True
    for name, text in texts:
        try:
            r = verify(text)
        except Exception as e:
            print(f"{name}: PARSE/VERIFY ERROR: {e}")
            ok = False
            continue
        flag = "OK " if (r["valid"] and r.get("honest_algebraic")) else "!! "
        if not (r["valid"] and r.get("honest_algebraic")):
            ok = False
        print(f"{flag}{name}: N={r['N']} score={r['score']} valid={r['valid']} "
              f"honest={r.get('honest_algebraic')} deg={r['field_degree']} "
              f"hist={r.get('line_size_histogram')}")
        if not r["valid"]:
            print("   reason:", r.get("reason"))
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
