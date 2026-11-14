> **Historical research note, not the current manuscript.** Some claims and record comparisons below are unverified or superseded, including statements of family optimality, completeness, novelty, and exhaustive global or local search. See [the revised paper](../paper/orchard.pdf) and [research notes](../research/README.md) for the corrected status.

# Literature audit, upper bounds, and the “at least four” variant

**Date:** 2026-11-14  
**Scope:** real planar orchards; contest variant = exactly 4 points per scored line and **no line with ≥5**; “at least four” variant = count every line with **≥4** points, long lines allowed.

Notation (Elkies):
- `t_4(n)` = max number of lines containing **exactly** 4 of `n` points, with no restriction on longer lines.
- `t_4^{(5)}(n)` = same, but **no 5 collinear**. This is the contest quantity.
- `T_4(n)` = max number of lines containing **at least** 4 points.

Clearly `t_4^{(5)}(n) ≤ t_4(n) ≤ T_4(n)`.

---

## 1. Novelty audit: what is already known

### 1.1 The six candidate algebraic records
Our six candidate no-five lower bounds are still:

| n | score | construction |
|---|---:|---|
| 28 | 42 | k4-ratio, deg 6 |
| 36 | 63 | celestial-k4, deg 6 |
| 37 | 63 | celestial-k4 + 1 |
| 40 | 70 | celestial-k4, deg 8 |
| 48 | 96 | celestial-k4, deg 8 |
| 49 | 97 | celestial-k4 + 1 |

Checked against:
- **OEIS A006065**: real exactly-4 sequence only to `n=20`; values match known optima there.
- **OEIS A172992**: integer-coordinate lower bounds to `n=60`; all six exceed it.
- **Friedman / MathWorld table**: bounds only to `n=25`; all six are beyond that range.
- **Du**: optimal catalogues for 13–18, 20-tree problem, and 22:28; no larger scores found.
- **Pegg / Wolfram blog**: 3-orchard methods and 4-orchard examples up to 25; no larger scores found.
- **Elkies math/0612749**: points-and-lines configurations, but aimed at `n^2` points with `n`-point lines and Brass transversals; not these scores.
- **Berman–Gévay–Richter-Gebert–Tabachnikov, arXiv:2408.09203**: “trivial celestial” `(n_4)` configurations are balanced point-line configurations (four lines through each point, four points on each line). They do **not** tabulate orchard scores of the form `score > n`, and the paper does not mention the orchard problem.

**Conclusion:** the six scores remain plausible new lower bounds for the no-five real orchard problem. They are now on the algebraic leaderboard, so priority should be tied to the submission timestamp. A final expert check (Pegg / Kamenetsky / Du / Grünbaum-school configuration literature) is still advisable before an arXiv claim.

### 1.2 Correction to the previous report
The general upper bound

\[
t_4(n) \le \left\lfloor \frac{n(n-1)-6}{14} \right\rfloor
\]

is **not novel**: it is already stated on MathWorld’s *Orchard-Planting Problem* page (Eq. 4), derived from Melchior’s inequality and pair-counting. MathWorld also records `t_4(21) ≤ 28` by deletion from `t_4(20)=23`.

What *may* be new is packaging this for the **at-least-four** variant `T_4`, and a few small corollaries below.

---

## 2. Upper bounds we can prove cleanly

### Proposition A (Melchior bound for exactly-four and at-least-four)
For `n ≥ 5`,

\[
t_4(n) \le T_4(n) \le \left\lfloor \frac{n(n-1)-6}{14} \right\rfloor .
\]

**Proof.** If all points are collinear then `T_4(n)=1` and the bound is trivial for `n≥5`. Otherwise let `t_j` be the number of lines containing exactly `j` points. Melchior’s inequality gives

\[
t_2 \ge 3 + \sum_{j\ge4}(j-3)t_j .
\]

Pair-counting gives

\[
\binom n2 = \sum_{j\ge2}\binom j2 t_j .
\]

For `T_4 = Σ_{j≥4} t_j`,

\[
\binom n2
\ge t_2 + \sum_{j\ge4}\binom j2 t_j
\ge 3 + \sum_{j\ge4}\left(\binom j2 + j - 3\right)t_j .
\]

For `j≥4`, `binom(j,2)+j-3 ≥ 7`, with equality only at `j=4`. Hence

\[
\binom n2 \ge 3 + 7T_4,
\]

which rearranges to the claim. The same proof with only `t_4` gives the exactly-four bound. ∎

**Status:** the exactly-four case is MathWorld Eq. 4. The at-least-four formulation is a direct corollary; I did not find it explicitly tabulated in the orchard sources checked, but it should be regarded as “Melchior folklore” unless a prior reference is located.

### Proposition B (exact slack identity)
Let `F ≥ 0` be the Melchior face-slack term. Then

\[
\binom n2 - 3
= 7T_4 + 3t_3 + F
  + \sum_{j\ge5}\left(\binom j2 + j - 3 - 7\right)t_j .
\]

The coefficients for `j=5,6,7,…` are `5,11,18,…`.

**Use:** any triple line, non-simplicial face, or long line strictly reduces the possible `T_4` below the Melchior ceiling. This is useful for ruling out near-extremal combinatorics.

### Proposition C (a concrete at-least-four upper bound: `T_4(12) ≤ 8`)
Melchior gives `T_4(12) ≤ 9`. If `T_4(12)=9`, equality in Proposition A forces

\[
t_4=9,\quad t_2=12,\quad t_3=0,\quad t_{\ge5}=0,\quad F=0.
\]

Dualizing, this is an arrangement of 12 pseudolines with 9 quadruple points and 12 double points. Bokowski–Pokora (arXiv:1607.05864, Cor. 1.7) prove that **no such pseudoline arrangement exists**. Therefore

\[
T_4(12) \le 8.
\]

This is a small but concrete upper bound for the at-least-four variant that is not in the MathWorld exactly-four table. It is a candidate “novel corollary” for a note.

### Proposition D (degree-convexity bound for the no-five contest)
For the contest quantity `L = t_4^{(5)}(n)`, let `r_i` be the number of scored 4-lines through point `i`. Then

\[
0 \le r_i \le \left\lfloor\frac{n-1}{3}\right\rfloor,\qquad
\sum_i r_i = 4L,\qquad
\sum_i \binom{r_i}{2} \le \binom{L}{2}.
\]

Convexity gives, for `n<16`,

\[
L \le \left\lfloor \frac{3n}{16-n} \right\rfloor .
\]

In particular `L(11) ≤ 6`, matching our construction. This proof is already in `N11_REPORT.md`; it is not new as a value, but it is a clean general inequality.

### Selected bound table
Full table: `bounds_table.md`.

| n | no-5 lower (best board) | no-5 upper | at-least-4 lower | at-least-4 upper |
|---:|---:|---:|---:|---:|
| 11 | 6 | 6 | 6 | 7 |
| 12 | 7 | 7 | 7 | **8** |
| 19 | 20 | 20 | 20 | 24 |
| 20 | 23 | 23 | 23 | 26 |
| 21 | 24 | 28 | 24 | 29 |
| 22 | 28 | 32 | 28 | 32 |
| 28 | 42 | 53 | 42 | 53 |
| 36 | 63 | 89 | 63 | 89 |
| 48 | 96 | 160 | 96 | 160 |
| 59 | 128 | 244 | 128 | 244 |
| 60 | 131 | 252 | **133** | 252 |

The upper bounds are far from the lower bounds for large `n`; closing that gap is a genuine research problem.

---

## 3. Progress on the “at least four points on a line” variant

### 3.1 A general lower bound from the four-axis construction
Xiang Li’s September 2026 construction (OEIS A172992 attachment) places `n` rational points on four concurrent supporting lines and guarantees

\[
\left\lfloor\frac{n^2+12}{28}\right\rfloor
\]

lines containing **exactly** four selected points, for `n ≥ 28`.

In the **at-least-four** variant, the four supporting lines themselves also count. For `n ≥ 28`, each supporting line contains at least 4 selected points. Therefore:

\[
\boxed{T_4(n) \ge \left\lfloor\frac{n^2+12}{28}\right\rfloor + 4 \qquad (n\ge28).}
\]

All coordinates are rational; a uniform dilation makes them integers.

### 3.2 Exact verification at n=60
Parameters:

- `A = {0,…,12}` on `y=0`: `(1/(180+x), 0)`
- `B = {0,…,12}` on `x=0`: `(0, 1/(60+y))`
- `C = {4,…,20}` on `y=x`: `(1/(240+z), 1/(240+z))`
- `D = {-8,…,8}` on `y=-x`: `(1/(120+w), -1/(120+w))`

A transversal through the `A`-point `x` and `B`-point `y` also contains the `C`-point `z=x+y` and the `D`-point `w=x-y`. The number of pairs `(x,y)` with `x+y∈[4,20]` and `x-y∈[-8,8]` is **129**.

Exact line-size histogram:

| line size | count |
|---:|---:|
| 2 | 400 |
| 3 | 56 |
| 4 | 129 |
| 13 | 2 |
| 17 | 2 |

Thus

\[
T_4(60) \ge 129 + 4 = 133.
\]

This beats the best no-five board value at `n=60` (integer record 131, algebraic 129). It is a concrete new lower bound for the at-least-four variant, assuming no better published `T_4(60)` is found.

Certificate files:
- `atleast4_N60.txt`
- `atleast4_N60_certificate.json`
- `atleast4_li_results.json` (four-axis lower bounds for all `n=11…60`)

### 3.3 What this does and does not say
- This is **not** a contest entry: the four supporting lines have 13 or 17 points, so they violate the contest’s no-≥5 rule.
- It **is** a valid integer-coordinate orchard for the classical “at least four” variant.
- For `n=59`, the same construction gives `124+4=128`, tying the current no-five integer record.
- For `n<59`, the existing no-five records are already better than this four-axis lower bound.

---

## 4. Updated publishability assessment

The paper case is now stronger than “six leaderboard records”:

1. **Six new no-five lower bounds** (28, 36, 37, 40, 48, 49), with exact algebraic coordinates and honesty certificates.
2. **At-least-four variant:** general lower bound `T_4(n) ≥ floor((n²+12)/28)+4` for `n≥28`, and exact `T_4(60) ≥ 133`.
3. **Upper-bound package:** Melchior bound for `T_4`, exact slack identity, and `T_4(12) ≤ 8` via pseudoline non-realizability.
4. **Negative structural results:** integer records are deletion/growth-optimal under the tested operations; honesty gap explains why algebraic scores can fall below integer scores.

This is still a modest, computational/constructive note, but it now has a clearer mathematical package. Suitable framing:

> “Four-point orchards: algebraic constructions, at-least-four bounds, and Melchior-type upper bounds.”

Before arXiv:
- Confirm the six no-five scores and `T_4(60)=133` against unpublished contest archives / experts.
- Check whether the `T_4` Melchior corollary or `T_4(12)≤8` already appears in the line-arrangement literature.
- Include exact coordinates / certificates and the verifier code.
- Do not claim optimality except where proved (e.g. `n=11`, and the known small `n` values).

---

## 5. Key references
- MathWorld, *Orchard-Planting Problem*: Melchior upper bound, small values, `n=21` deletion bound, `n=22` Grünbaum–Rigby construction.
- OEIS A006065 / A172992; Xiang Li, *A Four-Point-Line Construction for OEIS A172992* (Sept 2026).
- Du, eMath blog: optimal 13–18, 20-tree problem, 22:28.
- Pegg, Wolfram Blog (2018): orchard-planting polynomials and 4-orchard examples.
- Elkies, arXiv:math/0612749: points-and-lines configurations, `t_k` / `T_k` notation.
- Solymosi–Stojaković, arXiv:1107.0327: asymptotic `n^{2-o(1)}` lower bounds for no-`(k+1)` `k`-tuples.
- Berman–Gévay–Richter-Gebert–Tabachnikov, arXiv:2408.09203: movable celestial `(n_4)` configurations.
- Bokowski–Pokora, arXiv:1607.05864: pseudoline arrangements, nonexistence of 12 pseudolines with 9 quadruple and 12 double points.
