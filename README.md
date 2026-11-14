# Orchard

New algebraic constructions for the four-point orchard-planting problem, and a
**second-order Melchior bound** that proves several small cases exactly by hand.

Companion code and certificates for the paper
[*New constructions of four-point orchards*](paper/orchard.pdf).

---

## The problem

Place `n` points in the real plane to maximise the number of lines through four of
them. Three variants are usually conflated; we keep them separate:

| notation | meaning |
|---|---|
| `t_4(n)` | max number of lines containing **exactly** 4 points |
| `t_4^(5)(n)` | the same, but **no 5 points collinear** (the competition rule) |
| `T_4(n)` | max number of lines containing **at least** 4 points |

So `t_4^(5)(n) ≤ t_4(n) ≤ T_4(n)`. Constructions with long rows are excellent for
`T_4` but inadmissible for `t_4^(5)`.

## Main results

### 1. A second-order Melchior bound (new)

Melchior's inequality is normally applied to the *dual* of the point set. We apply
it a **second time, to the arrangement of the rich lines themselves**. Writing
`L` for the number of lines carrying at least `k ≥ 4` points and `r_p` for the number
of rich lines through `p`:

```
sum_p w(r_p)  ≤  C(L,2) - 3,      w(0,1,2)=0,  w(3)=3,  w(j)=C(j,2)+j-3  (j ≥ 4)
```

This holds for **every** `k ≥ 4` and needs no "no-long-lines" hypothesis, so it bounds
`T_k` too. Combined with a majorisation (balanced-sequence) argument it becomes a
closed-form feasibility test — three inequalities in `n`, `k`, `L`.

Consequences for `k = 4`:

* **Exact values by hand:** `t_4(8)=2`, `t_4(9)=3`, `t_4(10)=5`, `t_4(11)=6`, `t_4(13)=9`.
  These previously required exhaustive computer search (Burr–Grünbaum–Sloane; Du).
* **Elementary `t_4(12) ≤ 8`**, replacing Bokowski–Pokora's classification of
  pseudoline arrangements (the dual Hesse configuration).
* Improvements at `n = 14, 15`.
* The bound is **sharp**: the optimal 13-point configuration attains equality.
* Also `T_4(n) = t_4(n)` at `n = 8, 9, 10, 11, 13`.

For larger `k` the improvements are bigger: the test beats the classical bound up to
`n = 25` for `k=5`, `n = 37` for `k=6`, `n = 51` for `k=7`, `n = 67` for `k=8`.

### 2. Four new algebraic lower bounds for `t_4^(5)`

| n | field degree | previous record | new | gain | family |
|---:|---:|---:|---:|---:|---|
| 28 | 6 | 39 | **42** | +3 | k4-ratio, m=7 |
| 36 | 6 | 59 | **63** | +4 | celestial, m=9 |
| 40 | 8 | 69 | **70** | +1 | celestial, m=10 |
| 48 | 8 | 93 | **96** | +3 | celestial, m=12 |

These four counts are **proved optimal within their family**: the incidence conditions
are linear in the radii over a real cyclotomic field, which confines the useful radii to
a finite set of cosine ratios, and that family can be enumerated completely
(`ring_family.py`, paper Theorem 5.4). For even `m` at most four rings are possible at all
(Theorem 5.5). See [`docs/ring_structure.md`](docs/ring_structure.md).

All are four-ring *celestial* sets: four concentric regular `m`-gons with consecutive
rings offset by half a step, radii chosen in a real number field of degree 6 or 8 so
that many transversals meet all four rings. Every coordinate is exact; every count is
decided by a determinant in `Q(a)`, never by a floating-point tolerance. Each entry is
certified *honest* (not a projective image of a rational configuration) by exhibiting
4-point lines with irrational cross-ratio.

Diagrams: [`paper/figures/`](paper/figures).

### 3. The at-least-four variant

`T_4(n) ≥ floor((n²+12)/28) + 4` for `n ≥ 28` — Xiang Li's four-axis construction, plus
the four supporting lines, which count in this variant.

At `n = 60` an exact census gives **`T_4(60) ≥ 133`** (129 transversals + 4 axes), with
line-size histogram `t_2=400, t_3=56, t_4=129, t_13=2, t_17=2`. This beats the best
no-five value at `n=60` (integer record 131), but is *not* a competition entry.

We also show the constant `1/28` is **optimal within the four-axis family**: over all
231 pairs of distinct primitive axis directions with coordinates in `[-4,4]`, the
maximum asymptotic density is `1/28`, attained only by the original slopes and their
unimodular equivalents. Grid-shaped point sets (disk, square, hexagon, triangle) are
much worse, about `n²/67.7`.

## What is *not* claimed

* The general bound `t_4(n) ≤ floor((n(n-1)-6)/14)` is **classical** (Melchior +
  pair-counting), and is stated on MathWorld's orchard page.
* Four-ring celestial sets are a **known** class of configurations
  (Berman–Gévay–Richter-Gebert–Tabachnikov). Our contribution is the parameter choices
  and the record counts, not the shape.
* The four lower bounds are new **to the best of our knowledge** after checking OEIS
  A006065 / A172992, Friedman's tables, MathWorld, Du, Pegg, Elkies and the
  configuration literature. They are, however, proved optimal *within the four-ring cosine-ratio family*
  (paper, Theorem 5.4), and for even `m` no more than four rings are possible at all
  (Theorem 5.5).

## Layout

Scripts live at the repository root and use paths relative to it, so **run everything
from the root**.

```
paper/orchard.tex, orchard.pdf   the paper (LaTeX + compiled PDF)
paper/figures/                   construction diagrams (PDF + PNG)
entries/N11.txt … N60.txt        50 exact algebraic submissions
docs/                            literature audit, novelty assessment, bound tables
docs/ring_structure.md           the (2,2) mechanism, optimality + ring-cap theorems

algfield.py                      exact arithmetic in Q[a]/p(a)
orchard.py                       exact scoring, cross-ratios, honesty check
second_order_k.py                the new bound, general k  ← start here
second_order_full.py             k=4 table
verify_second_order.py           brute-force check of the majorisation lemma
validate_theorem.py              inequality vs 100 known configurations
validate_k.py                    inequality vs 158 grid/pencil/random tests, k=4..8
ring_family.py                   COMPLETE enumeration of the celestial family  ← new
ring_optimal.py                  ranked + mixed-modulus ring search
ring_collision.py                chord normal forms and (2,2) pair tables
multiring_search.py              continuous hill-climb (negative result, kept)
make_figures.py                  regenerate the four construction figures
fig_t4b.py                       the T_4(60) parameter-space figure
atleast4_li.py                   four-axis optimiser, n=11..60
slope_opt.py, slope_opt2.py      axis-slope family optimisation
gridrich.py                      rich-line counts for grid shapes
verify.py                        re-verify submitted entries
```

## Reproduce

Requires Python 3 with `numpy`, `sympy`, `matplotlib` (and `pypdf` only to inspect the PDF).

```bash
pip install numpy sympy matplotlib

python3 second_order_k.py        # the new bound for k=4..8, with all improvements
python3 verify_second_order.py   # exhaustive check: balanced == optimal (n=11..18)
python3 validate_theorem.py      # 0 violations over 100 configurations
python3 validate_k.py            # 0 violations over 158 tests, k=4..8
python3 make_figures.py          # four figures, scores asserted exactly
python3 fig_t4b.py               # T_4(60) = 133 certificate + figure
python3 atleast4_li.py           # four-axis table, n=11..60
```

To rebuild the paper (needs a LaTeX installation; [Tectonic](https://tectonic-typesetting.github.io/) works):

```bash
cd paper && tectonic -X compile orchard.tex
```

## Novelty status

The 13-point equality case and the general-`k` improvements were checked against
100 record configurations and 158 synthetic arrangements with no violations. The
second-order argument itself has not been located in the literature; if you know a
prior reference, please open an issue — we would like to cite it.

## Citation

```bibtex
@misc{orchard2026,
  author = {Kamenetsky, Dmitry},
  title  = {New constructions of four-point orchards},
  year   = {2026},
  url    = {https://github.com/dimkadimon/Orchard}
}
```

## AI assistance

This work was carried out by Dmitry Kamenetsky with substantial assistance from an AI
system (Arena.ai, Agent Mode). The AI assistant surveyed the literature and record
tables, wrote and ran the search and exact-arithmetic verification code, produced the
figures, assisted in developing and numerically checking the arguments, and drafted the
manuscript. The author posed the problem, directed the research, selected the results,
checked the statements and proofs, and takes full responsibility for the content.

Every claim here ships with a machine-checkable certificate — exact coordinates in a
number field, an exact collinearity census, and a cross-ratio honesty test — so results
can be verified independently of how they were found. Failed searches are reported
alongside successful ones (`docs/ring_structure.md`). See §2.5 of the paper.

## Licence

MIT — see [LICENSE](LICENSE).
