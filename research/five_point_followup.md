# Five-point follow-up: certificates, bounded searches, and long-line case

## Scope and result

For a real n-point set let `t5` count lines containing exactly five points, `T5` count lines containing at least five, and `s5` maximize this count subject to no six collinear. In the strict setting `t5=T5`; in general they differ.

**No new global record has been established.** This work does independently certify a published 24-point/18-line *strict* real construction, improving the prior repository's verified lower bound for n=24 from 15 to 18. It also constructs a rational 24-point/16-line strict configuration, an improvement on the repository's prior rational 15 but not a literature record: a rational 24/16 construction was already reported in [Du's post 24](https://bbs.emath.ac.cn/thread-50370-3-1.html). We have not proved that this particular arrangement is projectively inequivalent to prior ones.

## Algebraic 24/18 certificate

The projective coordinates and incidence pattern are from [Zhao Hui Du, post 23](https://bbs.emath.ac.cn/thread-50370-3-1.html). Set `a` to the real root of `a^3+a^2-1=0`. The polynomial is irreducible over Q and has a real root. `five_algebraic_24.py` transcribes all 24 projective triples, chooses the nonvanishing chart `z-5x-4y`, converts to affine coordinates modulo the polynomial, and writes `entries/five_Du24.txt`.

Two different exact implementations give histogram `{2: 75, 3: 7, 5: 18}`: normalized number-field line keys from `algfield.py`/`orchard.py`, and the independent integer-polynomial determinant census `verify_independent_algebraic.py`. Pair-accounting is `75+3*7+10*18=276=binom(24,2)`. Thus all points are distinct, exactly 18 lines contain five points, and no line contains six or more. Since the existing general upper bound is `T5(24)<=19`,

**18 <= s5(24) <= t5(24) <= T5(24) <= 19.**

The OEIS entry [A008997](https://oeis.org/A008997/internal) reports the exactly-five value 18, with an optimality claim; its classification/proof has not been independently audited here. Our upper bound does not prove that claim, and does not establish `T5(24)=18`.

The script exhaustively scores all subsets of this single 24-point configuration for n=18,...,23, with best strict counts 7,8,10,11,13,15 respectively. The outputs and independently checked coordinates are in `five_Du24_deletions.json` and `entries/five_Du24_subset_N*.txt`. This particular seed does **not** provide the sought 22/14 or 23/16 certificates. A deletion search on one seed is not a global bound.

## Rational extension

`five_rational_extension.py` adds one point to Du's published 23-point rational seed #0, testing intersections of its 4-point lines with all seed-spanned lines disjoint from them. It checks 29 distinct candidates, of which 27 have no six collinear. Its best 24-point result has histogram `{2:86, 3:10, 5:16}`; the added original-projective point is `(24:35:70)`. Both exact projective line normalization and an independent determinant census check `entries/five_rational_N24.txt`; the full projective point set is in `five_rational_N24.json`. This is a certificate, not a claim of global novelty or optimality.

## Exact beam and finite-grid MILP

`five_beam.py` uses primitive integer homogeneous coordinates. It deletes 3 points from the known rational 25/18 configuration, then regrows by intersections of spanned lines. It examines all 4/4, 4/3 and 3/3 pairs, samples 4/2, 3/2 and 2/2 pairs, screens no-six incidences exactly, and retains a bounded randomized beam. Two reproducible runs:

```
python five_beam.py --source 25 --delete 3 --rounds 20 --width 8 --branch 15 --samples 100 --seed 17 --baseline 18 --result research/five_beam25.json
python five_beam.py --source 25 --delete 3 --add 1 --rounds 20 --width 8 --branch 15 --samples 100 --seed 29 --baseline 21 --result research/five_beam26.json
```

At n=25: 69,058 point proposals, 5,100 exact child-state checks, best 18. At n=26: 100,209 proposals, 7,500 exact child-state checks, best 18; the target 21 reflects the published exactly-five 26/21 arrangement in [Elkies](https://arxiv.org/html/math/0612749v1), whose strictness we have not checked. The searches are incomplete because of the seed, finite rounds, sampling, width, branch and bit caps; neither failure is an upper bound. Counts are candidate occurrences, not necessarily distinct geometries.

`five_grid_milp.py` optimizes `T5` among 25-point subsets of the 40 points `{0,...,7} x {0,...,4}`. Every line supported by at least five grid points is represented; point variables and binary line variables enforce at least five selected points on each counted line. HiGHS reports a zero-gap optimum of 13 **for this finite pool only**; its selected set has histogram `{2:113,3:19,5:13}`, checked separately by exact determinant census. The solver's floating-point optimization report is not an independently checkable formal infeasibility certificate. See `five_grid_milp25.json`. This restricted construction is below the known 25/18 construction.

## Long-line upper-bound reduction

Let a configuration contain a line with `q>=6` points, and let `m=n-q`. Every other line containing at least five points intersects that line in at most one selected point, so contains at least four of the m outside points. Each such rich line therefore contains at least six distinct pairs of outside points, and no outside pair belongs to two lines. Consequently

`T5 <= 1 + floor(binomial(n-q,2)/6)`.

This covers all real or projective configurations including the fully collinear case. Applying it to the prior general upper-bound layers gives:

| n | candidate T5 layer | excludes a q-point line when q >= |
|---|---:|---:|
| 19 | 11 | 8 |
| 20 | 12 | 9 |
| 21 | 13 | 9 |
| 22 | 15 | 9 |
| 23 | 17 | 9 |
| 24 | 19 | 9 |
| 25 | 23 | 9 |

For example, any hypothetical 24-point configuration attaining 19 rich lines cannot have nine or more collinear points, but the 6-, 7-, and 8-point-line cases still require treatment. This reduction is rigorous but does **not** close the remaining one-unit gaps or audit Du's claimed complete classification.

## Reproduce

```
python five_algebraic_24.py
python five_rational_extension.py
OPENBLAS_NUM_THREADS=1 python five_grid_milp.py
python verify_independent_algebraic.py entries/five_Du24.txt entries/five_rational_N24.txt
python -m py_compile five_algebraic_24.py five_rational_extension.py five_beam.py five_grid_milp.py
```

The algebraic verification prints both `no_five` and `no_six` flags; earlier versions printed a single four-orchard-specific `strict` flag, which was potentially misleading for this five-point task.
