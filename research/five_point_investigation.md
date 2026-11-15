# Five-point orchards: literature and feasibility investigation

Checked 14 November 2026. This is a preliminary investigation, not a claim that the literature has been exhaustively surveyed or that the results are new.

## Variants and existing results

Use the manuscript's definitions: `s5(n)` counts exactly-five-point lines under the no-six restriction, `t5(n)` counts exactly-five-point lines without that restriction, and `T5(n)` counts lines with at least five points. Thus `s5(n) ≤ t5(n) ≤ T5(n)`. The 2026 [OEIS A008997 entry](https://oeis.org/A008997/internal) reports exactly-five values for n ≤ 24 and says the listed terms have been proved optimal, attributing updates to Zhao Hui Du. Du's [five-point discussion](https://bbs.emath.ac.cn/thread-50370-1-1.html) distinguishes real, rational, and complex results, and explicitly corrects errors in earlier posts. We have not independently audited the complete classification underlying the reported n=19–24 exactly-five upper bounds; nor should these exactly-five reports be treated as at-least-five bounds.

Older [MathWorld](https://mathworld.wolfram.com/Orchard-PlantingProblem.html) and [Friedman](https://erich-friedman.github.io/packing/trees/) tables give broader intervals for many of these cases. In particular, the 2026 OEIS report is more specific; those broader tables do not establish the newer reports or their proofs. Elkies describes a 25-point, 18-line five-orchard in [On some points-and-lines problems and configurations](https://arxiv.org/html/math/0612749v1); [Pegg's 2018 account](https://blog.wolfram.com/2018/02/02/cultivating-new-solutions-for-the-orchard-planting-problem/) supplies barycentric coordinates.

## What the present method proves

The manuscript's Theorem 2.1 and Corollary 2.3 already apply to k=5. For n ≥ 6, the classical Melchior pair bound is

    T5(n) ≤ floor((n(n − 1) − 6)/24).

Combine it with the balanced-degree pair and weight tests from Corollary 2.3. These apply to **at-least-five** lines, and therefore also bound the other two variants. For example, for n=25 and L=24, the sum of rich-line degrees is at least 120; the balanced weight is 5w(4)+20w(5)=275, greater than C(24,2)−3=273. The L=24 case is ruled out, so T5(25) ≤ 23 instead of the classical 24. For n=15 and 18, the sharper combined bound follows already from the balanced *pair-count* condition; attributing those improvements solely to the Melchior weight would be incorrect.

The following exact values follow from the combined upper test and the independently verified strict rational projective constructions included here:

| n | 5–8 | 9 | 10 | 11 | 12 | 13 | 14 | 15 | 16 | 17 | 18 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| s5(n) = t5(n) = T5(n) | 1 | 2 | 2 | 2 | 3 | 3 | 4 | 6 | 6 | 7 | 9 |

For n=5–8, two five-rich lines would require at least nine distinct points, and a strict one-line construction gives the lower bound. For n=9–18, reproduce the upper values with `second_order_k.bound(n, 5)` and the lower constructions using the scripts below. The n=15 construction is the set `(-(i+j), -ij)` for 0 ≤ i < j ≤ 5, i.e. the pairwise intersections of the six lines `y = i x + i²`. Its exact histogram is `{2: 42, 3: 1, 5: 6}`. The n=18 rational example is credited to Du and has histogram `{2: 63, 5: 9}`. These are alternative verifications of small values, not claims of new literature records. Rational projective configurations can be moved to a rational affine chart and scaled to integer coordinates, preserving incidences.

## Larger small cases: reports versus our verified bounds

The OEIS column reports exactly-five values and the entry's statement of optimality; the actual proofs for n=19–24 were not independently audited. The last two columns are results independently checked in this investigation. The strict lower constructions imply the same lower bound for both t5 and T5. The last column is an upper bound on **T5**, not necessarily the best published bound.

| n | OEIS reported t5(n) | Strict lower s5(n) ≥ | Combined-test upper T5(n) ≤ |
|---:|---:|---:|---:|
| 19 | 10 | 10 | 11 |
| 20 | 11 | 11 | 12 |
| 21 | 12 | 12 | 13 |
| 22 | 14 | 13 | 15 |
| 23 | 16 | 15 | 17 |
| 24 | 18 | 15 | 19 |
| 25 | not listed | 18 | 23 |

The n=19 and 20 examples are from Du's corrected forum posts; the n=21–23 lower bounds come from deleting points from Du's 23-point rational example. The n=24 and 25 lower bounds are obtained from the Elkies/Pegg configuration; the full n=25 histogram is `{2: 84, 3: 12, 5: 18}` with no line of six or more points. These are known constructions, not newly discovered records. Exhaustive deletion only optimizes subsets of the chosen fixed seed, not all real configurations. In particular, this investigation does not establish that 18 is the best possible at n=25 or that the OEIS exactly-five values for n=22–24 have no-six witnesses.

## Asymptotics and possible next steps

For unrestricted t5 and T5, a 5-by-m integer grid yields at least

    m + 2 Σ_{d=1}^{floor((m−1)/4)} (m−4d) = m²/4 + O(m)

five-point transversals, since each passes through one point on each of the five horizontal rows. Hence `t5(n), T5(n) ≥ n²/100 − O(n)` after adding generic points if necessary. The rows contain m points, so for m≥6 this is **not** an `s5` construction. The classical upper coefficient for all three is 1/24. These are explicit bounds, not a claim that the 1/100 lower constant or 1/24 upper constant is best known. For the strict variant, Solymosi and Stojaković's [construction](https://arxiv.org/abs/1107.0327) gives `s5(n) ≥ n^(2−O(1/sqrt(log n)))` for sufficiently large n. Erdős's question about whether strict k-orchards are o(n²) remains a substantial open direction; see also the [Erdős problem database](https://www.erdosproblems.com/588).

Promising follow-up work is (1) audit the 2026 exact-five classification and obtain checkable real/no-six certificates for n=22–24; (2) search for stronger strict and unrestricted n=25 constructions using exact projective growth and bounded deletion/regrowth; and (3) seek a complete real upper-bound argument for the one-unit gaps at n=19–24, treating long lines separately. These are proposals, not results.

## Reproduction

Run from the repository root:

```sh
python investigate_five.py
python five_point_more_seeds.py
python five_point_examples.py
python verify_five_independent.py
python -c 'from second_order_k import bound,melchior_k; [print(n,melchior_k(n,5),bound(n,5)) for n in range(6,26)]'
```

`investigate_five.py` checks Du's 20-point strict rational projective seed and exhaustively considers every subset of sizes 9–20. `five_point_more_seeds.py` does the same for selected sizes of Du's 23-point seed. `five_point_examples.py` verifies the additional n=15, 18, 19, 24 and 25 examples. The generated JSON files record homogeneous coordinates and complete line-size histograms. Two independent census methods are used: primitive-integer projective line coefficients, and direct 3-by-3 determinant membership with lines deduplicated by point sets. Both use exact integer arithmetic and check the pair-count identity. They establish scores of the specific constructions, not global novelty or optimality. The `n=k` boundary in `second_order_k.py` was corrected: the classical formula is applicable only for n≥k+1, whereas n=k has one rich line.
