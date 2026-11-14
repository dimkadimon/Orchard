> **Historical research note, not the current manuscript.** Some claims and record comparisons below are unverified or superseded, including statements of family optimality, completeness, novelty, and exhaustive global or local search. See [the revised paper](../paper/orchard.pdf) and [research notes](../research/README.md) for the corrected status.

# The (2,2) mechanism, family optimality, and the ring cap

This note records the structural results added after the initial record search, and the
new approaches that were tried. Everything here is reproducible from the repository root.

## 1. The observation

Classify each 4-point line of a four-ring celestial set by the partition of 4 giving how
many of its points lie on each ring. Census of the four records:

| n | m | (2,2) | (1,1,1,1) | total |
|---:|---:|---:|---:|---:|
| 28 | 7 | 21 | 21 | 42 |
| 36 | 9 | 54 | 9 | 63 |
| 40 | 10 | 70 | 0 | 70 |
| 48 | 12 | 96 | 0 | 96 |

For `n=40` and `n=48` **every** 4-point line is of type `(2,2)`: it meets two rings twice
each. No 3-point lines occur at `n=36`; `n=40` has 40 and `n=48` has 12.

## 2. Why (2,2) lines are easy to design

A chord of span `a` of a ring of radius `r`, modulus `m`, offset `phi` lies on

```
x . (cos t, sin t) = d,    d = r*cos(pi*a/m),    t = phi + pi*(2j+a)/m
```

(the foot of the perpendicular from the centre is the chord midpoint). So two rings
share a `(2,2)` line iff

```
phi_i - phi_j in (pi/m)Z      and      r_i*cos(pi*a/m) = r_j*cos(pi*b/m)
```

The second condition is **linear in the radii** over the real cyclotomic field
`Q(cos(pi/m))`. Two consequences:

* Rotation by `2pi/m` fixes a ring, so relative offsets reduce to `{0, pi/m}` — two classes.
* The useful radius ratios form the explicit finite set
  `R_m = { cos(pi a/m) / cos(pi b/m) }`.

The radii of all four records are exactly of this shape. For `m=12` they are
`cos(pi a/12)/cos(pi/12)` for `a=1,2,3,4`.

This is why the earlier *continuous* search failed: hill-climbing radii with a soft
collinearity penalty scored 0–7 against records of 15–43 over 96 ring signatures. The
conditions are codimension-one, so generic perturbation destroys incidences. They must
be solved, not searched.

## 3. Completeness of the family

Define `F(m,k)`: `k` concentric regular `m`-gons, `r_0 = 1`, `r_i in R_m` distinct,
offsets in `{0, pi/m}`.

**Lemma.** If every pair of rings carries at least one `(2,2)` line, the set lies in
`F(m,k)` up to rotation, reflection and scaling. (Apply the criterion to the pair
`(R_0, R_i)`.)

The hypothesis is a real restriction, but it holds for all four records, and a pair with
no `(2,2)` line forfeits its whole contribution.

## 4. Theorem: optimality within F(m,4)

Complete enumeration (`ring_family.py`), radii sorted since rings are unordered:

| m | n | \|R_m\| | configs | max t_4 | record |
|---:|---:|---:|---:|---:|---:|
| 5 | 20 | 3 | 0 | — | 23 |
| 6 | 24 | 3 | 0 | — | 31 |
| 7 | 28 | 7 | 160 | **42** | 42 |
| 8 | 32 | 7 | 160 | 40 | 48 |
| 9 | 36 | 13 | 1760 | **63** | 63 |
| 10 | 40 | 13 | 1760 | **70** | 70 |
| 11 | 44 | 21 | 9120 | 77 | 80 |
| 12 | 48 | 19 | 6528 | **96** | 96 |
| 13 | 52 | 31 | 32480 | 91 | 106 |
| 14 | 56 | 31 | 32480 | 98 | 119 |

**The four records are exactly the maxima of their family.** For `m=5,6` the family is
empty: `R_m` has only 3 elements, so four distinct radii are impossible.

Census method: at each point group directions to all other points modulo `pi` with a
tolerance rule, counting a line only at its least-indexed point. Cross-checked against
an independent determinant-based `O(n^3)` census on all four records (42/63/70/96 agree).

> Implementation note: a hash of rounded line coordinates **undercounts** (it splits
> lines at rounding boundaries — it gave 38 instead of 42 for `n=28`). Angular grouping
> with a tolerance is the correct approach.

## 5. Theorem: ring cap for even modulus

**For even `m`, a union of concentric regular `m`-gons with offsets in `(pi/m)Z` and no
five collinear points has at most 4 rings.**

*Proof.* For even `m` the span `a=m/2` gives `cos(pi a/m)=0`: each ring contains `m/2`
diameters, with directions `phi_i + pi/2 + (2pi/m)Z`. Two rings share a diameter iff
their offsets agree mod `2pi/m`. Offsets in `(pi/m)Z` reduce to two classes `0, pi/m`,
whose diameter sets are disjoint. Three rings in one class would share a diameter
carrying `2+2+2 = 6` points. So at most 2 rings per class, hence at most 4. ∎

Verified by exhaustion: for `m=6,8,10,12,14`, **every** triple of distinct radii from
`R_m` in a single class yields a ≥5-point line (1, 35, 286, 969, 4495 triples, no
exception). For odd `m` there are no diameters and valid same-class triples do exist.

So the `m=10` and `m=12` records use the maximum permitted number of rings.

## 6. New approaches tried

| approach | coverage | outcome |
|---|---|---|
| Continuous hill-climbing, 3–6 rings | 96 signatures | **failed**, scores 0–7 vs records 15–43 |
| `F(m,3)` complete | `6≤m≤20` | no improvement (e.g. 70 at n=60 vs 131) |
| `F(m,5)` complete | `5≤m≤12` | **no admissible member** (all have a 5-line) |
| `F(m,6)` complete | `5≤m≤10` | **no admissible member** |
| Mixed moduli, 2–3 rings, `lcm≤60` | 10 multisets completed | no improvement; not exhaustive |

Mixed moduli is a genuine enlargement: rings of moduli `m, m'` can share chord directions
when offsets differ by an element of `pi Z / lcm(m,m')`, and the condition becomes
`r cos(pi a/m) = r' cos(pi b/m')`. Multisets completed include `{15,15,15}`,
`{12,15,15}`, `{10,15,15}`, `{9,15,15}`, `{12,12,15}`, `{10,12,15}`, `{20,20}`, `{19,19}`,
`{14,14,14}`, `{13,13,13}`. Cases with `lcm>60` were skipped as combinatorially out of
reach. **This is the most promising direction left open.**

## 7. Scripts

```
ring_collision.py    chord normal forms, (2,2) pair tables, first restricted search
ring_optimal.py      part A: ranked search; part B: mixed-modulus search
ring_family.py       COMPLETE enumeration of F(m,k) + robust census  <- main result
multiring_search.py  the failed continuous hill-climb (kept as a record)
```

Outputs: `ring_family_complete.json`, `ring_optimality.json`, `ring_collision_results.json`,
and logs `ring_family_all_log.txt`, `ring_mixed_log.txt`, `multiring_log.txt`.

## 8. Honest summary

No new record came out of these attempts. What came out instead is stronger than another
record: a proof that the four existing records are the best possible in their family, and
a theorem explaining why the family has exactly four rings. The search space that remains
is now explicitly delimited rather than vaguely "unexplored".
