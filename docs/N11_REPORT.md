> **Historical research note, not the current manuscript.** Some claims and record comparisons below are unverified or superseded, including statements of family optimality, completeness, novelty, and exhaustive global or local search. See [the revised paper](../paper/orchard.pdf) and [research notes](../research/README.md) for the corrected status.

# N = 11 — optimal orchards for both boards

## Result

| board | previous best | new best | status |
|:--|--:|--:|:--|
| Algebraic | 5 | **6** | new record, and provably optimal |
| Integer / whole-number | 6 | **6** | already optimal; no score improvement is possible |

The same upper bound applies to both boards: **11 trees can make at most 6 lines of exactly four trees** (with no line of five or more).

---

## New algebraic N=11 orchard (6 lines)

```text
a^2 - 2 = 0
(0, 0), (2, 0), (-a + 1, 1), (-a + 2, -a), (1, -a - 1), (-a, a), (1, -1), (a, -a), (1, a - 1), (a + 2, a), (a + 1, 1)
```

Verified exactly in `Q[a]/(a^2 - 2)`:

- 11 distinct points
- 6 lines of exactly 4 points
- 0 lines of 5 or more points
- 6 certificate lines with irrational cross-ratio
- both real roots give the same line structure

A convenient way to see the construction is as two pencils of three lines, through `(0,0)` and `(2,0)`:

Through `(0,0)`:

```text
x + (a - 1)y = 0
x + y = 0
x - (a + 1)y = 0
```

Through `(2,0)`:

```text
x + (a + 1)y - 2 = 0
x - y - 2 = 0
x + (1 - a)y - 2 = 0
```

The 9 intersections of one line from the first pencil with one line from the second pencil, together with the two pencil centers, give the 11 trees. Each of the 6 pencil lines contains exactly 4 trees.

---

## Why 6 is the absolute maximum for N=11

Let `L` be the number of 4-tree lines, and let `r_p` be the number of such lines through tree `p`.

1. Through any one tree, each 4-tree line uses 3 other trees. Two different lines through that tree cannot share another tree. Since there are only 10 other trees,

   ```text
   r_p <= floor(10 / 3) = 3.
   ```

2. Counting tree-line incidences gives

   ```text
   sum_p r_p = 4L.
   ```

3. Counting pairs of 4-tree lines that meet at a tree gives

   ```text
   sum_p C(r_p, 2) <= C(L, 2),
   ```

   because two distinct straight lines meet in at most one point.

Now suppose `L = 7`. Then `sum_p r_p = 28` across 11 trees, with every `r_p <= 3`. The smallest possible value of `sum_p C(r_p, 2)` under those constraints is obtained by six trees on 3 lines and five trees on 2 lines:

```text
6*C(3,2) + 5*C(2,2) = 6*3 + 5*1 = 23.
```

But

```text
C(7,2) = 21.
```

So `L = 7` is impossible. Larger `L` is impossible as well. Therefore

```text
L <= 6.
```

The construction above reaches 6, so it is optimal. The existing integer record of 6 is therefore also optimal.

---

## Integer N=11 record (6 lines, optimal)

The current whole-number record is already optimal:

```text
(0, 0), (0, 6), (0, 9), (0, 10), (2, 8), (6, 4), (6, 6), (10, 0), (12, 2), (15, 0), (18, 0)
```

It has 6 lines of exactly four trees and no line of five or more.

---

## Files updated

- `entries/N11.txt` — new optimal algebraic entry
- `SUBMIT_batch1.txt` — regenerated with the new N=11 entry
- `results.json`, `final_rows.json` — updated score table
- `gallery.html` — updated visualization
