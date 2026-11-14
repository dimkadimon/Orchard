> **Historical research note, not the current manuscript.** Some claims and record comparisons below are unverified or superseded, including statements of family optimality, completeness, novelty, and exhaustive global or local search. See [the revised paper](../paper/orchard.pdf) and [research notes](../research/README.md) for the corrected status.

# Integer records adopted from ISS

As of 2026-11-14, Johan Roos holds the ISS integer records from n=24 to 60 (up to 131 at n=60), which exceed the previous algebraic lower bounds at 23 values.

With the author's permission ("You can take the best integer coordinates from ISS site and use them to find new records!"), we have:

- **Adopted the ISS integer configuration wherever it exceeds the algebraic one.** `entries/N*.txt` now contains the best known configuration for every n (max of integer ISS and algebraic). Degree-1 entries are integer (poly `a = 0`), valid with no 5-line, score = integer record. They are lower bounds for `T_4(n)` and `t_4(n)`; `t_4^{(5)}` honest lower bounds remain the four celestial values at n=28,36,40,48.

- **Attempted to improve beyond ISS by seeding search from ISS integers.** Two searches were run from every integer record (24≤n≤60):
  - integer hill-climbing (random point move ±6, keep valid, accept improvements; plateau moves; occasional large jump), 40k steps per n;
  - grow/prune via 3-point line intersections (intersections of pairs of 3-point lines as candidate swaps).
  No improvement beyond the ISS score was found for any n (the run reproduced each ISS score). The search is not exhaustive but suggests the ISS table is locally optimal for these neighbourhoods.

## Result

After adoption, `entries/` matches the `integer_solutions.json` best scores where they were higher, and retains the degree-6/8 celestial entries where they beat integers:

| n | old alg | new (adopted) | ISS integer | status |
|---|---|---|---|---|
| 24 | 29 (d2) | 31 (d1) | 31 | tie |
| 26 | 33 (d6) | 34 (d1) | 34 | tie |
| 30 | 43 (d2) | 44 (d1) | 44 | tie |
| 31 | 45 (d2) | 46 (d1) | 46 | tie |
| 32 | 47 (d2) | 48 (d1) | 48 | tie |
| 33 | 49 (d2) | 51 (d1) | 51 | tie |
| 34 | 52 (d2) | 53 (d1) | 53 | tie |
| 39 | 65 (d2) | 67 (d1) | 67 | tie |
| 43 | 76 (d8) | 77 (d1) | 77 | tie |
| 44 | 78 (d2) | 80 (d1) | 80 | tie |
| 45 | 81 (d2) | 83 (d1) | 83 | tie |
| 46 | 84 (d2) | 86 (d1) | 86 | tie |
| 50 | 98 (d8) | 99 (d1) | 99 | tie |
| 51 | 100 (d2) | 102 (d1) | 102 | tie |
| 52 | 103 (d2) | 106 (d1) | 106 | tie |
| 53 | 107 (d2) | 109 (d1) | 109 | tie |
| 54 | 110 (d2) | 113 (d1) | 113 | tie |
| 55 | 114 (d2) | 116 (d1) | 116 | tie |
| 56 | 117 (d2) | 119 (d1) | 119 | tie |
| 57 | 120 (d2) | 122 (d1) | 122 | tie |
| 58 | 123 (d2) | 125 (d1) | 125 | tie |
| 59 | 126 (d2) | 128 (d1) | 128 | tie |
| 60 | 129 (d2) | 131 (d1) | 131 | tie |

The four celestial records still beat integers: 16,17,18,22,28,36,37,40,48,49 (e.g. n=28 42 vs 39, n=36 63 vs 59, n=40 70 vs 69, n=48 96 vs 93).

Previous algebraic entries that were replaced are archived in `entries_backup/`.

Scripts: `iss_seeded_search.py` (exact scorer via `Field([0,1])`, hill climb + intersection grow).

Verification: `python3 verify.py entries/N*.txt` shows `valid=True` for all; `honest=False` for degree-1 integers is expected (they are rational). `python3 -c "from algfield import parse_entry; from orchard import analyze; ..."` confirms scores.
