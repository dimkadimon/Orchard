> **Historical research note, not the current manuscript.** Some claims and record comparisons below are unverified or superseded, including statements of family optimality, completeness, novelty, and exhaustive global or local search. See [the revised paper](../paper/orchard.pdf) and [research notes](../research/README.md) for the corrected status.

# Orchard Problem — Improvement Search, Novelty Audit, and arXiv Publishability

**Date:** 2026-11-14
**Scope:** N = 11…60, both the *integer* board (`/contests/orchard`) and the *algebraic* board (`/contests/orchard-algebraic`).

---

## 1. Bottom line

- **Improvement search:** I ran a thorough, exact improvement pass over **both** boards starting from the current records. **No score at any N = 11…60 improved on either board.** The records are deletion-optimal, growth-optimal, and (for the integer board) have no low-degree point that can be "algebraicized." They represent the current state of the art.
- **Novelty:** The genuinely new mathematical results are the **algebraic constructions already in our table**. At **six sizes — N = 28, 36, 37, 40, 48, 49 — they beat every published / integer-coordinate lower bound**, giving new lower bounds for the *no-five-collinear* 4-orchard problem. Everywhere else the values either match known optima (N ≤ 20) or fall below the integer records (a consequence of the algebraic board's *honesty* rule).
- **arXiv:** A short, carefully-framed note **is plausibly publishable**, but it is a *modest computational / constructive contribution*, not a breakthrough. It needs (a) a confirmed novelty check against the configuration literature, and (b) an honest framing. Details and a recommended outline are in §6.

---

## 2. What was tried (all exact, all reproducible)

| # | Method | Board | Result |
|---|--------|-------|--------|
| 1 | Deletion closure, depth ≤ 8, branch-and-bound on every integer record | integer | **0 improvements** |
| 2 | Integer growth (insert integer points on 3-lines / integral 3-line intersections) | integer | **0 improvements** |
| 3 | Deletion closure, depth ≤ 6, over **164** honest algebraic sources (records + all seed pools), honesty preserved | algebraic | **0 improvements** |
| 4 | "Delete min-degree point + add algebraic point" (`int(N) − mindeg(N) + 1`) | algebraic | **0 improvements** (reproduces our existing entries) |
| 5 | Rational growth + delete-and-regrow (any rational point clears to an integer point) | both | **0 in-range improvements** |
| 6 | Algebraic growth / padding of honest seeds over Q(√2) | algebraic | **0 improvements** |

Scripts: `int_closure.py`, `alg_closure.py`, `rational_grow.py`, `alg_grow_pass.py`, `alg_two_point.py`. Results: `*_results.json` (all empty in-range).

### Why the boards resist improvement
- **Integer records are dense and deletion-optimal.** Average 4-line degree is ≈ 4L/N (≈ 8.7 at N = 60), while consecutive records differ by only ≈ 3 lines. Deleting any point loses far more than the next record gains, so deletion closure cannot help.
- **No "algebraicizable" point.** To turn an integer orchard into an *honest* algebraic one without losing lines, you need a point lying on ≤ 1 four-line (so it can slide to an algebraic position along that line). **Only N = 12 has such a point**; every other integer record has minimum 4-degree ≥ 2 (mostly ≥ 3). A degree-d point is the rigid intersection of d lines and cannot move.
- **The honesty constraint is a real gap.** The algebraic board refuses any orchard that is a projective/affine image of an integer one (certified by an irrational cross-ratio). So at sizes where the integer record is strong but no good honest algebraic construction exists, our algebraic entry is `int(N−1) + 1`, which is 2–3 lines *below* `int(N)`. This is a structural limitation, not a search failure.

---

## 3. Current state of both boards

- **Algebraic board:** reset and re-populated **today (2026-11-14)**; the live records now equal our delivered table at **49/50** sizes (holder "Dmitry Kamenetsky"). **Only N = 11 differs: the board has 5, our entry has 6** (see §5).
- **Integer board:** unchanged; records by Burfitt, Hinrichs, Ripatti, Mammel, Kamenetsky, Mucha, Gehrke, Roos. These already **exceed OEIS A172992 at most N ≥ 24** (the site is ahead of the OEIS entry).

Our monotone algebraic table and the integer records are in `final_rows.json` / `records_int_now.json`; the full comparison is `novelty_table.md`.

---

## 4. Novelty audit

Baselines checked: **OEIS A006065** (real, exact only to N = 20), **OEIS A172992** (integer, to N = 60), **Friedman** (k = 4, to N = 25), **Du** (optimal 13–18; symmetric 22 : 28, Dec 2019), **Pegg / MSE** (22–25), the **2011 Infinite Search Space contest** (= A172992), and the current site boards.

### 4.1 New lower bounds (the publishable core)

For the **no-five-collinear** 4-orchard problem, our honest algebraic constructions exceed the best integer-coordinate construction at exactly these sizes:

| N | our algebraic | best integer / published | margin | field | construction |
|---|---|---|---|---|---|
| 28 | **42** | 39 | +3 | deg 6, `a⁶−7a⁴+14a²−7` | k4-ratio (m = 7) |
| 36 | **63** | 59 | +4 | deg 6, `a⁶−6a⁴+9a²−3` | celestial-k4 (m = 9) |
| 37 | **63** | 62 | +1 | deg 6 | celestial-k4 + 1 |
| 40 | **70** | 69 | +1 | deg 8, `a⁸−8a⁶+19a⁴−12a²+1` | celestial-k4 (m = 10) |
| 48 | **96** | 93 | +3 | deg 8, `a⁸−8a⁶+20a⁴−16a²+1` | celestial-k4 (m = 12) |
| 49 | **97** | 96 | +1 | deg 8 | celestial-k4 + 1 |

All six entries **re-verified exactly**: correct score, no ≥ 5-line, and an irrational cross-ratio on ≥ 1 four-line (honesty certificate). See `entries/N{28,36,37,40,48,49}.txt`.

### 4.2 Not novel
- **N ≤ 20:** our values equal OEIS A006065, which are *proved optimal* (Du). Nothing new.
- **N = 11 = 6:** the value is A006065(11), known since 2010. Our simple optimality proof (in `N11_REPORT.md`) is clean but the *value* is not new. It **is** a new *board* record (board currently 5).
- **N = 22 = 28:** Du's symmetric construction (Dec 2019); the board marks it `era: prior`. Not new.
- **Sizes where algebraic < integer record** (24, 26, 30–34, 39, 43–46, 50–60): the integer record is the better lower bound; our algebraic entry adds nothing there.

### 4.3 Important caveat — two different problems
- The **contest** (both boards) forbids **any line with ≥ 5 trees**. This is the *no-five-collinear* variant.
- **OEIS A006065 / A172992** are the *classical* orchard problem, where lines with > 4 trees are allowed (they just don't count). Classical constructions such as the grid/octagon dual (~n²/24) and **Xiang Li's Sept-2026 bound `a(n) ≥ ⌊(n²+12)/28⌋`** use **long lines** (parallel families dualize to collinear points) and are **invalid for the contest**.
- Consequently our contest orchards are valid lower bounds for the classical problem too, but the *cleanest* novelty claim is for the **no-five-collinear variant**, where the integer contest records are the benchmark and we beat them at the six sizes above.

### 4.4 Novelty still to confirm
The **celestial-k4 / k4-ratio** constructions are related to Grünbaum's "celestial" and Grünbaum–Rigby configuration families. The *orchard-optimization* (choosing radii to maximize 4-lines) and the specific scores appear new, but this **must be confirmed** against:
- Grünbaum & Rigby, *The real configuration (21₄)* (1990); Grünbaum's configuration/celestial papers;
- the Al Zimmermann "Orchard planting" contest archive (recmath.org, problem 22) and `gpj.connectfree.co.uk`;
- recent arXiv: Lu **2610.14099** ((n₄) configurations), Solymosi–Stojaković **1107.0327**, Elkies **math/0612749**.

---

## 5. Action item: submit N = 11 = 6

The algebraic board still shows **N = 11 = 5**. Our `entries/N11.txt` (score 6, honest, both ±√2 roots identical) is a **new board record** and is already in `SUBMIT_batch1.txt`. It is *proved optimal* (`N11_REPORT.md`). Submit it to claim the record. (The integer board already has 6 = optimal.)

---

## 6. arXiv publishability assessment

**Can it go on arXiv?** Technically yes — arXiv is a moderated preprint server (not peer-reviewed); new authors need an endorsement. The real question is whether it is a *meaningful, defensible* contribution.

### What would make it publishable
1. **New lower bounds** for the no-five-collinear 4-orchard problem at N = 28, 36, 37, 40, 48, 49 (concrete, exactly verifiable).
2. **The "honest algebraic orchard" framework**: certifying non-rational realizability via an **irrational cross-ratio** on a 4-line (a projective invariant). This is a clean, reusable methodological idea.
3. **Number-field constructions** (celestial-k4, k4-ratio) with exact coordinates and a reproducible exact-arithmetic verifier.
4. A clear **no-5 vs classical** discussion (avoids overclaiming against long-line constructions).

### Honest limitations
- Margins are small (+1 to +4 lines) and only at 6 sizes.
- N ≤ 20 is already optimal; the integer board was not improved.
- It is a **computational / constructive** result, not a new theorem or asymptotic bound.
- Novelty of the constructions vs the configuration literature is **not yet fully established**.

### Verdict
A **short note (≈ 6–10 pages) in `math.CO` / `math.MG` is reasonable**, framed as:
> *"Honest algebraic orchards: number-field constructions and new lower bounds for the four-point-line orchard problem with no five collinear."*

It is better suited to a problem/experimental venue (**Geombinatorics**, **Journal of Integer Sequences**, **Contributions to Discrete Mathematics**, **Experimental Mathematics**) than to a top journal. An arXiv preprint is appropriate to establish priority.

### Recommended outline
1. Introduction: orchard problem, the no-5 variant, prior records (OEIS, Friedman, Du, contests).
2. Honest algebraic orchards: definition, cross-ratio certificate, why projective images of integer orchards are excluded.
3. Constructions: celestial-k4 / k4-ratio families; exact coordinates in number fields.
4. New lower bounds: the six-size table, exact verification, honesty certificates.
5. Computational methodology: exact field arithmetic, deletion closure, growth; the negative results (boards are deletion/growth-optimal).
6. N = 11 optimality proof (short, self-contained).
7. Open problems: closing the honesty gap at large N; asymptotics of the no-5 variant.

### Before submitting
- **Confirm novelty** against the configuration literature and contest archives (§4.4); consider contacting experts (Kamenetsky, Du, Grünbaum school).
- **Contribute the new lower bounds to OEIS** (A006065 comments / a no-5 variant sequence) to establish priority.
- Obtain an **arXiv endorsement** if needed; prepare exact-coordinate appendices + verifier code.
- Do **not** overclaim: present as new *lower bounds* and a *framework*, not as optimal values or a breakthrough.

---

## 7. Files
- `novelty_table.md` — full N = 11…60 comparison vs all baselines.
- `novelty_data.json` — machine-readable comparison + novel-candidate list.
- `entries/N{28,36,37,40,48,49}.txt` — the six new-lower-bound orchards (exact, honest).
- `N11_REPORT.md` — N = 11 optimality proof + construction.
- `int_closure.py`, `alg_closure.py`, `rational_grow.py`, `alg_grow_pass.py` — the (negative) improvement searches.
- `final_rows.json`, `records_int_now.json`, `records_alg_now.json` — current tables / live boards.
