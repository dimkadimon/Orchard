# Research log: revised manuscript and record search

The source repository snapshot used for the revision was `a247d4ed227c54e18c36b1b048efec141ca363c2`; this work was applied on top of main commit `d30352fe7e960a398ca40779d3fdb191a3b74758`. Baseline comparison files `records_int_now.json`, `final_rows.json`, and `integer_solutions.json` are the pre-existing, pinned repository data; they are not a live or exhaustive literature survey. External context: [ISS rules and leaderboard](https://infinitesearchspace.com/contests/orchard), [OEIS A006065](https://oeis.org/A006065), and [OEIS A172992](https://oeis.org/A172992). The integer contest forbids five collinear; A172992's exactly-four definition permits longer lines. Exceeding an integer score is not proof of a new global real record.

## Methods

Run from the repository root, using Python 3 with numpy, sympy and scipy. Logs and JSON files in this directory record the measured runs; per-round logs lacking a final `RESULT` and JSON are interrupted/incomplete, not completed experiments. Times vary with machine load. Candidates, if produced, are in `research/candidates/` with coordinates and search metadata; For rational candidates use `python verify_rational_candidate.py research/candidates/FILENAME.txt` for an independent Fraction-based census; for algebraic candidates use `python verify.py research/candidates/FILENAME.txt`. The verifier CLI marks rational points as not algebraically honest, even if their strict score is valid. Independently audit any candidate before claiming a record.

* `beam_orchard.py`: exact integer homogeneous/projective, strict no-five delete-and-regrow beam search. It proposes 3/3, sampled 3/2 and sampled 2/2 line intersections, including points at infinity, and limits coordinate size, width and rounds. Example: `python beam_orchard.py --n 24 --depth 3 --width 8 --rounds 20 --pairs32 500 --pairs22 400 --seed 81024 --result research/beam/n24_d3.json`. The attempted sweeps use n=24–34 with seeds `81000+n` and 20 rounds; n=35–46 with `82000+n` and 16 rounds; and n=47–60 with `83000+n` and 12 rounds. All use depth 3, width 8 and 500/400 samples. Actual finished/partial rounds are in the logs; an attempt at a sweep is not a claim that all sizes finished. An additional n=60 run uses depth 4, width 4, 4 rounds, 300/200 samples, seed 85660.
* `beam_algebraic.py`: from the exact algebraic ring configurations n=28,36,40,48, delete 2 points then regrow a width-2 beam. Floating-point arithmetic only proposes intersections; completed states are checked by exact number-field census. Each has 3 rounds, samples 120/100, tolerance 1e-8 and seed `94000+n`. Example: `python beam_algebraic.py --n 28 --depth 2 --width 2 --rounds 3 --pairs32 120 --pairs22 100 --seed 94028 --result research/algebraic/n28.json`. The tolerance and candidate pruning make this a heuristic. The same algebraic arithmetic is used to generate and verify proposals, so independent review of any result remains important.
* `milp_orchard.py`: selects exactly n points from a finite candidate pool built from an integer seed and 28 projective intersection proposals. The model enumerates all lines of the pool exactly, forbids five selected points on any line, and maximizes strict four-lines. For n=24,30,39,49,60, run `python milp_orchard.py --n N --extra 28 --delete 4 --trials 5 --seed $((91000+N)) --seconds 60 --result research/milp/nN.json` (substitute the number for N). All five returned solver status optimal **within their particular pools** and matched their source scores, respectively 31,44,67,96,131. That is not a global upper bound.
* `milp_grid.py`: a distinct finite-pool construction search selecting n lattice points from a w×h integer grid, with the same strict model. Runs use `(n,w,h)=(24,7,7),(28,8,8),(32,9,9),(36,10,10)`, 75 seconds per instance. Example: `python milp_grid.py --n 24 --width 7 --height 7 --seconds 75 --result research/grid/n24_7x7.json`. Check JSON status: a time-limited incumbent is not a pool optimum; if no JSON exists the run did not finish.
* `search_deletions.py 6`: exactly enumerated all 21,686,018 deletions of one through six points from the n=28,36,40,48 certified sets, scoring from their exact incidence lists. No deletion improved the larger pinned algebraic/integer baseline. Exhaustive only for these source configurations and deletion sizes.
* `search_five_axes.py`: 16 seeded annealing restarts, n=60, for one five-axis direction family and interval parameters; its best at-least-four score was 128, versus the existing 133 from four axes. An independent rational census gave `{2:482,3:52,4:121,5:3,12:2,14:1,18:1}`. This is not a strict no-five set. It is not exhaustive.
* `search_extend_T4.py`: exact greedy extensions of the n=60 four-axis set to n=61–72. At 61 and 62 they tied the existing Li lower bound and subsequently fell below it. Greedy growth cannot rule out other extensions.

For logs marked optimal, numerical HiGHS termination and the exact line enumeration/selected incumbent are evidence for the finite model, not a formally verified solver proof. No solver optimality statement applies to all rational or real point sets. Search space, random choices, and time/coordinate caps limit all other runs.

## Verification and corrections

`python verify.py entries/N11.txt entries/N12.txt entries/N13.txt entries/N28.txt entries/N36.txt entries/N40.txt entries/N48.txt` verifies exact number-field line counts and no-five. `python verify_determinants.py` performs a different determinant-based census for the four large certificates, but shares the field parser and arithmetic. `python verify_rational_censuses.py` independently checks the rational n=12 and n=60 configurations using `Fraction`; the latter has 129 exactly-four lines and four longer axes, hence 133 at-least-four lines.

`slope_opt2.py` fixes the prior lost-offset and empty-initial-state errors, but remains a numerical heuristic. `ring_family.py` uses rounded floats and tolerance-based incidence grouping and cannot establish family completeness or optimality. See the revised manuscript for the counterexample and the conditional even-modulus ring-cap statement. The old `docs/` notes and pre-existing results files are retained as historical material and may contain superseded claims. The seven-page PDF was built from the revised source in the earlier revision workspace and copied here unchanged. The author has not yet approved this AI-assisted revision; independent mathematical review and a fuller novelty audit remain necessary.

## Recorded outcome of this bounded attempt

No candidate exceeded the pinned comparison scores. In the exact projective beam,
11 runs completed: n=24–29, 35–37, 47, and 60 (the latter used depth 4). These
runs evaluated 124,220 child states from 2,094,654 proposed points. Runs at
n=30, 38 and 48 were deliberately interrupted after 3, 1 and 2 completed rounds;
their partial logs are included, and remaining planned sizes were not run. The
n=60 beam regenerated a strict 131-line set, rather than beating it. The algebraic
beam completed all four seeded configurations, with 45 exact final-state checks
each, reproducing the source best scores. These counts include repeated or
seed-reproducing configurations, not necessarily distinct configurations.

The grid MILP found scores 25 (7×7, n=24; solver-optimal within that grid), 31
(8×8, n=28), 36 (9×9, n=32), and 40 (10×10, n=36). The latter three reached the
75-second limit and their optimality gaps remained open. Their selected point sets
in the JSON were separately verified by `verify_rational_candidate.py`: all had
no five collinear and the recorded counts. They did not exceed the pinned baselines.
All model optimality claims are subject to the numerical solver qualification above.
