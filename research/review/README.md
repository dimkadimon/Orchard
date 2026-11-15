# Review follow-up: verification and restricted small-case search

Run commands from the repository root. The verification run used Python 3.13.18, SymPy 1.14.0 and NumPy 2.3.5; these are recorded in `verification.log`. Install SymPy and NumPy in your environment if needed. Parsers expect trusted repository inputs, not untrusted text.

```sh
python verify.py entries/N{11,13,14,15,28,36,40,48}.txt entries/extended_N{42,43,44,45,46}.txt
python verify_independent_algebraic.py entries/N{11,13,14,15,28,36,40,48}.txt entries/extended_N{42,43,44,45,46}.txt
python verify_rational_censuses.py
python second_order_k.py
python search_small_longline.py --width 5 --height 5 --result research/review/n14_grid_5x5.json
(cd paper && pdflatex -interaction=nonstopmode -halt-on-error orchard.tex && pdflatex -interaction=nonstopmode -halt-on-error orchard.tex)
```

The primary verifier uses its exact algebraic-field implementation and checks cross-ratios as well as the line census; the second census uses separate polynomial-remainder/determinant arithmetic but relies on SymPy for parsing, irreducibility, and real-root existence. It does **not** independently verify the cross-ratios. The rational census uses exact `Fraction` arithmetic. These tests are reproducible finite checks, not independent peer review or proofs of optimality.

## Restricted grid test

For the conditional five-point-line branch of a hypothetical 14-point, 11-rich-line configuration, the search enumerates each 9-point subset of a 5-by-5 integer grid. It retains subsets having exactly ten triple lines and no line with four or more of those nine points. For each retained set, it considers every line through at least two pairwise intersections of the ten triple lines that avoids the nine points. It checks whether at most five such intersection points on that line cover all ten triple lines. Arithmetic is exact in integer projective coordinates, including intersections at infinity.

The completed run examines all 2,042,975 subsets; eight meet the seed condition, 72 distinct seed/axis pairs are considered, and at most seven of the ten triple lines are covered by off-seed intersection points on any tested axis. There is no candidate in this restricted search. The script also supports adding generic points if a cover uses fewer than five and mapping any discovered candidate to an integer affine chart; independently census any future candidate before making claims.

The grid is a projectively restricted sample, not an exhaustive enumeration of all nine-point real sets. The no-four restriction and ten-triple seed condition omit other relevant five-point-line cases, including cases where an outside line contains four or more points. The negative result therefore does not settle T4(14). Likewise it says nothing conclusive about T4(15), real-coordinate optimality, or asymptotic bounds. The 9-line/10-triple classification at https://arxiv.org/abs/1306.6105 (Proposition 4.8) may be relevant to classifying the dual of this particular no-four outside-point subcase, but its applicability to a global proof has not been established.
