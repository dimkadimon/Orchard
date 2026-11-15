# Four-point orchards

This repository contains code, certificates, search results and a revised draft of
[*Algebraic four-point orchards and a dual Melchior bound*](paper/orchard.pdf).
The revision needs the author's review and approval before submission.

## Results and scope

For an n-point set, distinguish exactly-four lines `t₄`, the exactly-four count
with **no five collinear** (the strict competition variant), and at-least-four
lines `T₄`. These are different optimisation problems. The revised paper centers on four exact ring constructions, giving
strict configurations with (n, score) = (28,42), (36,63), (40,70), (48,96).
A further nested algebraic construction gives strict scores **75, 78, 81, 84, 87**
for n = **42, 43, 44, 45, 46**, respectively. Each is one above the cited
integer-coordinate competition score. These comparisons do **not** establish
new global real-coordinate records; the irrational constructions cannot be
submitted to the integer-only contest.

For the at-least-four problem, an exact rational census establishes `T₄(60) ≥ 133`.
That configuration contains long lines and is **not** admissible in the strict
no-five variant. It is based on Xiang Li's four-axis construction; see the
research files for details. This older observation is not a result claimed in the
focused arXiv draft.

The earlier parity-only chord condition is incorrect. The earlier claims of
complete family enumeration, family-optimal scores, and a proven optimum for the
four-axis slope family are withdrawn; the relevant numerical searches are
heuristics or finite restricted computations. Novelty has not been fully
established and neither the proof nor code has undergone independent peer review.
The files in `docs/` are historical research notes, **not** the revised claims.

The focused [arXiv source package](arxiv_submission/orchard_arxiv_source.tar.gz)
and [submission notes](arxiv_submission/SUBMISSION_NOTES.md) include all four
ring figures, the five extensions, and their exact verification certificates.

## Search and verification

See [research/README.md](research/README.md) for methods, seeds, commands,
comparisons, recorded results and limitations. The new search includes exact
projective and stochastic beams, algebraic-field growth and annealing, finite-pool/grid MILPs,
and earlier deletion and extension searches. These searches are incomplete and
cannot prove that no improvements exist.

```bash
python -m pip install numpy sympy scipy
python verify.py entries/N11.txt entries/N12.txt entries/N13.txt entries/N28.txt entries/N36.txt entries/N40.txt entries/N48.txt
python verify_rational_censuses.py
python verify_determinants.py
python verify.py entries/extended_N{42,43,44,45,46}.txt
python verify_independent_algebraic.py entries/extended_N{42,43,44,45,46}.txt
```

The earlier determinant verifier uses a different census but shares number-field
arithmetic. The new independent algebraic verifier uses separate integer-polynomial
arithmetic and triple determinants; the rational checks use Python `Fraction`. Generate figures with `make_figures.py` (requires
matplotlib); compile the manuscript from `paper/` with a LaTeX installation:

```bash
cd paper && pdflatex orchard.tex && pdflatex orchard.tex
# or: tectonic -X compile orchard.tex
```

## Attribution

Dmitry Kamenetsky conceived the research idea, verified the manuscript, and is
responsible for errors. AI assistance contributed to literature investigation,
search and verification scripts, drafting and editing; see the paper attribution.

For the exploratory extension to five-point lines, see
[`research/five_point_investigation.md`](research/five_point_investigation.md).

MIT licence; see [LICENSE](LICENSE).
