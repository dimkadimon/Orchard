Reproducibility of the algebraic constructions
==============================================

Files N28.txt, N36.txt, N40.txt and N48.txt contain the four ring
configurations. Files extended_N42.txt through extended_N46.txt contain
the 40-point configuration's nested extensions. Each file gives ordered
affine coordinates in Q(a); the defining polynomial is on its first line.
Incidences hold at every real root of the corresponding irreducible
polynomial. The 46-point diagram uses its root in (1.97,1.98).

Requirements: Python 3, SymPy and NumPy. Run both exact implementations:

  python verify_independent_algebraic.py N*.txt extended_N*.txt
  python verify.py N*.txt extended_N*.txt

Expected histograms (numbers of lines of each size):
  N=28: {2:126, 4:42}
  N=36: {2:252, 4:63}
  N=40: {2:240, 3:40, 4:70}
  N=42: {2:297, 3:38, 4:75}
  N=43: {2:324, 3:37, 4:78}
  N=44: {2:349, 3:37, 4:81}
  N=45: {2:384, 3:34, 4:84}
  N=46: {2:420, 3:31, 4:87}
  N=48: {2:516, 3:12, 4:96}

The independent verifier uses exact integer-polynomial determinant computations modulo
the defining polynomial and reports no_five=True for all nine entries. It
uses SymPy for parsing, irreducibility and existence of a real root. Its
parser is intended for trusted inputs, not untrusted text. It checks the
constructions, not global optimality or priority. The primary verifier groups normalized field line equations and checks irrational
cross-ratios. The source repository is
https://github.com/dimkadimon/Orchard/.
