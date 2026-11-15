Reproducibility of the algebraic constructions
==============================================

Files extended_N42.txt through extended_N46.txt give ordered affine
coordinates in Q(a), where a^8-8*a^6+19*a^4-12*a^2+1 = 0.
All the incidences hold at every real root; the diagram uses the unique
root in (1.97,1.98). Each file's first line specifies the polynomial.

Requirement: Python 3 and SymPy (tested with SymPy 1.14). Run:

  python verify_independent_algebraic.py extended_N*.txt

Expected histograms (entries give the number of lines of each size):
  N=42: {2:297, 3:38, 4:75}
  N=43: {2:324, 3:37, 4:78}
  N=44: {2:349, 3:37, 4:81}
  N=45: {2:384, 3:34, 4:84}
  N=46: {2:420, 3:31, 4:87}

The verifier uses exact integer-polynomial determinant computations modulo
the defining polynomial and reports no_five=True for all five entries. It
uses SymPy for parsing, irreducibility and existence of a real root. Its
parser is intended for trusted inputs, not untrusted text. It checks the
construction, not global optimality or priority. The source repository is
https://github.com/dimkadimon/Orchard/.
