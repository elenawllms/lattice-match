"""Generation parameters for the superlattice search.

These bound the combinatorial enumeration in :mod:`pipeline.superlattices`.
Values are carried over unchanged from ``update_substrate_list.ipynb``.
"""

#: Largest coincident interface area to enumerate, in square angstroms.
#: Superlattices above this are considered too large to be physically useful.
MCIA_MAX = 200.0

#: Smallest permitted lattice parameter along either axis, in angstroms.
#: Subdividing a cell below this produces nets no real film would match.
MIN_PARAM = 3.0

#: Most unit cells that may be stacked along a single axis.
MAX_AXIS_RATIO = 5

#: Primes used to reject reducible numerator/denominator pairs, e.g. 2a/2,
#: which would duplicate a cell already enumerated as a/1. Must cover every
#: prime up to MAX_AXIS_RATIO.
FIRST_FEW_PRIMES = (2, 3, 5, 7)
