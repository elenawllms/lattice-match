"""Coincidence-cell enumeration.

Given a substrate surface net, enumerate the film lattice parameters that can
be made commensurate with it, together with the minimum coincident interface
area (MCIA) each registry costs.

This is the standard coincidence-site-lattice / domain-matching-epitaxy
construction; see Zur & McGill, J. Appl. Phys. 55, 378 (1984).

Ported from ``update_substrate_list.ipynb`` cells 17-19. The numerics are
unchanged; only names, typing and documentation differ.
"""

from __future__ import annotations

from math import floor, sqrt
from typing import NamedTuple

from .constants import FIRST_FEW_PRIMES, MAX_AXIS_RATIO, MCIA_MAX, MIN_PARAM


class Superlattice(NamedTuple):
    """A coincidence cell spanning ``numer_a`` x ``numer_b`` substrate cells."""

    numer_a: int
    numer_b: int
    mcia: float


class Sublattice(NamedTuple):
    """A film net commensurate with a given superlattice.

    ``denom_a`` film cells span ``numer_a`` substrate cells along a, and
    likewise for b. ``label`` is the human-readable recipe, e.g. ``(2a x 3/2 b)``.
    """

    a: float
    b: float
    label: str


def format_fraction(numer: int, denom: int, s: int, axis: str) -> str:
    """Render one axis of a superlattice recipe, e.g. ``2a``, ``b/3``, ``2√5 a``.

    ``s`` is the rotational supercell multiplier: the cell is (√s x √s)R-theta.
    """
    d = "" if denom == 1 else ("/" + str(denom))
    if numer != 1:
        n = str(numer)
    elif s == 1:
        n = axis
    else:
        n = ""
    root = "" if s == 1 else ("√" + str(s))
    if axis not in n:
        if denom != 1:
            d += " " + axis
        else:
            d += axis
    return n + root + d


def get_superlattices(a: float, b: float, s: int) -> list[Superlattice]:
    """Enumerate coincidence cells over an ``a`` x ``b`` net in a (√s x √s) supercell.

    Returns every ``numer_a`` x ``numer_b`` stacking whose area stays under
    :data:`~pipeline.constants.MCIA_MAX`, capped by
    :data:`~pipeline.constants.MAX_AXIS_RATIO` along each axis.
    """
    area_ratio = MCIA_MAX / (a * b * s)
    const = sqrt(s)

    superlattices: list[Superlattice] = []
    a_limit = min(floor(area_ratio) + 1, floor(MAX_AXIS_RATIO / const) + 1)
    for numer_a in range(1, a_limit):
        b_limit = min(
            floor(area_ratio) + 1,
            floor(MAX_AXIS_RATIO / const) + 1,
            floor(area_ratio / numer_a) + 1,
        )
        for numer_b in range(1, b_limit):
            superlattices.append(
                Superlattice(numer_a, numer_b, numer_a * a * numer_b * b * s)
            )
    return superlattices


def get_sublattices(
    a: float,
    b: float,
    s: int,
    superlattice: Superlattice,
    is_square: bool,
) -> list[Sublattice]:
    """Subdivide one coincidence cell into the film nets that tile it.

    Rejects reducible fractions (2a/2 duplicating a/1) via prime divisors, and
    refuses to shrink either axis below :data:`~pipeline.constants.MIN_PARAM`.

    For non-square nets the axis-swapped variant is also emitted, representing
    the film rotated 90 degrees relative to the substrate.
    """
    const = sqrt(s)
    numer_a, numer_b, _mcia = superlattice
    sublattices: list[Sublattice] = []

    a_primes = [p for p in FIRST_FEW_PRIMES if numer_a % p == 0]
    b_primes = [p for p in FIRST_FEW_PRIMES if numer_b % p == 0]

    a_limit = min(floor(const * a * numer_a / MIN_PARAM) + 1, MAX_AXIS_RATIO + 1)
    for denom_a in range(1, a_limit):
        if any(denom_a % p == 0 for p in a_primes):
            continue

        b_limit = min(floor(const * b * numer_b / MIN_PARAM) + 1, MAX_AXIS_RATIO + 1)
        for denom_b in range(1, b_limit):
            if any(denom_b % p == 0 for p in b_primes):
                continue

            new_a = a * const * numer_a / denom_a
            new_b = b * const * numer_b / denom_b
            a_frac = format_fraction(numer_a, denom_a, s, "a")
            b_frac = format_fraction(numer_b, denom_b, s, "b")

            sublattices.append(
                Sublattice(new_a, new_b, f"({a_frac} × {b_frac})")
            )
            if not is_square:
                sublattices.append(
                    Sublattice(new_b, new_a, f"({b_frac} × {a_frac})")
                )

    return sublattices
