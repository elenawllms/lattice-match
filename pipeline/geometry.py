"""Map a bulk crystal plus a Miller plane onto its 2D surface net.

Ported from ``update_substrate_list.ipynb`` cells 3-7, with two corrections:

1. ``new_hexagonal_plane`` now receives BOTH ``a`` and ``c``. The notebook
   declared it as ``newHexagonalPlane(name, c, plane)`` but called it as
   ``newHexagonalPlane(name, a, plane)``, so the parameter named ``c`` was
   bound to ``a`` while the body separately read a module-global ``a``. Every
   hexagonal non-basal plane was computed as though c == a. For Sapphire
   (a=4.763, c=13.003) this shipped M-plane as 4.763 x 4.763 instead of
   4.763 x 13.003.

2. Trigonal structures are supported, and unknown structures or planes now
   raise instead of returning None. The notebook returned None and the caller
   skipped it, which silently dropped LiNbO3, LiTaO3 and several individual
   faces with no warning.
"""

from __future__ import annotations

import math

from .nets import Rectangle, RectangleRoot2, Square, SurfaceNet, Triangle


class UnsupportedPlane(ValueError):
    """Raised when a structure/plane pair has no defined surface net."""


# Hexagonal and trigonal planes are written either as Crystec's letter codes or
# as 4-index Miller-Bravais notation. Normalise both to the letter code.
HEXAGONAL_PLANE_ALIASES: dict[str, str] = {
    "C": "C", "0001": "C",
    "A": "A", "1120": "A", "11-20": "A",
    "M": "M", "1010": "M", "1-100": "M",
    "R": "R", "1012": "R", "1-102": "R",
}


def new_cubic_plane(name: str, a: float, plane: str) -> SurfaceNet:
    if plane in ("001", "100", "010"):
        return Square(name, a)
    if plane == "111":
        # The (111) face is triangular. Stored as √2·a rather than the true
        # nearest-neighbour spacing a/√2, because get_sublattices divides by
        # integer denominators and recovers a/√2 = (√2·a)/2.
        return Triangle(name, math.sqrt(2) * a)
    if plane in ("011", "110", "101"):
        return RectangleRoot2(name, a)
    raise UnsupportedPlane(f"cubic plane {plane!r} is not supported")


def new_orthorhombic_plane(
    name: str, a: float, b: float, c: float, plane: str
) -> SurfaceNet:
    nets = {
        "001": (a, b),
        "010": (a, c),
        "100": (b, c),
        "110": (math.hypot(a, b), c),
        "101": (math.hypot(a, c), b),
        "011": (math.hypot(b, c), a),
    }
    if plane not in nets:
        raise UnsupportedPlane(f"orthorhombic plane {plane!r} is not supported")
    return Rectangle(name, *nets[plane])


def new_tetragonal_plane(name: str, a: float, c: float, plane: str) -> SurfaceNet:
    if plane == "001":
        return Square(name, a)
    if plane in ("010", "100"):
        return Rectangle(name, a, c)
    if plane == "110":
        return Rectangle(name, math.sqrt(2) * a, c)
    if plane in ("101", "011"):
        return Rectangle(name, a, math.hypot(a, c))
    raise UnsupportedPlane(f"tetragonal plane {plane!r} is not supported")


def new_monoclinic_plane(name: str, a: float, b: float, c: float, plane: str) -> SurfaceNet:
    """Surface net for a monoclinic plane in the standard setting.

    Standard setting has alpha = gamma = 90 and beta != 90, with b the unique
    axis. Verified empirically against Materials Project conventional cells,
    where beta is the off-90 angle in ~93% of monoclinic entries and the
    remainder are within rounding of 90.

    Only the two faces whose in-plane vectors are mutually orthogonal can be
    represented:

      (001)  spanned by a and b, with gamma = 90  ->  Rectangle(a, b)
      (100)  spanned by b and c, with alpha = 90  ->  Rectangle(b, c)

    (010) is spanned by a and c at angle beta, and the {110} faces mix the
    a-b plane with c, which is inclined to a. Both are oblique nets, which the
    superlattice enumeration cannot express -- it assumes orthogonal axes.

    Callers must confirm alpha and gamma really are 90 for the specific
    material; see pipeline.films.
    """
    if plane == "001":
        return Rectangle(name, a, b)
    if plane == "100":
        return Rectangle(name, b, c)
    raise UnsupportedPlane(
        f"monoclinic plane {plane!r} is an oblique net (only 001 and 100 are "
        f"orthogonal in the standard setting)"
    )


def new_hexagonal_plane(name: str, a: float, c: float, plane: str) -> SurfaceNet:
    """Surface net for a hexagonal or trigonal (hexagonal-setting) plane.

    Both ``a`` and ``c`` are used. See the module docstring for the bug this
    replaces.

    The A- and R-plane formulas are carried over from the notebook unchanged
    apart from now receiving a real ``c``. They do not match the textbook
    derivation and are flagged for Falson Lab review in docs/OPEN-QUESTIONS.md;
    they are deliberately NOT altered here, because choosing a different
    convention is a scientific decision, not a bug fix.
    """
    key = HEXAGONAL_PLANE_ALIASES.get(plane)
    if key is None:
        raise UnsupportedPlane(
            f"hexagonal/trigonal plane {plane!r} is not supported "
            f"(expected one of {sorted(set(HEXAGONAL_PLANE_ALIASES))})"
        )
    if key == "C":
        return Triangle(name, a)
    if key == "A":
        return Rectangle(name, a, math.sqrt(3 * a**2 + c**2))
    if key == "R":
        return Rectangle(name, c, math.sqrt(3) * a)
    return Rectangle(name, a, c)  # M


def new_plane(
    name: str, structure: str, a: float, b: float, c: float, plane: str
) -> SurfaceNet:
    """Dispatch on crystal structure. Raises rather than returning None."""
    structure = structure.strip().lower()
    if structure == "cubic":
        return new_cubic_plane(name, a, plane)
    if structure == "orthorhombic":
        return new_orthorhombic_plane(name, a, b, c, plane)
    if structure == "tetragonal":
        return new_tetragonal_plane(name, a, c, plane)
    if structure == "monoclinic":
        return new_monoclinic_plane(name, a, b, c, plane)
    # Trigonal substrates here (LiNbO3, LiTaO3) are R3c in the hexagonal
    # setting, so their basal/prismatic faces use the hexagonal construction.
    if structure in ("hexagonal", "trigonal"):
        return new_hexagonal_plane(name, a, c, plane)
    raise UnsupportedPlane(f"structure {structure!r} is not supported")
