"""Regression tests for the surface-net construction.

The Sapphire cases pin the two bugs described in geometry.py: the notebook
bound the parameter named ``c`` to ``a``, so every hexagonal non-basal plane
was computed as though c == a; and the A- and R-plane expressions were
transposed onto each other's labels.

Values assume a primitive hexagonal lattice. Sapphire is actually R-3c, whose
centring changes A and R -- see docs/OPEN-QUESTIONS.md item 1.
"""

from __future__ import annotations

import math

import numpy as np
import pandas as pd
import pytest

from pipeline.geometry import (
    UnsupportedPlane,
    new_hexagonal_plane,
    new_plane,
)
from pipeline.nets import Rectangle, RectangleRoot2, Square, Triangle

# Sapphire, from crystec.csv
SAPPHIRE_A = 4.763
SAPPHIRE_C = 13.003


# The surface net of each plane, for a primitive hexagonal lattice. Derived by
# enumerating the lattice points lying in the plane through the origin and
# reducing to a shortest basis; all three come out rectangular at exactly 90
# degrees. See new_hexagonal_plane's docstring for the in-plane vectors.
SAPPHIRE_NETS = {
    "M": (SAPPHIRE_A, SAPPHIRE_C),                                        # a x c
    "A": (math.sqrt(3) * SAPPHIRE_A, SAPPHIRE_C),                         # sqrt(3)a x c
    "R": (SAPPHIRE_A, math.sqrt(3 * SAPPHIRE_A**2 + SAPPHIRE_C**2)),      # a x sqrt(3a^2+c^2)
}

# What the originally shipped CSV contained, when c was never passed and A/R
# were transposed. Every one of these must now be wrong.
SAPPHIRE_SHIPPED = {"M": (4.763, 4.763), "A": (4.763, 9.526), "R": (4.763, 8.249757996450562)}


@pytest.mark.parametrize("plane", ["M", "A", "R"])
def test_sapphire_non_basal_planes(plane):
    """A, R and M must depend on c, and A/R must be on the correct labels."""
    net = new_hexagonal_plane(f"Sapphire ({plane})", SAPPHIRE_A, SAPPHIRE_C, plane)
    assert isinstance(net, Rectangle)
    lo, hi = sorted(SAPPHIRE_NETS[plane])          # Rectangle normalises a <= b
    assert (net.a, net.b) == pytest.approx((lo, hi))
    assert (net.a, net.b) != pytest.approx(SAPPHIRE_SHIPPED[plane])


def test_a_and_r_are_not_transposed():
    """Regression guard for the specific bug: each formula was individually
    correct but attached to the other plane's label."""
    a, c = SAPPHIRE_A, SAPPHIRE_C
    a_net = new_hexagonal_plane("x", a, c, "A")
    r_net = new_hexagonal_plane("x", a, c, "R")

    # A-plane is spanned by |a1 - a2| = sqrt(3)a and by c.
    assert sorted((a_net.a, a_net.b)) == pytest.approx(sorted((math.sqrt(3) * a, c)))
    # R-plane is spanned by |a1 + a2| = a and by |-a1 + a2 + c|.
    assert sorted((r_net.a, r_net.b)) == pytest.approx(
        sorted((a, math.sqrt(3 * a**2 + c**2)))
    )
    # They must not be each other.
    assert (a_net.a, a_net.b) != pytest.approx((r_net.a, r_net.b))


def test_hexagonal_nets_are_rectangular_for_a_primitive_lattice():
    """The superlattice enumeration assumes orthogonal axes. Verify that
    assumption directly, by checking each net against the in-plane vectors."""
    a, c = 3.19, 5.19  # GaN, a genuine primitive hexagonal (P6_3mc) lattice
    a1 = np.array([a, 0.0, 0.0])
    a2 = np.array([-a / 2, a * math.sqrt(3) / 2, 0.0])
    cv = np.array([0.0, 0.0, c])

    for plane, (u, v) in {
        "M": (a2, cv),               # normal a1
        "A": (a1 - a2, cv),          # normal a1 + a2
        "R": (a1 + a2, -a1 + a2 + cv),
    }.items():
        assert abs(float(u @ v)) < 1e-9, f"{plane}: in-plane vectors are not orthogonal"
        net = new_hexagonal_plane("x", a, c, plane)
        expected = sorted((float(np.linalg.norm(u)), float(np.linalg.norm(v))))
        assert (net.a, net.b) == pytest.approx(expected), f"{plane} net mismatch"


def test_sapphire_basal_plane_unchanged():
    """C-plane never depended on c, so it must not move."""
    net = new_hexagonal_plane("Sapphire (C)", SAPPHIRE_A, SAPPHIRE_C, "C")
    assert isinstance(net, Triangle)
    assert net.a == pytest.approx(SAPPHIRE_A)


def test_triangle_recasts_as_centred_rectangular_cell():
    """A triangular net also enters the 2D table as the a x sqrt(3)a cell."""
    tri = Triangle("Sapphire (C)", SAPPHIRE_A)
    rect = Rectangle("Sapphire (C)", SAPPHIRE_A, math.sqrt(3) * SAPPHIRE_A)
    pd.testing.assert_frame_equal(tri.coords_2d(), rect.coords_2d())


def test_c_actually_changes_the_result():
    """Guard against the parameter being accepted but ignored."""
    m1 = new_hexagonal_plane("x", 4.763, 13.003, "M")
    m2 = new_hexagonal_plane("x", 4.763, 20.0, "M")
    assert m1.b != pytest.approx(m2.b)


def test_trigonal_is_supported():
    """LiNbO3 and LiTaO3 vanished entirely because there was no trigonal branch."""
    net = new_plane("LiNbO3 (0001)", "trigonal", 5.15, float("nan"), 13.86, "0001")
    assert isinstance(net, Triangle)
    assert net.a == pytest.approx(5.15)


@pytest.mark.parametrize(
    "plane, letter", [("0001", "C"), ("1010", "M"), ("1120", "A"), ("1-102", "R")]
)
def test_four_index_planes_alias_to_letters(plane, letter):
    """Crystec writes hexagonal planes as letters, trigonal as 4-index."""
    a, c = 5.15, 13.86
    assert type(new_hexagonal_plane("x", a, c, plane)) is type(
        new_hexagonal_plane("x", a, c, letter)
    )


def test_unknown_plane_raises_rather_than_returning_none():
    """The notebook returned None here and the caller skipped it silently."""
    with pytest.raises(UnsupportedPlane, match="not supported"):
        new_plane("CdS (111)", "hexagonal", 4.14, float("nan"), 6.76, "111")


def test_unknown_structure_raises():
    with pytest.raises(UnsupportedPlane, match="structure"):
        new_plane("X (100)", "triclinic", 5.0, 5.0, 5.0, "100")


@pytest.mark.parametrize(
    "plane, expected",
    [("100", Square), ("001", Square), ("111", Triangle), ("110", RectangleRoot2)],
)
def test_cubic_planes(plane, expected):
    assert isinstance(new_plane("Si", "cubic", 5.43, float("nan"), float("nan"), plane), expected)


def test_cubic_111_uses_root_two_convention():
    """(111) is stored as sqrt(2)*a; get_sublattices recovers a/sqrt(2) by dividing."""
    net = new_plane("Si (111)", "cubic", 5.43, float("nan"), float("nan"), "111")
    assert net.a == pytest.approx(math.sqrt(2) * 5.43)


def test_rectangle_orders_axes():
    assert Rectangle("x", 9.0, 4.0).a == pytest.approx(4.0)
    assert Rectangle("x", 9.0, 4.0).b == pytest.approx(9.0)


# --- monoclinic ------------------------------------------------------------
# Standard setting: alpha = gamma = 90, beta != 90, unique axis b. Verified
# against Materials Project conventional cells, where beta is the off-90 angle.


def test_monoclinic_001_is_the_a_by_b_rectangle():
    """(001) is spanned by a and b, and gamma = 90, so it is orthogonal."""
    net = new_plane("X (001)", "monoclinic", 5.0, 9.0, 7.0, "001")
    assert isinstance(net, Rectangle)
    assert (net.a, net.b) == pytest.approx((5.0, 9.0))


def test_monoclinic_100_is_the_b_by_c_rectangle():
    """(100) is spanned by b and c, and alpha = 90, so it is orthogonal."""
    net = new_plane("X (100)", "monoclinic", 5.0, 9.0, 7.0, "100")
    assert isinstance(net, Rectangle)
    assert (net.a, net.b) == pytest.approx((7.0, 9.0))


@pytest.mark.parametrize("plane", ["010", "110", "1-10", "011"])
def test_monoclinic_oblique_planes_are_rejected(plane):
    """(010) and the {110} faces are oblique; the enumeration assumes
    orthogonal axes, so they must be refused rather than approximated."""
    with pytest.raises(UnsupportedPlane, match="oblique"):
        new_plane(f"X ({plane})", "monoclinic", 5.0, 9.0, 7.0, plane)
