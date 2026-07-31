"""Regression tests for the surface-net construction.

The Sapphire cases pin the bug described in geometry.py: the notebook bound
the parameter named ``c`` to ``a``, so every hexagonal non-basal plane was
computed as though c == a.
"""

from __future__ import annotations

import math

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


@pytest.mark.parametrize(
    "plane, expected_a, expected_b, shipped_b",
    [
        # plane      expected a       expected b                          what the old CSV shipped
        ("M", SAPPHIRE_A, SAPPHIRE_C, 4.763),
        ("A", SAPPHIRE_A, math.sqrt(3 * SAPPHIRE_A**2 + SAPPHIRE_C**2), 9.526),
        ("R", math.sqrt(3) * SAPPHIRE_A, SAPPHIRE_C, 8.249757996450562),
    ],
)
def test_sapphire_non_basal_planes_use_c(plane, expected_a, expected_b, shipped_b):
    """A, R and M must depend on c. Previously all three collapsed to c == a."""
    net = new_hexagonal_plane(f"Sapphire ({plane})", SAPPHIRE_A, SAPPHIRE_C, plane)
    assert isinstance(net, Rectangle)
    # Rectangle normalises so that a <= b.
    lo, hi = sorted((expected_a, expected_b))
    assert net.a == pytest.approx(lo)
    assert net.b == pytest.approx(hi)
    # And confirm we actually moved off the old, wrong value.
    assert net.b != pytest.approx(shipped_b)


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
