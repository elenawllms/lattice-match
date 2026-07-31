"""Tests for the coincidence-cell enumeration."""

from __future__ import annotations

import math

import pytest

from pipeline.constants import MAX_AXIS_RATIO, MCIA_MAX, MIN_PARAM
from pipeline.superlattices import (
    format_fraction,
    get_sublattices,
    get_superlattices,
)


def test_superlattice_area_never_exceeds_cap():
    for s in (1, 2, 3, 5):
        for spl in get_superlattices(4.0, 5.5, s):
            assert spl.mcia <= MCIA_MAX


def test_superlattice_axis_stacking_is_bounded():
    for spl in get_superlattices(3.2, 3.2, 1):
        assert spl.numer_a <= MAX_AXIS_RATIO
        assert spl.numer_b <= MAX_AXIS_RATIO


def test_superlattice_mcia_matches_its_definition():
    a, b, s = 4.1, 6.3, 2
    for spl in get_superlattices(a, b, s):
        assert spl.mcia == pytest.approx(spl.numer_a * a * spl.numer_b * b * s)


def test_superlattice_count_shrinks_monotonically_with_s():
    """Square.coords_2d breaks out of its rotation loop on the first empty
    result. That is only valid if the count is monotonically non-increasing."""
    counts = [len(get_superlattices(4.0, 4.0, s)) for s in (1, 2, 5, 10, 13, 17, 26)]
    assert counts == sorted(counts, reverse=True)


def test_sublattices_respect_minimum_parameter():
    a, b, s = 5.43, 7.68, 1
    for spl in get_superlattices(a, b, s):
        for sbl in get_sublattices(a, b, s, spl, is_square=False):
            assert sbl.a >= MIN_PARAM - 1e-9
            assert sbl.b >= MIN_PARAM - 1e-9


def test_reducible_fractions_are_rejected():
    """2a/2 must not be emitted; it duplicates a/1."""
    a = 5.0
    spl = next(s for s in get_superlattices(a, a, 1) if s.numer_a == 2)
    labels = [s.label for s in get_sublattices(a, a, 1, spl, is_square=True)]
    assert not any("2a/2" in l.replace(" ", "") for l in labels)


def test_square_nets_do_not_emit_axis_swapped_duplicates():
    a = 4.0
    spl = get_superlattices(a, a, 1)[0]
    square = get_sublattices(a, a, 1, spl, is_square=True)
    rect = get_sublattices(a, a, 1, spl, is_square=False)
    assert len(rect) == 2 * len(square)


def test_sublattice_dimensions_follow_the_recipe():
    a, b, s = 5.0, 6.0, 1
    spl = next(s2 for s2 in get_superlattices(a, b, s) if s2.numer_a == 2)
    for sbl in get_sublattices(a, b, s, spl, is_square=False):
        # every emitted axis is (integer * const * base) / integer
        ratio = sbl.a / (math.sqrt(s) * a)
        assert ratio == pytest.approx(round(ratio * 60) / 60, rel=1e-9)


@pytest.mark.parametrize(
    "numer, denom, s, axis, expected",
    [
        # When numer == 1 and s == 1 the axis letter carries the numerator
        # slot, giving "a/2". Otherwise the axis is appended after the
        # denominator, giving "3/2 b". Both forms appear in the shipped CSV.
        (1, 1, 1, "a", "a"),
        (2, 1, 1, "a", "2a"),
        (1, 2, 1, "a", "a/2"),
        (3, 2, 1, "b", "3/2b"),
        (1, 1, 5, "a", "√5a"),
        (2, 3, 5, "a", "2√5/3a"),
    ],
)
def test_format_fraction(numer, denom, s, axis, expected):
    assert format_fraction(numer, denom, s, axis).replace(" ", "") == expected
