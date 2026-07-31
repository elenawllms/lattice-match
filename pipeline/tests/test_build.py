"""Golden-diff the regenerated tables against the originally shipped CSVs.

Every difference must trace to a documented fix. This test is what makes it
safe to replace the production data: it proves the physics corrections changed
exactly what they were supposed to and nothing else.
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest

from pipeline.build import BuildReport, build_nets, build_tables
from pipeline.colors import assign_colors

REPO = Path(__file__).resolve().parents[2]
SHIPPED = REPO / "src" / "assets" / "data"
CRYSTEC = REPO / "pipeline" / "data" / "crystec.csv"

# Faces whose geometry is expected to move, because the hexagonal handler now
# receives a real c. Sapphire is the only included substrate with non-basal
# hexagonal faces.
EXPECTED_CHANGED_FACES = {"Sapphire (A)", "Sapphire (M)", "Sapphire (R)"}

# Faces that did not exist before, because there was no trigonal branch.
EXPECTED_NEW_FACES_2D = {
    "LiNbO3 (0001)", "LiNbO3 (1010)", "LiNbO3 (1120)",
    "LiTaO3 (0001)", "LiTaO3 (1010)", "LiTaO3 (1120)",
}
EXPECTED_NEW_FACES_1D = {"LiNbO3 (0001)", "LiTaO3 (0001)"}

KEY = ["substrate", "dimensions", "angle"]


@pytest.fixture(scope="module")
def regenerated():
    report = BuildReport()
    nets = build_nets(pd.read_csv(CRYSTEC), report)
    assign_colors(nets)
    sub_2d, sub_1d = build_tables(nets)
    return sub_2d, sub_1d, report


def _geometry_diff(old: pd.DataFrame, new: pd.DataFrame, cols: list[str]) -> set[str]:
    """Faces present in both tables whose lattice parameters moved."""
    shared_faces = set(old.substrate) & set(new.substrate)
    o = old[old.substrate.isin(shared_faces)].set_index(KEY)[cols].sort_index()
    n = new[new.substrate.isin(shared_faces)].set_index(KEY)[cols].sort_index()
    common = o.index.intersection(n.index)
    o = o.loc[common].groupby(level=[0, 1, 2]).first()
    n = n.loc[common].groupby(level=[0, 1, 2]).first()
    delta = (o - n).abs().max(axis=1)
    return {idx[0] for idx in delta[delta > 1e-9].index}


def test_2d_only_expected_faces_changed(regenerated):
    sub_2d, _, _ = regenerated
    old = pd.read_csv(SHIPPED / "sublattices_2d.csv")
    assert _geometry_diff(old, sub_2d, ["a", "b", "mcia"]) == EXPECTED_CHANGED_FACES


def test_1d_geometry_is_untouched(regenerated):
    _, sub_1d, _ = regenerated
    old = pd.read_csv(SHIPPED / "sublattices_1d.csv")
    assert _geometry_diff(old, sub_1d, ["a", "mcia"]) == set()


def test_only_trigonal_faces_were_added(regenerated):
    sub_2d, sub_1d, _ = regenerated
    old_2d = pd.read_csv(SHIPPED / "sublattices_2d.csv")
    old_1d = pd.read_csv(SHIPPED / "sublattices_1d.csv")
    assert set(sub_2d.substrate) - set(old_2d.substrate) == EXPECTED_NEW_FACES_2D
    assert set(sub_1d.substrate) - set(old_1d.substrate) == EXPECTED_NEW_FACES_1D


def test_no_faces_were_lost(regenerated):
    """The corrections must not drop anything the old pipeline produced."""
    sub_2d, sub_1d, _ = regenerated
    old_2d = pd.read_csv(SHIPPED / "sublattices_2d.csv")
    old_1d = pd.read_csv(SHIPPED / "sublattices_1d.csv")
    assert not set(old_2d.substrate) - set(sub_2d.substrate)
    assert not set(old_1d.substrate) - set(sub_1d.substrate)


def test_sapphire_m_plane_now_uses_c(regenerated):
    """The headline correction, asserted on the actual output table.

    Compares the 1:1 registry, labelled "(a × b)" -- not the max, since
    superlattices stack several cells along each axis.
    """
    sub_2d, _, _ = regenerated
    one_to_one = sub_2d[
        (sub_2d.substrate == "Sapphire (M)") & (sub_2d.dimensions == "(a × b)")
    ]
    assert len(one_to_one) == 1
    assert one_to_one["a"].iloc[0] == pytest.approx(4.763, rel=1e-9)
    assert one_to_one["b"].iloc[0] == pytest.approx(13.003, rel=1e-9)

    # What the shipped CSV had: c collapsed onto a.
    old = pd.read_csv(SHIPPED / "sublattices_2d.csv")
    old_one_to_one = old[
        (old.substrate == "Sapphire (M)") & (old.dimensions == "(a × b)")
    ]
    assert old_one_to_one["b"].iloc[0] == pytest.approx(4.763, rel=1e-9)


def test_every_exclusion_is_recorded_with_a_reason(regenerated):
    """Nothing may be dropped silently, which is how four substrates vanished."""
    _, _, report = regenerated
    assert {e["substrate"] for e in report.substrates_excluded} == {"CdS", "CdSe"}
    for entry in report.substrates_excluded + report.faces_skipped:
        assert entry["reason"] and entry["reason"] != "no reason recorded"


def test_skipped_faces_are_only_oblique_111_planes(regenerated):
    """(111) of a tetragonal/orthorhombic crystal is an oblique net."""
    _, _, report = regenerated
    assert {f["face"] for f in report.faces_skipped} == {
        "LiAlO2 (111)", "LiGaO2 (111)", "TiO2 (111)", "YAlO3 (111)"
    }


def test_output_schema_matches_the_app(regenerated):
    """app.py reads these column names directly."""
    sub_2d, sub_1d, _ = regenerated
    assert list(sub_2d.columns) == list(pd.read_csv(SHIPPED / "sublattices_2d.csv").columns)
    assert list(sub_1d.columns) == list(pd.read_csv(SHIPPED / "sublattices_1d.csv").columns)


def test_color_column_format_is_preserved(regenerated):
    sub_2d, _, _ = regenerated
    assert sub_2d["color_column"].iloc[0].startswith("rgb(")
    channels = sub_2d[["R", "G", "B"]]
    assert channels.ge(0).all().all() and channels.le(1).all().all()


# --- film catalogue --------------------------------------------------------


def test_film_nets_use_the_corrected_hexagonal_geometry():
    """The shipped catalogue listed GaN (0001) at 5.192 A, which is GaN's c."""
    from pipeline.films import build_film_nets

    materials = pd.DataFrame([
        dict(formula="GaN", elements="Ga, N", num_elements=2,
             crystal_system="Hexagonal", point_group="6mm",
             a=3.189, b=3.189, c=5.192, alpha=90.0, beta=90.0, gamma=120.0),
    ])
    _films_2d, films_1d, _skipped = build_film_nets(materials)
    assert films_1d["a"].iloc[0] == pytest.approx(3.189)

    shipped = pd.read_csv(SHIPPED / "stable_films_1d.csv")
    assert shipped[shipped.formula == "GaN"]["a"].iloc[0] == pytest.approx(5.192, abs=1e-3)


def test_film_validation_catches_the_c_for_a_regression():
    """validate_known_films must fail loudly if c is ever used as a again."""
    from pipeline.films import validate_known_films

    bad = pd.DataFrame([dict(formula="GaN", a=5.192)])
    problems = validate_known_films(bad, pd.DataFrame())
    assert any("GaN" in p and "c is being used" in p for p in problems)

    good = pd.DataFrame([
        dict(formula="GaN", a=3.19), dict(formula="ZnO", a=3.25), dict(formula="AlN", a=3.11),
    ])
    assert validate_known_films(good, pd.DataFrame()) == []


def test_monoclinic_materials_off_standard_setting_are_skipped():
    """new_monoclinic_plane assumes alpha = gamma = 90; anything else is oblique."""
    from pipeline.films import build_film_nets

    materials = pd.DataFrame([
        dict(formula="Bad", elements="X", num_elements=1, crystal_system="Monoclinic",
             point_group="2/m", a=5.0, b=9.0, c=7.0, alpha=97.0, beta=90.0, gamma=90.0),
    ])
    films_2d, _films_1d, skipped = build_film_nets(materials)
    assert films_2d.empty
    assert "standard setting" in skipped[0]["reason"]
