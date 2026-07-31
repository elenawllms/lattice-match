"""Build the film catalogue from the Materials Project.

    export MP_API_KEY=...
    python -m pipeline.films --out src/assets/data

This restores what commit ffb6432 ("removed unnecessary API call") deleted. The
original had the key hardcoded in app.py; it is read from the environment here
and never reaches the shipped data.

Films go through exactly the same plane -> surface-net construction as
substrates (:mod:`pipeline.geometry`), so the hexagonal c-vs-a correction
applies to the film catalogue too. The shipped stable_films_1d.csv lists
GaN (0001) with a = 5.192 A, which is GaN's c; its a is 3.19 A. Regenerating
through the corrected geometry is expected to fix that -- see
:func:`validate_known_films`, which asserts it.

NOTE: this module has not been executed, because doing so requires an API key
that the repo does not have (and the leaked one must be revoked). Treat the
first real run as needing review: check the build report, then
validate_known_films.
"""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path

import pandas as pd

from .geometry import UnsupportedPlane, new_plane
from .nets import Rectangle, RectangleRoot2, Square, Triangle

#: Planes offered per crystal system, matching the checklists in the app UI.
PLANES_BY_SYSTEM: dict[str, tuple[str, ...]] = {
    "Cubic": ("100", "110", "111"),
    "Tetragonal": ("001", "100", "101", "110"),
    "Orthorhombic": ("001", "010", "100", "011", "101", "110"),
    "Hexagonal": ("0001", "1-100", "1-102", "11-20"),
    "Trigonal": ("0001", "1-100", "11-20"),
    # Only the two monoclinic faces whose in-plane vectors are orthogonal in
    # the standard setting. (010) and the {110} faces are oblique nets; the
    # app offers (1-10)/(110) checkboxes that will simply match nothing.
    # See new_monoclinic_plane and docs/OPEN-QUESTIONS.md item 6.
    "Monoclinic": ("001", "100"),
}

#: Tolerance, in degrees, for treating a lattice angle as a right angle.
RIGHT_ANGLE_TOL = 0.5

#: Columns a cached Materials Project fetch must carry to be reusable.
REQUIRED_MATERIAL_COLUMNS = {
    "formula", "elements", "num_elements", "crystal_system", "point_group",
    "centring", "a", "b", "c", "alpha", "beta", "gamma",
}

#: Films larger than this along any axis are excluded, matching the note in the
#: app's film modal ("lattice parameters exceeding 16A are excluded").
MAX_PARAM = 16.0

#: Matching the app's note: "films with more than 3 elements are excluded".
MAX_ELEMENTS = 3

OUTPUT_COLUMNS_2D = [
    "name", "a", "b", "elements", "crystal_system",
    "point_group", "formula", "plane", "num_elements",
]
OUTPUT_COLUMNS_1D = [
    "name", "a", "elements", "crystal_system",
    "point_group", "formula", "plane", "num_elements",
]


def fetch_stable_materials(
    api_key: str,
    max_elements: int = MAX_ELEMENTS,
    cache: Path | None = None,
) -> pd.DataFrame:
    """Pull thermodynamically stable materials and their CONVENTIONAL cells.

    Requires the ``mp-api`` package.

    The ``structure`` field of a summary document is the PRIMITIVE cell. Using
    it directly would be wrong: primitive fcc silicon has a = 3.849 A, while
    the value this tool needs -- and the one the original catalogue used -- is
    the conventional 5.444 A. Every cubic entry would have been off by a factor
    of sqrt(2). Each structure is therefore reduced to its conventional
    standard form, which costs ~7 ms per material.
    """
    if cache is not None and cache.exists():
        cached = pd.read_parquet(cache)
        missing = REQUIRED_MATERIAL_COLUMNS - set(cached.columns)
        if missing:
            # A cache written before a schema change would otherwise be used
            # silently, and the defaults are not safe: a missing `centring`
            # would treat every R-centred material as primitive, quietly
            # producing the wrong hexagonal A- and R-plane meshes.
            print(f"cache {cache} is missing {sorted(missing)}; re-fetching")
        else:
            return cached

    from mp_api.client import MPRester  # imported lazily; build-time only
    from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

    with MPRester(api_key) as mpr:
        docs = mpr.materials.summary.search(
            energy_above_hull=(0, 0.0),
            num_elements=(1, max_elements),
            fields=[
                "material_id", "formula_pretty", "elements",
                "symmetry", "structure", "nelements",
            ],
        )

    rows: list[dict] = []
    for doc in docs:
        try:
            lattice = (
                SpacegroupAnalyzer(doc.structure)
                .get_conventional_standard_structure()
                .lattice
            )
        except Exception:  # noqa: BLE001 - symmetry analysis can fail on odd cells
            continue
        rows.append(
            {
                "material_id": str(doc.material_id),
                "formula": doc.formula_pretty,
                "elements": ", ".join(sorted(str(e) for e in doc.elements)),
                "num_elements": doc.nelements,
                "crystal_system": str(doc.symmetry.crystal_system).capitalize(),
                "point_group": doc.symmetry.point_group,
                # Bravais centring, needed for the hexagonal A/R meshes. The
                # Hermann-Mauguin symbol's first letter is the lattice type,
                # so R-3c and R3c give "R" and P6_3mc gives "P".
                "centring": "R" if str(doc.symmetry.symbol).startswith("R") else "P",
                "a": lattice.a,
                "b": lattice.b,
                "c": lattice.c,
                "alpha": lattice.alpha,
                "beta": lattice.beta,
                "gamma": lattice.gamma,
            }
        )

    df = pd.DataFrame(rows)
    if cache is not None:
        cache.parent.mkdir(parents=True, exist_ok=True)
        df.to_parquet(cache)
    return df


def build_film_nets(materials: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, list[dict]]:
    """Expand each material into its surface nets, one row per (material, plane)."""
    rows_2d: list[dict] = []
    rows_1d: list[dict] = []
    skipped: list[dict] = []

    has_angles = {"alpha", "beta", "gamma"}.issubset(materials.columns)

    for _, m in materials.iterrows():
        planes = PLANES_BY_SYSTEM.get(m["crystal_system"])
        if not planes:
            skipped.append(
                {"formula": m["formula"], "system": m["crystal_system"],
                 "reason": "no planes defined for this crystal system"}
            )
            continue

        # new_monoclinic_plane assumes the standard setting, alpha = gamma = 90.
        # A handful of Materials Project conventional cells do not satisfy that;
        # their (001)/(100) faces would be oblique, so skip rather than guess.
        if m["crystal_system"] == "Monoclinic" and has_angles:
            if (
                abs(m["alpha"] - 90.0) > RIGHT_ANGLE_TOL
                or abs(m["gamma"] - 90.0) > RIGHT_ANGLE_TOL
            ):
                skipped.append(
                    {"formula": m["formula"], "system": "Monoclinic",
                     "reason": f"not in standard setting "
                               f"(alpha={m['alpha']:.2f}, gamma={m['gamma']:.2f})"}
                )
                continue

        for plane in planes:
            name = f"{m['formula']} ({plane})"
            try:
                net = new_plane(
                    name, m["crystal_system"], m["a"], m["b"], m["c"], plane,
                    m.get("centring", "P"),
                )
            except UnsupportedPlane as exc:
                skipped.append({"formula": m["formula"], "plane": plane, "reason": str(exc)})
                continue

            common = {
                "name": name,
                "elements": m["elements"],
                "crystal_system": m["crystal_system"],
                "point_group": m["point_group"],
                "formula": m["formula"],
                "plane": f"({plane})",
                "num_elements": m["num_elements"],
            }

            if isinstance(net, Triangle):
                if net.a <= MAX_PARAM:
                    rows_1d.append({**common, "a": net.a})
                # A triangular net also matches rectangular substrates via its
                # centred rectangular cell, mirroring Triangle.coords_2d.
                rect_b = net.a * 3**0.5
                if max(net.a, rect_b) <= MAX_PARAM:
                    rows_2d.append({**common, "a": net.a, "b": rect_b})
                continue

            if isinstance(net, Square):
                a = b = net.a
            elif isinstance(net, RectangleRoot2):
                a, b = net.a, net.a * 2**0.5
            elif isinstance(net, Rectangle):
                a, b = net.a, net.b
            else:  # pragma: no cover - defensive
                skipped.append({"formula": m["formula"], "plane": plane,
                                "reason": f"unhandled net type {type(net).__name__}"})
                continue

            if max(a, b) <= MAX_PARAM:
                rows_2d.append({**common, "a": a, "b": b})

    return (
        pd.DataFrame(rows_2d, columns=OUTPUT_COLUMNS_2D),
        pd.DataFrame(rows_1d, columns=OUTPUT_COLUMNS_1D),
        skipped,
    )


def validate_known_films(films_1d: pd.DataFrame, films_2d: pd.DataFrame) -> list[str]:
    """Sanity-check the catalogue against literature values.

    The shipped data had GaN (0001) at 5.192 A -- GaN's c, not its a. These
    checks exist so that regression is caught on the first run rather than
    shipped again.
    """
    expected = {"GaN": 3.19, "ZnO": 3.25, "AlN": 3.11}
    problems = []
    for formula, a_expected in expected.items():
        hit = films_1d[films_1d["formula"] == formula]
        if hit.empty:
            problems.append(f"{formula}: no basal-plane entry found")
            continue
        a = hit["a"].iloc[0]
        if abs(a - a_expected) / a_expected > 0.05:
            problems.append(
                f"{formula} (0001): a = {a:.3f} A, expected ~{a_expected} A "
                f"-- looks like c is being used in place of a"
            )
    return problems


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--max-elements", type=int, default=MAX_ELEMENTS)
    parser.add_argument(
        "--cache",
        type=Path,
        default=None,
        help="parquet file to cache the Materials Project fetch in, so that "
             "re-runs do not re-download ~24k materials",
    )
    args = parser.parse_args(argv)

    api_key = os.environ.get("MP_API_KEY")
    if not api_key:
        parser.error(
            "MP_API_KEY is not set. Get a key from https://next-gen.materialsproject.org/api "
            "and export it. The key previously hardcoded in app.py is in git history "
            "and must be treated as compromised."
        )

    materials = fetch_stable_materials(api_key, args.max_elements, cache=args.cache)
    films_2d, films_1d, skipped = build_film_nets(materials)
    problems = validate_known_films(films_1d, films_2d)

    args.out.mkdir(parents=True, exist_ok=True)
    films_2d.to_csv(args.out / "stable_films_2d.csv", index=False)
    films_1d.to_csv(args.out / "stable_films_1d.csv", index=False)
    (args.out / "films_report.json").write_text(
        json.dumps(
            {
                "materials": len(materials),
                "rows_2d": len(films_2d),
                "rows_1d": len(films_1d),
                "skipped": skipped[:200],
                "skipped_total": len(skipped),
                "validation_problems": problems,
            },
            indent=2,
        )
    )

    print(f"materials      : {len(materials)}")
    print(f"rows (2d / 1d) : {len(films_2d)} / {len(films_1d)}")
    print(f"skipped        : {len(skipped)}")
    for p in problems:
        print(f"  VALIDATION: {p}")
    return 1 if problems else 0


if __name__ == "__main__":
    raise SystemExit(main())
