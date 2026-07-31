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
    # Monoclinic films appear in the shipped CSVs and in the UI filter, but the
    # generator that produced them was never committed and pipeline.geometry
    # has no monoclinic branch. Deriving one is a scientific decision, not a
    # port -- see docs/OPEN-QUESTIONS.md. Left empty so the build reports the
    # gap loudly instead of silently emitting wrong nets.
    "Monoclinic": (),
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


def fetch_stable_materials(api_key: str, max_elements: int = MAX_ELEMENTS) -> pd.DataFrame:
    """Pull thermodynamically stable materials and their conventional cells.

    Requires the ``mp-api`` package.
    """
    from mp_api.client import MPRester  # imported lazily; build-time only

    with MPRester(api_key) as mpr:
        docs = mpr.materials.summary.search(
            energy_above_hull=(0, 0.0),
            num_elements=(1, max_elements),
            fields=[
                "material_id", "formula_pretty", "elements",
                "symmetry", "structure", "nelements",
            ],
        )

    rows = []
    for doc in docs:
        lattice = doc.structure.lattice
        rows.append(
            {
                "material_id": str(doc.material_id),
                "formula": doc.formula_pretty,
                "elements": ", ".join(sorted(str(e) for e in doc.elements)),
                "num_elements": doc.nelements,
                "crystal_system": str(doc.symmetry.crystal_system).capitalize(),
                "point_group": doc.symmetry.point_group,
                "a": lattice.a,
                "b": lattice.b,
                "c": lattice.c,
            }
        )
    return pd.DataFrame(rows)


def build_film_nets(materials: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, list[dict]]:
    """Expand each material into its surface nets, one row per (material, plane)."""
    rows_2d: list[dict] = []
    rows_1d: list[dict] = []
    skipped: list[dict] = []

    for _, m in materials.iterrows():
        planes = PLANES_BY_SYSTEM.get(m["crystal_system"])
        if not planes:
            skipped.append(
                {"formula": m["formula"], "system": m["crystal_system"],
                 "reason": "no planes defined for this crystal system"}
            )
            continue

        for plane in planes:
            name = f"{m['formula']} ({plane})"
            try:
                net = new_plane(
                    name, m["crystal_system"], m["a"], m["b"], m["c"], plane
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
    args = parser.parse_args(argv)

    api_key = os.environ.get("MP_API_KEY")
    if not api_key:
        parser.error(
            "MP_API_KEY is not set. Get a key from https://next-gen.materialsproject.org/api "
            "and export it. The key previously hardcoded in app.py is in git history "
            "and must be treated as compromised."
        )

    materials = fetch_stable_materials(api_key, args.max_elements)
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
