"""Regenerate the substrate superlattice tables from crystec.csv.

    python -m pipeline.build --out src/assets/data

Writes ``sublattices_2d.csv`` and ``sublattices_1d.csv``, plus a
``build_report.json`` recording every substrate face that could not be
represented. The notebook this replaces dropped such faces silently.
"""

from __future__ import annotations

import argparse
import json
from dataclasses import dataclass, field
from pathlib import Path

import pandas as pd

from .colors import assign_colors, rgb_string
from .geometry import UnsupportedPlane, new_plane
from .nets import SurfaceNet, Triangle

DEFAULT_CRYSTEC = Path(__file__).resolve().parent / "data" / "crystec.csv"


@dataclass
class BuildReport:
    """What made it into the output, and what did not."""

    substrates_in_source: int = 0
    substrates_excluded: list[dict] = field(default_factory=list)
    faces_built: int = 0
    faces_skipped: list[dict] = field(default_factory=list)
    rows_2d: int = 0
    rows_1d: int = 0

    def to_dict(self) -> dict:
        return {
            "substrates_in_source": self.substrates_in_source,
            "substrates_excluded": self.substrates_excluded,
            "faces_built": self.faces_built,
            "faces_skipped": self.faces_skipped,
            "rows_2d": self.rows_2d,
            "rows_1d": self.rows_1d,
        }


def build_nets(crystec: pd.DataFrame, report: BuildReport) -> list[SurfaceNet]:
    """Turn each included substrate/plane pair into a surface net."""
    report.substrates_in_source = len(crystec)

    if "include" in crystec.columns:
        excluded = crystec[~crystec["include"].astype(bool)]
        for _, row in excluded.iterrows():
            report.substrates_excluded.append(
                {"substrate": row["Name"], "reason": str(row.get("notes", "") or "no reason recorded")}
            )
        crystec = crystec[crystec["include"].astype(bool)]

    nets: list[SurfaceNet] = []
    for _, row in crystec.iterrows():
        for plane in str(row["Plane"]).split(", "):
            name = f"{row['Name']} ({plane})"
            try:
                nets.append(
                    new_plane(name, row["Structure"], row["a"], row["b"], row["c"], plane)
                )
            except UnsupportedPlane as exc:
                # Recorded, never silent. A (111) face of a tetragonal or
                # orthorhombic crystal is an oblique net, which none of the four
                # supported net types can represent.
                report.faces_skipped.append({"face": name, "reason": str(exc)})

    report.faces_built = len(nets)
    return nets


def build_tables(nets: list[SurfaceNet]) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Assemble the 2D and 1D sublattice tables from coloured nets."""
    frames_2d = [f for f in (n.coords_2d() for n in nets) if not f.empty]
    sub_2d = pd.concat(frames_2d).reset_index(drop=True)
    triangles = [n for n in nets if isinstance(n, Triangle)]
    sub_1d = pd.concat([t.coords_1d() for t in triangles]).reset_index(drop=True)

    for df in (sub_2d, sub_1d):
        df["color_column"] = [
            rgb_string(r, g, b) for r, g, b in zip(df["R"], df["G"], df["B"])
        ]
    return sub_2d, sub_1d


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--crystec", type=Path, default=DEFAULT_CRYSTEC)
    parser.add_argument("--out", type=Path, required=True, help="output directory")
    parser.add_argument(
        "--strict",
        action="store_true",
        help="exit non-zero if any substrate face could not be represented",
    )
    args = parser.parse_args(argv)

    crystec = pd.read_csv(args.crystec)
    report = BuildReport()

    nets = build_nets(crystec, report)
    assign_colors(nets)
    sub_2d, sub_1d = build_tables(nets)
    report.rows_2d, report.rows_1d = len(sub_2d), len(sub_1d)

    args.out.mkdir(parents=True, exist_ok=True)
    sub_2d.to_csv(args.out / "sublattices_2d.csv", index=False)
    sub_1d.to_csv(args.out / "sublattices_1d.csv", index=False)
    (args.out / "build_report.json").write_text(json.dumps(report.to_dict(), indent=2))

    print(f"substrates in source : {report.substrates_in_source}")
    print(f"faces built          : {report.faces_built}")
    print(f"rows (2d / 1d)       : {report.rows_2d} / {report.rows_1d}")
    for item in report.substrates_excluded:
        print(f"  EXCLUDED  {item['substrate']}: {item['reason']}")
    for item in report.faces_skipped:
        print(f"  SKIPPED   {item['face']}: {item['reason']}")

    if args.strict and (report.faces_skipped or report.substrates_excluded):
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
