"""2D surface nets and the film lattice parameters they can match.

A bulk crystal cut along a Miller plane presents one of four net types.
Each knows which rotational supercells are commensurate with it, and can
enumerate the film nets it accepts.

Ported from ``update_substrate_list.ipynb`` cells 21-24.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field

import pandas as pd

from .superlattices import format_fraction, get_sublattices, get_superlattices
from .constants import FIRST_FEW_PRIMES, MAX_AXIS_RATIO, MCIA_MAX, MIN_PARAM

COORD_COLUMNS_2D = ["a", "b", "dimensions", "mcia", "angle"]
COORD_COLUMNS_1D = ["a", "dimensions", "mcia", "angle"]

# Commensurate rotations of a square net: s = m^2 + n^2 for the (m, n)
# direction, giving a (√s x √s)R-theta supercell. Ordered by increasing s,
# which the enumeration relies on to stop early (see Square.coords_2d).
SQUARE_ROTATIONS: tuple[tuple[int, tuple[int, int]], ...] = (
    (1, (0, 1)), (2, (1, 1)), (5, (1, 2)), (10, (1, 3)), (13, (2, 3)),
    (17, (1, 4)), (26, (1, 5)), (29, (2, 5)), (34, (3, 5)), (37, (1, 6)),
    (41, (4, 5)), (53, (2, 7)), (58, (3, 7)), (61, (5, 6)),
)

# Commensurate rotations of a triangular net: the classic √3 R30, √7 R19.1,
# √13 R13.9 supercells.
TRIANGLE_ROTATIONS: tuple[tuple[int, str], ...] = (
    (1, "0º"), (3, "30º"), (7, "19.1º"), (13, "13.9º"),
)


def _angle_label(m: int, n: int) -> str:
    return "{:.1f}º".format(math.degrees(math.atan2(m, n)))


@dataclass
class SurfaceNet:
    """Base class carrying the substrate label and its plot colour."""

    name: str
    R: float | None = field(default=None, init=False)
    G: float | None = field(default=None, init=False)
    B: float | None = field(default=None, init=False)

    def add_color(self, r: float, g: float, b: float) -> None:
        self.R, self.G, self.B = r, g, b

    def _finalize(self, rows: list, columns: list[str]) -> pd.DataFrame:
        df = pd.DataFrame(rows, columns=columns)
        df["substrate"] = self.name
        df["R"], df["G"], df["B"] = self.R, self.G, self.B
        return df

    def _concat(self, frames: list[pd.DataFrame]) -> pd.DataFrame:
        # Drop empty frames before concatenating; pandas warns about dtype
        # inference when all-NA frames are included.
        frames = [f for f in frames if not f.empty]
        df = pd.concat(frames).reset_index(drop=True) if frames else pd.DataFrame(
            columns=COORD_COLUMNS_2D
        )
        df["substrate"] = self.name
        df["R"], df["G"], df["B"] = self.R, self.G, self.B
        return df

    def coords_2d(self) -> pd.DataFrame:
        raise NotImplementedError


@dataclass
class Rectangle(SurfaceNet):
    """A rectangular net. No rotational supercells are enumerated."""

    a: float = 0.0
    b: float = 0.0

    def __post_init__(self) -> None:
        if self.a > self.b:
            self.a, self.b = self.b, self.a

    def coords_2d(self) -> pd.DataFrame:
        s, angle = 1, "0º"
        rows = [
            [sbl.a, sbl.b, sbl.label, spl.mcia, angle]
            for spl in get_superlattices(self.a, self.b, s)
            for sbl in get_sublattices(self.a, self.b, s, spl, is_square=False)
        ]
        return self._finalize(rows, COORD_COLUMNS_2D)


@dataclass
class Square(SurfaceNet):
    """A square net. Enumerates the commensurate square-lattice rotations."""

    a: float = 0.0

    def coords_2d(self) -> pd.DataFrame:
        frames = []
        for s, (m, n) in SQUARE_ROTATIONS:
            superlattices = get_superlattices(self.a, self.a, s)
            # SQUARE_ROTATIONS is ordered by increasing s, and the enumeration
            # bounds shrink monotonically with s, so once one rotation yields
            # nothing every larger one will too. This break is an early exit,
            # not a skipped case.
            if not superlattices:
                break
            angle = _angle_label(m, n)
            rows = [
                [sbl.a, sbl.b, sbl.label, spl.mcia, angle]
                for spl in superlattices
                for sbl in get_sublattices(self.a, self.a, s, spl, is_square=True)
            ]
            frames.append(pd.DataFrame(rows, columns=COORD_COLUMNS_2D))

        return self._concat(frames)


@dataclass
class RectangleRoot2(SurfaceNet):
    """The a x √2·a net presented by a cubic {110} face."""

    a: float = 0.0

    def coords_2d(self) -> pd.DataFrame:
        b = self.a * math.sqrt(2)
        frames = []
        for s, angle in ((1, "0º"), (3, "35.3º")):
            rows = [
                [sbl.a, sbl.b, sbl.label, spl.mcia, angle]
                for spl in get_superlattices(self.a, b, s)
                for sbl in get_sublattices(self.a, b, s, spl, is_square=False)
            ]
            frames.append(pd.DataFrame(rows, columns=COORD_COLUMNS_2D))

        return self._concat(frames)


@dataclass
class Triangle(SurfaceNet):
    """A triangular net, as presented by cubic (111) or hexagonal (0001).

    Appears in both output tables: ``coords_1d`` gives the triangular matches,
    while ``coords_2d`` recasts it as the centred rectangular cell a x √3·a so
    that rectangular films can match a hexagonal substrate.
    """

    a: float = 0.0

    def coords_2d(self) -> pd.DataFrame:
        rect = Rectangle(self.name, self.a, self.a * math.sqrt(3))
        rect.add_color(self.R, self.G, self.B)
        return rect.coords_2d()

    def coords_1d(self) -> pd.DataFrame:
        # Area of the primitive triangular cell.
        #
        # NOTE: the notebook comment claimed "2 times the simple triangle area,
        # so Area = sqrt(3)*a^2/2" but the code used /4. /4 is preserved here to
        # keep 1D output identical; see docs/OPEN-QUESTIONS.md, which tracks the
        # resulting factor-of-2 offset between 1D and 2D MCIA scales.
        unit_area = math.sqrt(3) * self.a**2 / 4

        rows = []
        for s, angle in TRIANGLE_ROTATIONS:
            area = unit_area * s
            for numer in range(1, math.floor(math.sqrt(MCIA_MAX / area)) + 1):
                mcia = numer**2 * area
                denom_limit = min(
                    math.floor(math.sqrt(s) * self.a * numer / MIN_PARAM) + 1,
                    MAX_AXIS_RATIO + 1,
                )
                for denom in range(1, denom_limit):
                    if any(
                        denom % p == 0 and numer % p == 0 for p in FIRST_FEW_PRIMES
                    ):
                        continue
                    rows.append(
                        [
                            math.sqrt(s) * self.a * numer / denom,
                            format_fraction(numer, denom, s, "a"),
                            mcia,
                            angle,
                        ]
                    )

        return self._finalize(rows, COORD_COLUMNS_1D)
