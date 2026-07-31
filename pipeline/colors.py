"""Per-substrate plot colours.

Each substrate face gets a distinct colour, used for the scatter markers and
for the Voronoi partition fill.

The HSV sweep is carried over from the notebook so that regenerated data is
directly comparable to the shipped CSVs. It is not a good palette: at ~90
faces, adjacent hues are visually indistinguishable, which undermines the
Voronoi's whole purpose. Replacing it is tracked in docs/OPEN-QUESTIONS.md
and deliberately kept separate from the physics corrections.
"""

from __future__ import annotations

from typing import Sequence

from matplotlib import colormaps

from .nets import SurfaceNet


def assign_colors(nets: Sequence[SurfaceNet]) -> None:
    """Assign each net an evenly spaced hue, in place.

    Uses ``matplotlib.colormaps`` rather than the notebook's ``cm.get_cmap``,
    which was removed in matplotlib 3.9.
    """
    cmap = colormaps["hsv"]
    n = len(nets)
    for i, net in enumerate(nets):
        r, g, b, _alpha = cmap(i / n)
        net.add_color(r, g, b)


def rgb_string(r: float, g: float, b: float) -> str:
    """Format a colour the way the CSV and Plotly expect: ``rgb(255.00, 0.00, 0.00)``."""
    return f"rgb({255 * r:,.2f}, {255 * g:,.2f}, {255 * b:,.2f})"
