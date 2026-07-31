"""Properties the substrate palette must hold.

The previous scheme swept HSV across ~100 faces. Measured against the data-viz
method's validator, a categorical palette carries three reliably distinguishable
colours on an all-pairs form like a scatter -- the fourth slot already fails the
normal-vision floor. So identity by hue alone was never achievable here, and
these tests pin what the replacement actually promises instead: hue groups a
material, lightness orders its faces, and every mark is legible on the surface.
"""

from __future__ import annotations

import colorsys

import pandas as pd
import pytest

from pipeline.colors import (
    CHROMA,
    FACE_LIGHTNESS,
    _oklch_to_srgb,
    assign_colors,
    relative_luminance,
    rgb_string,
)
from pipeline.nets import Square

SURFACE = (252 / 255, 252 / 255, 251 / 255)
MIN_CONTRAST = 3.0


def contrast(rgb, other=SURFACE) -> float:
    a, b = relative_luminance(*rgb), relative_luminance(*other)
    hi, lo = max(a, b), min(a, b)
    return (hi + 0.05) / (lo + 0.05)


@pytest.fixture(scope="module")
def shipped():
    df = pd.read_csv("src/assets/data/sublattices_2d.csv")
    return df.groupby("substrate")[["R", "G", "B"]].first()


def test_every_face_clears_the_contrast_floor(shipped):
    """Marks below 3:1 on the surface are the one thing the method will not
    let you dismiss without a relief channel."""
    worst = min(contrast(tuple(row)) for row in shipped.values)
    assert worst >= MIN_CONTRAST, f"worst face contrast {worst:.2f}:1"


def test_faces_of_one_material_share_a_hue(shipped):
    """Hue names the material; it must not drift across its own faces."""
    by_material: dict[str, list[float]] = {}
    for name, row in shipped.iterrows():
        h, _s, _v = colorsys.rgb_to_hsv(row.R, row.G, row.B)
        by_material.setdefault(name.split(" (")[0], []).append(h * 360)
    multi = {m: hs for m, hs in by_material.items() if len(hs) > 1}
    assert multi, "expected at least one material with several faces"
    for material, hues in multi.items():
        assert max(hues) - min(hues) < 5.0, f"{material} hues spread {max(hues) - min(hues):.1f} deg"


def test_faces_of_one_material_differ_in_lightness(shipped):
    """Lightness is what separates faces once hue is spent on the material."""
    by_material: dict[str, list[float]] = {}
    for name, row in shipped.iterrows():
        by_material.setdefault(name.split(" (")[0], []).append(
            relative_luminance(row.R, row.G, row.B)
        )
    for material, lums in by_material.items():
        if len(lums) < 2:
            continue
        lums.sort()
        gaps = [b - a for a, b in zip(lums, lums[1:])]
        assert min(gaps) > 0.02, f"{material} faces separated by only {min(gaps):.3f} luminance"


def test_lightness_steps_stay_in_the_light_mode_band():
    assert all(0.43 <= step <= 0.77 for step in FACE_LIGHTNESS)
    assert list(FACE_LIGHTNESS) == sorted(FACE_LIGHTNESS, reverse=True), "steps must descend"


def test_chroma_clears_the_identity_floor():
    """Below ~0.10 a hue reads as grey and stops doing identity work."""
    assert CHROMA >= 0.10


def test_lightest_step_is_the_binding_constraint():
    """Documents why the lightest step is 0.62 and not the band ceiling."""
    too_light = min(contrast(_oklch_to_srgb(0.70, CHROMA, h)) for h in range(0, 360, 15))
    chosen = min(contrast(_oklch_to_srgb(FACE_LIGHTNESS[0], CHROMA, h)) for h in range(0, 360, 15))
    assert too_light < MIN_CONTRAST <= chosen


def test_colour_follows_the_substrate_not_its_position():
    """Adding a substrate must not repaint the others' hue ordering rule: a
    material's hue depends on its name's rank, so the mapping is stable for a
    fixed catalogue and never depends on row order."""
    a = [Square(f"{m} (100)") for m in ("Silicon", "MgO", "SrTiO3")]
    b = [Square(f"{m} (100)") for m in ("SrTiO3", "Silicon", "MgO")]
    assign_colors(a)
    assign_colors(b)
    by_name_a = {n.name: (n.R, n.G, n.B) for n in a}
    by_name_b = {n.name: (n.R, n.G, n.B) for n in b}
    assert by_name_a == by_name_b


def test_rgb_string_format_is_unchanged():
    assert rgb_string(1.0, 0.0, 0.5) == "rgb(255.00, 0.00, 127.50)"
