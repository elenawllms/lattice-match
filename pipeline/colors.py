"""Per-substrate plot colours.

The old scheme swept HSV evenly across ~100 substrate faces. That is not a bad
palette so much as an impossible goal: on a scatter or a map, where any two
marks can sit side by side, a validated categorical palette carries **three**
reliably distinguishable colours. The fourth already fails the normal-vision
floor (yellow vs orange, OKLab ΔE 13.7 against a floor of 15). A hundred hues
at 3.6-degree spacing cannot encode identity, and pretending otherwise is what
made the Voronoi hard to read.

So colour here does a narrower, honest job: **hue names the substrate, lightness
separates its faces**. Silicon's three faces are three steps of one hue rather
than three unrelated colours. That is composite encoding -- the reader gets
"which material" from hue and "which cut" from lightness, and exact identity
from hover, which is the only channel that actually scales to 100 items.

Generated in OKLCH so the steps are perceptually even, unlike HSV, where equal
hue steps are wildly uneven in perceived difference and lightness swings by a
factor of three across the wheel.
"""

from __future__ import annotations

import math
from typing import Sequence

from .nets import SurfaceNet

#: Chroma for generated hues. The method's floor is ~0.10, below which a hue
#: reads as grey and stops doing identity work.
CHROMA = 0.125

#: Lightness steps for the faces of one substrate, inside the light-mode band
#: (OKLCH L 0.43-0.77). The lightest step is 0.62, not the band's 0.77 ceiling:
#: measured across all hues at this chroma, 0.64 is where the worst case first
#: clears 3:1 against the #fcfcfb surface, so anything lighter puts marks below
#: the contrast floor. Steps are spaced ~0.06 apart, which keeps four faces of
#: one material clearly ordered by lightness.
FACE_LIGHTNESS = (0.62, 0.56, 0.50, 0.44)

#: Hue rotation, in degrees, so the first substrate does not land on pure red.
HUE_OFFSET = 25.0


def _oklch_to_srgb(lightness: float, chroma: float, hue_deg: float) -> tuple[float, float, float]:
    """OKLCH -> linear sRGB -> gamma-encoded sRGB, clipped to gamut.

    Ottosson's OKLab matrices. Kept inline rather than pulling in a colour
    library, since this is the only place the pipeline needs them.
    """
    h = math.radians(hue_deg)
    a = chroma * math.cos(h)
    b = chroma * math.sin(h)

    l_ = lightness + 0.3963377774 * a + 0.2158037573 * b
    m_ = lightness - 0.1055613458 * a - 0.0638541728 * b
    s_ = lightness - 0.0894841775 * a - 1.2914855480 * b
    l3, m3, s3 = l_**3, m_**3, s_**3

    r = +4.0767416621 * l3 - 3.3077115913 * m3 + 0.2309699292 * s3
    g = -1.2684380046 * l3 + 2.6097574011 * m3 - 0.3413193965 * s3
    bl = -0.0041960863 * l3 - 0.7034186147 * m3 + 1.7076147010 * s3

    def encode(c: float) -> float:
        c = max(0.0, min(1.0, c))
        return 12.92 * c if c <= 0.0031308 else 1.055 * (c ** (1 / 2.4)) - 0.055

    return encode(r), encode(g), encode(bl)


def assign_colors(nets: Sequence[SurfaceNet]) -> None:
    """Colour each net in place: hue by substrate, lightness by face.

    Substrates are ordered by name so a given material keeps its colour as the
    catalogue grows -- colour follows the entity, not its position in the list.
    """
    # "Silicon (110)" -> "Silicon". Faces of one material share a hue.
    def material(net: SurfaceNet) -> str:
        return net.name.split(" (")[0]

    materials = sorted({material(n) for n in nets})
    hue_step = 360.0 / max(1, len(materials))
    hue_of = {m: (HUE_OFFSET + i * hue_step) % 360.0 for i, m in enumerate(materials)}

    seen: dict[str, int] = {}
    for net in nets:
        m = material(net)
        face_index = seen.get(m, 0)
        seen[m] = face_index + 1
        lightness = FACE_LIGHTNESS[face_index % len(FACE_LIGHTNESS)]
        net.add_color(*_oklch_to_srgb(lightness, CHROMA, hue_of[m]))


def rgb_string(r: float, g: float, b: float) -> str:
    """Format a colour the way the CSV and Plotly expect: ``rgb(255.00, 0.00, 0.00)``."""
    return f"rgb({255 * r:,.2f}, {255 * g:,.2f}, {255 * b:,.2f})"


def relative_luminance(r: float, g: float, b: float) -> float:
    """WCAG relative luminance, for contrast checks in tests."""
    def lin(c: float) -> float:
        return c / 12.92 if c <= 0.04045 else ((c + 0.055) / 1.055) ** 2.4

    return 0.2126 * lin(r) + 0.7152 * lin(g) + 0.0722 * lin(b)
