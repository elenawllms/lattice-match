# Open scientific questions

Things found while porting `update_substrate_list.ipynb` into `pipeline/` that
require a decision from the Falson Lab rather than a code fix. Each was left
**unchanged** so that the port is a faithful one; changing any of them alters
published results and needs a deliberate call.

---

## 1. The A- and R-plane formulas do not match the textbook derivation

**Status:** preserved as-is. Highest priority of the items here.

`pipeline/geometry.py::new_hexagonal_plane` builds these nets:

| Plane | Code produces | Standard derivation |
|---|---|---|
| M (10-10) | `Rectangle(a, c)` | `a × c` — agrees |
| A (11-20) | `Rectangle(a, √(3a² + c²))` | `√3·a × c` |
| R (1-102) | `Rectangle(c, √3·a)` | `√3·a × √(a² + c²)` |

For A-plane, the in-plane directions are [0001] (length `c`) and [1-100]
(length `√3·a`), so the net should be `√3·a × c`. The quantity
`√(3a² + c²)` is the **diagonal** of that rectangle, not an edge.

The R-plane formula `Rectangle(c, √3·a)` is what the A-plane net should be,
which suggests the two may have been transposed at some point.

This is **separate from** the `c`-was-never-passed bug that has already been
fixed (§1.1 of the overhaul plan). Fixing that bug makes A/R/M depend on `c`
for the first time, but if the formulas themselves are wrong then A and R are
still wrong — differently. Because Sapphire A/R/M are heavily used substrates,
this should be resolved before the regenerated data is published.

**Needed:** confirmation of the intended convention, then a test in
`pipeline/tests/test_geometry.py` pinning the agreed values.

---

## 2. CdS and CdSe carry cubic Miller indices on a hexagonal structure

**Status:** excluded, with the reason recorded in `pipeline/data/crystec.csv`.

```
CdS   hexagonal   100, 111, 110   4.14   —   6.76
CdSe  hexagonal   100, 111, 110   4.30   —   7.01
```

`100`, `111` and `110` are cubic indices; a hexagonal entry needs C/M/A/R or
4-index Miller–Bravais notation. The plane column looks copy-pasted from the
cubic rows directly above them in the file.

The notebook removed these with a bare `crystec.drop([11, 12, 21, 22])`, which
also removed LiNbO₃ and LiTaO₃. That drop was a workaround for exactly the four
rows the old code could not process. LiNbO₃/LiTaO₃ are now handled properly
(trigonal support); CdS/CdSe still need real plane data.

**Needed:** the actual planes Crystec sells these in. Then set `include` to
`True` in `pipeline/data/crystec.csv` and the build will pick them up.

---

## 3. (111) faces of tetragonal and orthorhombic crystals are unrepresentable

**Status:** skipped, and reported by `pipeline/build.py`.

Affects `LiAlO2 (111)`, `LiGaO2 (111)`, `TiO2 (111)`, `YAlO3 (111)`.

The (111) face of a tetragonal or orthorhombic crystal is an **oblique** net,
which none of the four supported net types (square, rectangle, √2-rectangle,
triangle) can express. The notebook returned `None` for these and the caller
skipped them silently; they are now reported in `build_report.json`.

**Needed:** either drop these planes from `crystec.csv` as genuinely
unsupported, or add an `Oblique` net type. The latter is a real extension —
the whole superlattice enumeration assumes orthogonal axes.

---

## 4. 1D and 2D MCIA are on different scales

**Status:** preserved. Behaviour-affecting.

`pipeline/nets.py::Triangle.coords_1d` computes the primitive cell area as:

```python
unit_area = math.sqrt(3) * self.a**2 / 4
```

while the notebook's own comment directly above it says:

> Take the unit cell area to be 2 times the simple triangle area.
> Then Area = (sqrt(3)*a^2)/2

The code uses `/4`, the comment says `/2`. So triangular MCIA values are **2×
smaller** than rectangular ones for an equivalent registry.

This matters because the cost function penalises MCIA:
`cost ∝ mcia^(1 - large_superlattice)`. A triangular match therefore gets a
systematically smaller penalty than a rectangular one. Within each table the
ranking is self-consistent, but the two tables cannot be compared, and the
"Allow large superlattices" slider means something slightly different in each.

**Needed:** decide which is intended. If `/2`, the 1D MCIA column changes by a
factor of 2 and previously published 1D rankings shift relative to 2D.

---

## 5. `mismatch()` argument order differs between the table and the map

**Status:** preserved in `src/app.py`; to be unified in the static rewrite.

`mismatch(sub, film)` returns `(film - sub) / sub`.

- `update_table` calls `mismatch(row["a"], a)` — substrate first. Strain is
  normalised by the substrate. This is the conventional definition.
- `update_best_2d_matches` calls `mismatch(a, sp["a"])` — grid point first.
  Strain is normalised by the film, and the sign flips.

So the ranked table and the Voronoi/heatmap compute different quantities. For
small mismatches the argmin is usually the same, so the maps look plausible;
they diverge as mismatch grows.

**Recommendation:** standardise on `(film - sub) / sub`, matching the table
(and the "measure of strain" description in the app's own Learn More text).
This is assumed in the rewrite unless the lab says otherwise.

---

## 6. Monoclinic films support only 2 of the 4 offered planes

**Status:** partially resolved. (001) and (100) are implemented; (1-10) and
(110) are refused.

Materials Project conventional cells for monoclinic materials are in the
standard setting — α = γ = 90°, β the unique angle — confirmed on a 400-material
sample (372 have β off-90; the rest are within rounding of 90°). That makes two
faces rigorously derivable:

| Plane | Spanned by | Angle | Net |
|---|---|---|---|
| (001) | a, b | γ = 90° | `Rectangle(a, b)` |
| (100) | b, c | α = 90° | `Rectangle(b, c)` |
| (010) | a, c | β ≠ 90° | **oblique** |
| (110), (1-10) | a±b, c | c inclined to a | **oblique** |

`pipeline/films.py` additionally verifies α and γ per material and skips the
few that are not in the standard setting (22 of 3,256 in the real run).

The app still offers (1-10) and (110) checkboxes in the monoclinic filter.
They now match nothing. Monoclinic film rows dropped from 9,544 to 5,427 as a
result.

**Needed:** either remove those two checkboxes from the UI, or add an oblique
net type — which is a real extension, since the superlattice enumeration
assumes orthogonal axes (same blocker as item 3).

---

## 7. The substrate colour palette is not distinguishable

**Status:** preserved (HSV sweep), so regenerated data stays comparable.

`pipeline/colors.py` assigns hues evenly around the HSV wheel across ~100
substrate faces. Adjacent faces differ by ~3.6° of hue and are not visually
separable, which directly undermines the Voronoi diagram, whose entire purpose
is telling substrates apart by colour.

**Suggested:** a categorical palette that varies lightness and saturation as
well as hue, ideally grouping faces of the same substrate into one hue family
with different lightness. Purely cosmetic — no effect on rankings.

---

## 8. Triangular films now appear in the 2D table as well as the 1D table

**Status:** changed behaviour. Introduced by the regenerated film catalogue.

A triangular net can match a rectangular substrate through its centred
rectangular cell a × √3·a. `pipeline/nets.py::Triangle.coords_2d` already does
this for *substrates*, and `pipeline/films.py` now does it for *films* too.

The original catalogue was inconsistent about this:

| | old 2D | old 1D | new 2D | new 1D |
|---|---|---|---|---|
| cubic (111) | **0** | 4,785 | 2,066 | 4,753 |
| hexagonal (0001) | 2,393 | 2,820 | 4,181 | 4,961 |

Hexagonal basal films appeared in the 2D table but cubic (111) films did not,
with no apparent reason. The new pipeline treats both the same way.

Consequence: triangular films now show as black dots in *both* the
"Square/Rectangular" and "Triangular" plots. That is physically reasonable —
you can grow a triangular film on a rectangular substrate — but it is a visible
change, and worth confirming it is what you want.

**Needed:** confirm that triangular films belong in both plots. If not, drop
the `rows_2d.append` for `Triangle` in `pipeline/films.py::build_film_nets`.
