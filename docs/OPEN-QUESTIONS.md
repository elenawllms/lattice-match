# Open scientific questions

Things found while porting `update_substrate_list.ipynb` into `pipeline/` that
require a decision from the Falson Lab rather than a code fix. Each was left
**unchanged** so that the port is a faithful one; changing any of them alters
published results and needs a deliberate call.

---

## 1. The A- and R-plane formulas were transposed

**Status: RESOLVED for primitive hexagonal lattices.** The R-centring caveat
below is still open.

Both formulas in `new_hexagonal_plane` were individually correct but attached
to the wrong plane labels. Derived by enumerating the lattice points that lie
in each plane through the origin and reducing to a shortest basis, for a
**primitive hexagonal** lattice with a = 4.763, c = 13.003 (sapphire's
hexagonal cell):

| Plane | In-plane vectors | True net | Old code |
|---|---|---|---|
| M (10-10) | a₂, c | `a × c` = 4.763 × 13.003 | correct |
| A (11-20) | a₁−a₂, c | `√3·a × c` = **8.250 × 13.003** | gave 4.763 × 15.399 — the R net |
| R (1-102) | a₁+a₂, −a₁+a₂+c | `a × √(3a²+c²)` = **4.763 × 15.399** | gave 8.250 × 13.003 — the A net |

All three come out at exactly 90°, confirming they are genuine rectangular
meshes, so the enumeration's orthogonal-axes assumption is sound here.

The branches are now swapped. Pinned by `test_a_and_r_are_not_transposed` and
`test_hexagonal_nets_are_rectangular_for_a_primitive_lattice`, the latter
checking each net against the in-plane vectors directly rather than against a
hardcoded formula.

This was **separate from** the `c`-was-never-passed bug (§1.1 of the overhaul
plan). That bug made A/R/M collapse to c = a; fixing it alone would have left A
and R still wrong — simply each other.

### The R-centred case, resolved

The derivation above assumes a **primitive hexagonal** lattice — correct for
wurtzites like GaN, ZnO and AlN (P6₃mc). Sapphire, LiNbO₃ and LiTaO₃ are
R-3c/R3c, whose lattice is **rhombohedrally centred**, with extra points at
(2/3,1/3,1/3) and (1/3,2/3,2/3).

Recomputed with the centring included, in a single Cartesian frame, and checked
against `area × d_spacing = V_primitive` — an identity any correct 2D mesh must
satisfy. It holds exactly for every plane below.

| Plane | Wurtzite (P6₃mc) | R-centred (R-3c / R3c) |
|---|---|---|
| C (0001) | `a`, triangular | **same** |
| M (10-10) | `a × c` | **same** |
| A (11-20) | `√3·a × c` | **truly oblique**, 84–86° |
| R (1-102) | `a × √(3a²+c²)` | **`a × √(3a²+c²) / 3`** |

Measured, for the three R-centred substrates:

| | a | c | R-plane true | `a × √(3a²+c²)/3` | A-plane true |
|---|---|---|---|---|---|
| Al₂O₃ | 4.805 | 13.116 | 4.805 × 5.178, 90° | 4.805 × 5.178 ✓ | 5.178 × 7.064, 84.2° |
| LiNbO₃ | 5.269 | 13.903 | 5.269 × 5.544, 90° | 5.269 × 5.544 ✓ | 5.544 × 7.649, 86.0° |
| LiTaO₃ | 5.134 | 13.816 | 5.134 × 5.477, 90° | 5.134 × 5.477 ✓ | 5.477 × 7.507, 84.9° |

So, concretely:

**R-plane** — the current formula overstates the long axis by exactly 3×. The
true mesh is rectangular, so this is fixable exactly: divide by 3 when the
lattice is R-centred.

**A-plane** — the true mesh is oblique and therefore cannot be represented at
all by the superlattice enumeration, which assumes orthogonal axes. The
smallest *rectangular* sublattice is `√3·a × c`, exactly 3× the primitive area
— which is what the code already produces. That is a legitimate, if
conservative, choice: a valid sublattice, just not the primitive one.

### Why the 3× matters

Both A and R currently use a cell 3× larger in area than the true surface mesh.
Two consequences, neither fatal:

1. **MCIA is overstated 3×**, so those faces are over-penalised in ranking. At
   the default slider (`q = 0.5`) the penalty goes as `MCIA^0.5`, so they are
   ranked ~1.7× worse than they deserve.
2. **Coverage is lost.** `get_superlattices` caps enumeration at
   `MCIA_MAX = 200 Å²`. For sapphire R-plane the true mesh is 24.9 Å², allowing
   stackings up to the `MAX_AXIS_RATIO = 5` limit; the 3× cell is 74.6 Å²,
   allowing only 2. Valid matches are never enumerated.

The geometry itself is still *found* — `get_sublattices` divides by integer
denominators up to 5, so the `(a × b/3)` entry does appear with the right
dimensions. It just carries the wrong MCIA.

### Implemented

`crystec.csv` now carries a `centring` column (`P`/`R`), and
`new_hexagonal_plane` takes a `centring` argument. Films derive it from the
Materials Project Hermann-Mauguin symbol, whose first letter is the lattice
type. Sapphire, LiNbO₃ and LiTaO₃ are marked `R`; ZnO is `P`.

Effect on Sapphire (R), at the 1:1 registry: 4.763 × 15.399 → **4.763 × 5.133**,
matching the measured mesh. Its superlattice count rose from 36 to 112, because
the true 24.9 Å² cell admits far more stackings under `MCIA_MAX = 200 Å²` than
the 74.6 Å² one did — the coverage loss predicted above, recovered.

**Still outstanding:** the A-plane of R-centred crystals keeps the 3× rectangular
sublattice, because its true mesh is oblique and cannot be represented. Those
faces still report 3× the true coincident area and are under-ranked accordingly.
Fixing it properly requires oblique-net support in the superlattice enumeration
— the same blocker as item 3. Affects Sapphire (A), LiNbO₃ (1120), LiTaO₃ (1120).

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
