# Refactor log

A running record of the 2026 overhaul: what changed, why, and what it measured.
Newest entries at the bottom. Companion documents:

- `docs/OPEN-QUESTIONS.md` — scientific questions still needing a lab decision
- `README.md` — how to run and regenerate everything

**Branch:** `stage0-emergency-perf-fix` (not yet merged or deployed)

---

## Why this started

The tool rendered extremely slowly, sometimes crashed the browser, and needed a
$25/mo Render instance. The diagnosis found that the slowness was entirely
architectural — the actual math is tiny — and, separately, that the tool had
been returning **scientifically wrong answers** for several heavily used
substrates.

Plan: five stages. Stage 0 stops the bleeding on the live app; stage 1 rebuilds
the science as a tested package; stages 2–4 replace the Dash server with a
static site.

---

## Stage 0 — Emergency patch to the live Dash app

*Commit `55fb26a`. Goal: make the free tier viable without refactoring.*

### The boot bug (the expensive one)

`app.run_server()` was called at **module import time** with no `__main__`
guard. `render.yaml` starts the app with `gunicorn app:server`, and gunicorn
imports the module to find `server` — an import that never returned, because it
blocked forever inside Flask's development server. The site was very likely
being served by Werkzeug, single-threaded, with gunicorn workers looping. That
is the main reason a paid instance was needed. Fixed with a `__main__` guard.

### The first-paint killer

`dcc.Dropdown(films, ...)` seeded **78,417 film names — 5.89 MB of JSON** into
the initial page payload, mounted into React-Select before the user clicked
anything. That dropdown was never a callback `Input`, so the options did
nothing. `pull_materials` now supplies them alongside the values it already set.

### The wasted computation

`update_best_2d_matches` had no `display-mode` input, so the Voronoi grid solve
ran on every substrate change and every slider release — including on the
default Scatter tab, where the result is discarded. Now gated.

The solve itself used `np.vectorize`, which is a Python for-loop: 10,000 calls,
each running ~8 pandas Series ops and evaluating the cost function **twice**.
Replaced with a chunked broadcast (chunked to bound peak memory on a 512 MB
instance). Verified **bit-identical argmin** against the old implementation;
`minCost` agrees to within 1 ULP.

### Everything else in this commit

| Change | Reason |
|---|---|
| `go.Scatter` → `go.Scattergl` ×4 | default view drew 7,201 SVG nodes |
| heatmap `x`/`y` `150` → `RESOLUTION` | `z` was 100×100; axes were misregistered |
| `update_table` returns `rows[:200]` | returned all 7,201 rows to a 10-row table |
| unified the `>1000` threshold | guard and warning used opposite comparisons |
| element filters test falsiness | cleared multi-select gives `[]`, not `None`, so searches silently returned nothing |
| guarded `int(num_elements)` | field renders as free text (`type="numeric"` is not valid Dash) |
| `pathlib` instead of `/opt/render/project/src/` | app could not run locally at all |
| pinned all dependencies | `dash>=2.0.0` now resolves to Dash 3.x, which removed `run_server` |
| `compress=True` | Dash defaults it to **False**; installing Flask-Compress alone does nothing |
| gunicorn `--bind $PORT` | gunicorn's own default is `127.0.0.1:8000` |
| dropped `stable_films.csv`, `.ipynb_checkpoints`, `costToColor`, unused `ctx` | dead weight |

### Measured, worst case (no substrate filter, all 7,201 superlattices)

| | before | after |
|---|---|---|
| initial layout | 5.89 MB | **4.2 KB** gzipped |
| 2D figure | 2.0 MB | **80 KB** gzipped, 249 ms |
| matches table | 1.8 MB | **3 KB** gzipped, 26 ms |
| Voronoi solve | 1.8 s | **34 ms**, skipped entirely on Scatter |
| gunicorn | worker restart loop | boots once per worker, 145 MB total |

**Known remaining:** each callback still uploads ~1.7 MB, because the
superlattice catalogue round-trips through `dcc.Store`. That is architectural
and is what the static rewrite removes.

---

## Stage 1 — The science, extracted and tested

*Commit `b9191a3`. `pipeline/` replaces `update_substrate_list.ipynb`.*

The notebook could not be run at all — cell 39 is a bare `sublattices_2d[]`, a
`SyntaxError`. Ported to real modules with type hints and tests.

### Physics bugs found and fixed

**1. Hexagonal planes ignored `c`.** The notebook declared
`newHexagonalPlane(name, c, plane)` but called it as `(name, a, plane)`, so the
parameter named `c` was bound to `a` while the body read a module-global `a`.
Every hexagonal non-basal plane was computed as though `c == a`.

**2. Four substrates were silently missing.** `new_plane` returned `None` for
unhandled cases and the caller skipped it. Combined with a bare
`crystec.drop([11, 12, 21, 22])` — a hardcoded positional drop with no comment
— this removed LiNbO₃, LiTaO₃, CdS and CdSe. 41 substrates in the source, 37
reaching the app, no warning anywhere.

Now: trigonal is supported (recovering LiNbO₃/LiTaO₃), `new_plane` **raises**
instead of returning `None`, and `build.py` records every skipped face in
`build_report.json`. The positional drop is replaced by explicit
`include`/`notes` columns. CdS and CdSe stay excluded with a stated reason —
their `Plane` column holds cubic Miller indices on a hexagonal structure.

**3. `cm.get_cmap`** was removed in matplotlib 3.9.

### Golden diff

`pipeline/tests/test_build.py` regenerates the tables and asserts that **only**
the documented corrections changed anything:

- 2D: 7,201 → 7,307 rows. Geometry changed on exactly Sapphire A, M, R. Six
  faces added, all LiNbO₃/LiTaO₃. Nothing lost.
- 1D: 336 → 384 rows. Zero geometry changes; only the two trigonal basal faces.

This is what makes it safe to replace production data.

Also added `pipeline/films.py`, restoring the Materials Project fetch deleted in
`ffb6432`, reading `MP_API_KEY` from the environment instead of hardcoding it.

---

## Film catalogue regenerated

*Commit `a5eb503`. First real run against the Materials Project: 24,241 stable
materials with ≤3 elements.*

### A bug caught before it shipped

MP's `summary.structure` field is the **primitive** cell. Using it directly
would have made every cubic film wrong by √2 — primitive fcc silicon is
a = 3.849 Å, but the value this tool needs is the conventional 5.444 Å. Now
reduced via `SpacegroupAnalyzer.get_conventional_standard_structure()`. Cubic
output is **bit-identical** to the original catalogue, which is the strongest
evidence the port is faithful.

### The hexagonal bug, confirmed on real data

| | old catalogue | new | literature |
|---|---|---|---|
| GaN | 5.192 | **3.189** | 3.19 |
| ZnO | 5.222 | **3.237** | 3.25 |
| AlN | 5.017 | **3.129** | 3.11 |
| Mg | 5.141 | **3.172** | 3.21 |
| Zn | 4.873 | **2.614** | 2.66 |

Every one had been using `c`. `validate_known_films` now asserts these so the
regression cannot ship again.

### Validation

**48,553 of 48,675 non-hexagonal pairs identical** (99.75%) once axis ordering
is normalised. The 122 that differ are MP data drift since the Jan 2025
snapshot. Raw a/b ordering differs more widely, but only because `Rectangle`
normalises to a ≤ b and the original did not — same net.

| system | old rows | new rows | why |
|---|---|---|---|
| Hexagonal | 9,672 | 9,913 | the c-for-a fix |
| Trigonal | **0** | **3,227** | absent entirely before |
| Monoclinic | 9,544 | 5,427 | 2 supported planes instead of 4 |
| Cubic | 10,075 | 12,077 | more stable entries in MP now |

### Monoclinic support added

MP conventional cells use the standard setting — α = γ = 90°, β unique —
confirmed on a 400-material sample (372 with β off-90). That makes two faces
rigorously derivable: (001) → `Rectangle(a, b)` and (100) → `Rectangle(b, c)`.
(010) and the {110} faces are genuinely oblique and are refused. `films.py`
verifies α and γ per material and skips the 22 of 3,256 not in standard setting.

---

## Hexagonal A- and R-planes: transposed, then centred

*Commits `2558958` and `1a0ceae`.*

### The transposition

Both expressions were individually correct but attached to the **wrong plane
labels**. Derived by enumerating the lattice points lying in each plane through
the origin and reducing to a shortest basis:

| Plane | In-plane vectors | True net |
|---|---|---|
| M (10-10) | a₂, c | `a × c` |
| A (11-20) | a₁−a₂, c | `√3·a × c` |
| R (1-102) | a₁+a₂, −a₁+a₂+c | `a × √(3a²+c²)` |

Sapphire, at the 1:1 registry:

```
face           shipped         after swap       derived truth
Sapphire (M)   4.763 ×  4.763  4.763 × 13.003   4.763 × 13.003
Sapphire (A)   4.763 ×  9.526  8.250 × 13.003   8.250 × 13.003
Sapphire (R)   4.763 ×  8.250  4.763 × 15.399   4.763 × 15.399
```

`test_hexagonal_nets_are_rectangular_for_a_primitive_lattice` checks each net
against the in-plane vectors themselves rather than a hardcoded formula, so it
verifies the derivation instead of restating it.

### Rhombohedral centring

The above assumes a **primitive** hexagonal lattice — right for wurtzites, wrong
for sapphire, LiNbO₃ and LiTaO₃, which are R-3c/R3c and carry extra lattice
points at (2/3,1/3,1/3) and (1/3,2/3,2/3).

Two dead ends before getting this right, both worth recording:

- pymatgen's `SlabGenerator` given the *conventional* cell returns the
  hexagonal-lattice answer, because in that cell the centring translations are
  basis atoms rather than lattice translations.
- Taking the plane normal from the conventional cell and lattice vectors from
  the primitive cell mixes two different Cartesian orientations.

Settled by computing in a single frame and cross-checking against
`area × d_spacing = V_primitive`, an identity any correct 2D mesh must satisfy.

| Plane | Wurtzite (P6₃mc) | R-centred (R-3c / R3c) |
|---|---|---|
| C (0001) | `a`, triangular | same |
| M (10-10) | `a × c` | same |
| A (11-20) | `√3·a × c` | **truly oblique**, 84–86° |
| R (1-102) | `a × √(3a²+c²)` | **`a × √(3a²+c²)/3`** |

| | a | c | R-plane measured | `a×√(3a²+c²)/3` |
|---|---|---|---|---|
| Al₂O₃ | 4.805 | 13.116 | 4.805 × 5.178 | 4.805 × 5.178 ✓ |
| LiNbO₃ | 5.269 | 13.903 | 5.269 × 5.544 | 5.269 × 5.544 ✓ |
| LiTaO₃ | 5.134 | 13.816 | 5.134 × 5.477 | 5.134 × 5.477 ✓ |

**Implemented:** `crystec.csv` gained a `centring` column (`P`/`R`);
`new_hexagonal_plane` takes a `centring` argument; films derive it from the MP
Hermann-Mauguin symbol's first letter.

Sapphire (R) went from 4.763 × 15.399 to **4.763 × 5.133**, and its superlattice
count from **36 to 112** — the true 24.9 Å² cell admits far more stackings under
`MCIA_MAX = 200 Å²` than the 74.6 Å² one did. That is real coverage that had
been silently unavailable.

**Not fixed:** the A-plane of R-centred crystals. Its true mesh is oblique and
cannot be represented, because the superlattice enumeration assumes orthogonal
axes. `√3·a × c` is the smallest rectangular sublattice, at exactly 3× the
primitive area — valid, but it reports 3× the true coincident area, so those
faces are under-ranked by ~1.7× at the default slider. Needs oblique-net
support; tracked as open question 1.

---

## Status

| Stage | State |
|---|---|
| 0 — emergency patch | done, verified, **not deployed** |
| 1 — pipeline package | done, 61 tests passing |
| — film catalogue | regenerated and validated, **not shipped** |
| 2 — binary data format | not started |
| 3 — static TypeScript app | not started |
| 4 — deploy, retire Render | not started |

### Production data has deliberately not been regenerated

`src/assets/data/*.csv` is still the original. The corrections are ready to
apply, but shipping them changes published scientific results, so it should be
a deliberate act:

```bash
python -m pipeline.build --out src/assets/data
MP_API_KEY=... python -m pipeline.films --out src/assets/data --cache .mp_cache.parquet
```

### Outstanding decisions

See `docs/OPEN-QUESTIONS.md`. The ones that block or change data:

1. A-plane of R-centred substrates — needs oblique-net support (item 1)
2. CdS/CdSe plane data from Crystec (item 2)
3. (111) of tetragonal/orthorhombic — oblique, unrepresentable (item 3)
4. 1D vs 2D MCIA differ by a factor of 2 (item 4)
5. `mismatch()` argument order differs between table and Voronoi (item 5)
6. Monoclinic (1-10)/(110) checkboxes now match nothing (item 6)
7. Triangular films now appear in the 2D plot as well as 1D (item 8)
