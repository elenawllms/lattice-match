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

## Shipping the corrections

*Commit `758938a`. First regeneration of `src/assets/data` from `pipeline/`.*

Substrates, 7,201 → 7,343 rows:

| | before | after |
|---|---|---|
| Sapphire (M) | 4.763 × 4.763 | 4.763 × 13.003 |
| Sapphire (A) | 4.763 × 9.526 | 8.250 × 13.003 |
| Sapphire (R) | 4.763 × 8.250 | 4.763 × 5.133 |
| LiNbO₃, LiTaO₃ | absent | 6 faces |

Films, 73,421 → 74,459 rows. GaN 5.192 → 3.189, ZnO 5.222 → 3.237,
AlN 5.017 → 3.129. Trigonal films 0 → 3,497.

Open question 5 resolved: `mismatch()` is called substrate-first everywhere.

Two testing hazards this created, both closed. The golden diff compared against
`src/assets/data`, which this commit overwrote — left alone it would have
compared the pipeline against itself and passed vacuously forever, so the
pre-correction tables became gzipped fixtures. And `fetch_stable_materials`
would have reused a cache written before the `centring` column existed,
silently treating every R-centred film as primitive; it now validates the
cached schema.

---

## Stage 2 — Binary bundle format

*Commit `d0d1642`.*

A JSON manifest plus a flat `.bin` of concatenated typed arrays, decoded as
zero-copy views. Per-table sidecar metadata, so films stay unfetched until the
film panel opens.

| | CSV | bin | gzipped |
|---|---|---|---|
| superlattices_2d | 0.95 MB | 0.23 MB | **28.9 KB** |
| films_2d | 6.49 MB | 2.83 MB | 945 KB (lazy) |

Initial data load: **34 KB gzipped**, against 5.89 MB of dropdown options alone
before.

Element membership is a 128-bit mask, so the film filters are four bitwise ops
per row instead of six row-wise Python string scans over 81,000 rows.

The round-trip tests caught a bug before it shipped: colour channels were named
by lowercasing R/G/B, so `"B"` collided with the `b` lattice parameter and
overwrote its manifest entry. The browser would have decoded colour bytes as
lattice constants. `BundleWriter` now rejects duplicate column names.

---

## Stage 3 — The browser app

*Commit `5f467f2`.* Plain ES modules with JSDoc types — no build step, no
`node_modules`. Plotly.js gl2d vendored (517 KB gzipped) rather than hotlinked.

The Voronoi is now a single WebGL draw call at 512×512, with the catalogue
uploaded as a float texture and the shader looping over it. **The
1,000-superlattice cap does not exist in this implementation.** Rasters reach
Plotly as a `layout.images` entry backed by an offscreen canvas, which is why
the gl2d partial bundle suffices — no image or heatmap trace module is needed.

New: CSV export of the ranked table, and a film filter that works (it was never
wired as a callback input in the Dash app).

### Two bugs the differential test caught

**float32 was not adequate.** `mismatch()` is `(film − sub) / sub`, a difference
of two nearby numbers, so float32's ~6×10⁻⁸ relative error on `a` is amplified
about 47,000× — roughly 10⁻⁴ on the mismatch itself. Small mismatches, which
are exactly the good matches this tool ranks, lose the most. Two adjacent
entries had already swapped rank. Now float64; the eager payload went 24 → 34 KB
gzipped.

Worth noting the earlier test asserting "float32 is adequate" checked the wrong
quantity: it verified `a` round-trips, not the derived mismatch.

**The ranking was not reproducible.** Binary insertion placed equal-cost entries
*before* existing ones, breaking ties opposite to Python's stable sort. Ties are
ordinary here — InSb and CdTe are both cubic at 6.48 Å, so their (110) faces are
numerically identical.

`pipeline/tests/test_web_differential.py` now pins the two implementations,
running 132 assertions through macOS's built-in JavaScriptCore (no Node needed)
and asserting the assertion count so the harness cannot silently check nothing.

---

## Stage 4 — GitHub Pages

*Commit `8002f38`.* `deploy.yml` uploads `web/` on push to master; nothing is
built. It refuses to deploy on a root-absolute path (which would break under
the `/lattice-match/` base path) or an external CDN reference.

`ci.yml` installs Node so the differential test runs on Linux and **fails if it
skips**, and byte-compares freshly encoded bundles against the committed ones so
`web/data` cannot drift from `src/assets/data`.

Hosting: **$25/mo → $0**.

---

## Status

| Stage | State |
|---|---|
| 0 — emergency patch to the Dash app | done; superseded by the rewrite |
| 1 — pipeline package | done |
| — corrections shipped to production data | done |
| 2 — binary bundle format | done |
| 3 — browser app | done, renders |
| 4 — GitHub Pages + CI | done, **awaiting first push** |

77 tests passing. Nothing has been pushed yet; `origin/master` is still at
`68361ec`.

### To go live

```bash
git checkout master && git merge stage0-emergency-perf-fix && git push origin master
```

Then in the repo settings, set **Pages → Source → GitHub Actions**. The site
lands at `https://elenawllms.github.io/lattice-match/`.

### Where things stand overall

| | before | after |
|---|---|---|
| initial page payload | 5.89 MB of dropdown options alone | 34 KB of data + 517 KB vendored Plotly |
| per-interaction traffic | ~8 MB of JSON | zero — nothing leaves the browser |
| Voronoi solve | 1.8 s, capped at 1,000 superlattices | one GPU draw call, no cap |
| hosting | $25/mo | $0 |
| Sapphire A/M/R | wrong | correct (A still a 3× rectangular approximation) |
| hexagonal films | using `c` instead of `a` | correct |
| trigonal substrates and films | absent | present |
| tests | none | 77 |

### Outstanding decisions

See `docs/OPEN-QUESTIONS.md`. The ones that block or change data:

1. A-plane of R-centred substrates — needs oblique-net support (item 1)
2. CdS/CdSe plane data from Crystec (item 2)
3. (111) of tetragonal/orthorhombic — oblique, unrepresentable (item 3)
4. 1D vs 2D MCIA differ by a factor of 2 (item 4)
5. `mismatch()` argument order differs between table and Voronoi (item 5)
6. Monoclinic (1-10)/(110) checkboxes now match nothing (item 6)
7. Triangular films now appear in the 2D plot as well as 1D (item 8)
