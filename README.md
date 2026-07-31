# Lattice Match

Find commercially available substrates that lattice-match a thin film you want
to grow — including rotated and domain-matched registries, not just 1:1 epitaxy.

Built for the [Falson Lab](https://epitaxy.caltech.edu/) at Caltech.

## What it does

Given a film's in-plane lattice parameters, the tool ranks every commercially
available substrate face by how well its surface net can be made commensurate
with that film.

For each substrate face it precomputes **coincidence cells**: integer stackings
of the substrate net, optionally in a (√s × √s)R θ rotated supercell, whose area
stays under a cap. Each is then subdivided into the film nets that tile it, so a
row of the table reads "3 film cells span 2 substrate cells along a".

This is the coincidence-site-lattice / domain-matching-epitaxy construction of
Zur & McGill, *J. Appl. Phys.* **55**, 378 (1984).

### Scoring

Two quantities determine match quality:

- **Mismatch** — fractional strain per axis, `(film − substrate) / substrate`.
- **MCIA** — minimum coincident interface area, in Å². Larger coincidence cells
  mean a longer repeat distance before registry recurs, so smaller is better.

They combine as

```
cost₂ᴅ = (|εₐ|^p + |ε_b|^p)^(1/p) · MCIA^(1−q)
cost₁ᴅ = |εₐ| · MCIA^(1−q)
```

- `p` = the **single vs. double axis** slider, in [0.05, 1]. At `p → 0.05` the
  quasi-norm blows up unless one axis is nearly perfect, so small `p` rewards a
  single very well-matched axis. At `p = 1` it is a plain sum, rewarding both
  axes being decent.
- `q` = the **allow large superlattices** slider, in [0.05, 1]. At `q = 1` the
  MCIA penalty vanishes entirely; at `q = 0.05` it is close to linear in area.

### Views

- **Scatter** — substrate superlattices (coloured squares, sized by MCIA) and
  selected films (black circles) in (a, b) space. Click a film to score it.
- **Voronoi** — which substrate wins at every point in (a, b) space.
- **Mismatch heatmap** — how good that best match is, everywhere.

The Voronoi and heatmap are the interesting part scientifically: they describe
how to grow a *class* of films, and show which substrates are worth having.

## Layout

```
pipeline/     Python. The science. Runs offline; nothing here serves requests.
  geometry.py       bulk crystal + Miller plane -> 2D surface net
  nets.py           the four net types and their commensurate rotations
  superlattices.py  coincidence-cell enumeration
  films.py          Materials Project fetch -> film catalogue
  build.py          CLI: regenerate the substrate tables
  encode.py         CLI: pack the tables into binary bundles for the web app
  data/crystec.csv  substrate source of truth (hand-maintained)
  tests/
web/          The site. No build step; this directory IS what gets served.
  index.html
  src/core/         cost, ranking, film filtering, the WebGL Voronoi solve
  src/data/         binary bundle decoder
  src/ui/           Plotly figures, styles
  data/             generated binary bundles
  vendor/           Plotly.js gl2d, vendored so there is no CDN dependency
src/          The old Dash app, kept as the reference implementation that the
              differential test checks the browser code against.
docs/         REFACTOR-LOG.md (running change record) + OPEN-QUESTIONS.md
```

## Running

The site is static. Any file server will do:

```bash
cd web && python3 -m http.server 8899     # http://127.0.0.1:8899
```

There is no build step and no `npm install`. Editing a file under `web/src/`
and reloading is the whole development loop. `web/jsconfig.json` gives VS Code
full type checking of the JSDoc annotations without a compiler.

The old Dash app still runs, if you need to compare against it:

```bash
.venv/bin/python src/app.py               # http://127.0.0.1:8050
```

## Regenerating the data

Substrate superlattices, from `pipeline/data/crystec.csv`:

```bash
.venv/bin/python -m pipeline.build --out src/assets/data
```

This writes `sublattices_2d.csv`, `sublattices_1d.csv` and a
`build_report.json` listing every substrate face that could not be represented.
Read that report — it is how you find out that something was dropped.

Film catalogue, from the Materials Project:

```bash
export MP_API_KEY=...        # https://next-gen.materialsproject.org/api
.venv/bin/python -m pipeline.films --out src/assets/data
```

The key is needed only at build time and never ships. A key was previously
hardcoded in `src/app.py` and remains recoverable from git history — treat it
as compromised.

## Deployment

Pushing to `master` publishes `web/` to GitHub Pages via
`.github/workflows/deploy.yml`. Nothing is built; the directory is uploaded
as-is. The workflow refuses to deploy if the page gains a root-absolute path
(which would break under the `/lattice-match/` base path) or a reference to an
external CDN.

Hosting cost is zero.

## Tests

```bash
.venv/bin/python -m pytest -q
```

Three groups matter:

- `test_build.py` golden-diffs regenerated tables against the originally
  shipped CSVs, in `pipeline/tests/fixtures/`, and asserts that **only** the
  documented corrections changed anything. That is what made it safe to
  replace production data.
- `test_encode.py` round-trips the binary bundles, checking byte alignment,
  non-overlapping columns and the element-bitmask semantics -- all things that
  fail silently rather than loudly in that format.
- `test_web_differential.py` runs the browser's cost and ranking code against
  the Python implementation over the real catalogue. The physics exists twice
  now, and two implementations drifting apart is exactly how the old app ended
  up with a table and a map that disagreed. It runs under macOS's built-in
  JavaScriptCore locally and Node in CI, and CI fails if it skips.

## Documentation

- `docs/REFACTOR-LOG.md` — running record of the 2026 overhaul: what changed,
  why, and what it measured. Start here.
- `docs/OPEN-QUESTIONS.md` — items needing a scientific decision rather than a
  code change, most importantly the A-plane of rhombohedrally centred
  substrates, whose true surface mesh is oblique and cannot currently be
  represented.
