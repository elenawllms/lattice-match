// @ts-check
/**
 * Wiring: data in, state, plots and table out.
 *
 * Everything runs in the browser. There is no server, and no interaction
 * sends a byte anywhere. The Dash app moved roughly 8 MB of JSON per slider
 * drag to compute exactly what this file computes in under a millisecond.
 */

import { loadEager } from './data/bundle.js';
import { best1d, best2d } from './core/match.js';
import { filterFilms } from './core/filter.js';
import { VoronoiSolver, solveOnCpu } from './core/voronoi.js';
import { VIEW_2D, render1d, render2d } from './ui/plots.js';

const DATA_BASE = 'data';
const TABLE_LIMIT = 200;
const PAGE_SIZE = 10;
const MAX_FILMS = 1000;

const PLANES_BY_SYSTEM = {
  Cubic: ['(100)', '(110)', '(111)'],
  Tetragonal: ['(001)', '(100)', '(101)', '(110)'],
  Orthorhombic: ['(001)', '(010)', '(100)', '(011)', '(101)', '(110)'],
  Monoclinic: ['(001)', '(100)'],
  Hexagonal: ['(0001)', '(1-100)', '(1-102)', '(11-20)'],
  Trigonal: ['(0001)', '(1-100)', '(11-20)'],
};

const $ = (/** @type {string} */ id) => /** @type {HTMLElement} */ (document.getElementById(id));

const state = {
  /** @type {Set<string>} */ substrates: new Set(),
  dimension: 2,
  a: 5,
  b: 5,
  singleDouble: 0.5,
  largeSuperlattice: 0.5,
  mode: 'Scatter',
  page: 0,
  /** @type {any[]} */ matches: [],
  /** @type {number[]} */ filmRows2d: [],
  /** @type {number[]} */ filmRows1d: [],
  selectedFilm: 'Manual Entry',
};

/** @type {any} */ let data;
/** @type {any} */ let solver;
/** @type {boolean|undefined} */ let webglOk;
/** Row indices of the currently selected substrates, per table. */
const visible = { '2d': /** @type {number[]} */ ([]), '1d': /** @type {number[]} */ ([]) };
let rasterDirty = true;

async function main() {
  try {
    data = await loadEager(DATA_BASE);
  } catch (err) {
    const el = $('loading');
    el.classList.add('error');
    el.textContent = `Could not load data: ${/** @type {Error} */ (err).message}`;
    return;
  }

  buildSubstrateList();
  buildFilmFilters();
  wireControls();

  selectAllSubstrates();
  $('loading').hidden = true;
  refresh({ substrates: true });
}

/**
 * WebGL contexts are a limited resource, and Plotly's scattergl needs its own
 * for each plot. Creating the solver's context up front -- before Plotly had
 * rendered anything -- risked exhausting that budget and making the scatter
 * traces fail with "WebGL is not supported". Build it only when a raster mode
 * is actually selected.
 */
function getSolver() {
  if (solver === undefined) solver = new VoronoiSolver(512);
  return solver;
}

/** True when the browser can drive Plotly's WebGL scatter traces. */
function webglAvailable() {
  if (webglOk === undefined) {
    try {
      const probe = document.createElement('canvas');
      webglOk = Boolean(probe.getContext('webgl') || probe.getContext('experimental-webgl'));
    } catch {
      webglOk = false;
    }
  }
  return webglOk;
}

/* ------------------------------------------------------------------ data */

function superlattices(dim) {
  return data.tables[dim === 1 ? 'superlattices_1d' : 'superlattices_2d'];
}

/** Recompute which rows belong to the selected substrates. */
function recomputeVisible() {
  for (const dim of ['2d', '1d']) {
    const table = data.tables[`superlattices_${dim}`];
    const rows = [];
    for (let i = 0; i < table.rows; i++) {
      if (state.substrates.has(table.label('substrate', i))) rows.push(i);
    }
    visible[dim] = rows;
  }
  rasterDirty = true;
}

/* --------------------------------------------------------------- compute */

function recomputeMatches() {
  const dim = state.dimension;
  const table = superlattices(dim);
  const mask = new Uint8Array(table.rows);
  for (const i of visible[dim === 1 ? '1d' : '2d']) mask[i] = 1;

  state.matches = dim === 1
    ? best1d(table, state.a, state.largeSuperlattice, TABLE_LIMIT, mask)
    : best2d(table, state.a, state.b, state.singleDouble, state.largeSuperlattice, TABLE_LIMIT, mask);
  state.page = 0;
}

function updateRaster() {
  if (state.mode === 'Scatter') return null;
  const table = data.tables.superlattices_2d;
  const rows = visible['2d'];
  if (!rows.length) return null;

  const limits = /** @type {[number,number,number,number]} */ (
    [VIEW_2D.aMin, VIEW_2D.aMax, VIEW_2D.bMin, VIEW_2D.bMax]);
  const mode = state.mode === 'Voronoi' ? 'voronoi' : 'heatmap';

  const gpu = getSolver();
  if (gpu.supported) {
    if (rasterDirty) {
      const a = new Float32Array(rows.length);
      const b = new Float32Array(rows.length);
      const mcia = new Float32Array(rows.length);
      const rgb = new Uint8Array(rows.length * 3);
      const ca = table.col('a'), cb = table.col('b'), cm = table.col('mcia');
      const cr = table.col('red'), cg = table.col('green'), cbl = table.col('blue');
      rows.forEach((i, k) => {
        a[k] = ca[i]; b[k] = cb[i]; mcia[k] = cm[i];
        rgb[k * 3] = cr[i]; rgb[k * 3 + 1] = cg[i]; rgb[k * 3 + 2] = cbl[i];
      });
      gpu.upload(a, b, mcia, rgb);
      rasterDirty = false;
    }
    // The heatmap ramp needs a scale; the worst cost in the current view is a
    // stable choice and keeps the colouring comparable as sliders move.
    const scale = mode === 'heatmap' ? heatmapScale(table, rows) : 1;
    return gpu.render(limits, state.singleDouble, state.largeSuperlattice, mode, scale);
  }

  // CPU fallback: lower resolution, identical maths.
  const a = new Float32Array(rows.length), b = new Float32Array(rows.length);
  const mcia = new Float32Array(rows.length), rgb = new Uint8Array(rows.length * 3);
  const ca = table.col('a'), cb = table.col('b'), cm = table.col('mcia');
  const cr = table.col('red'), cg = table.col('green'), cbl = table.col('blue');
  rows.forEach((i, k) => {
    a[k] = ca[i]; b[k] = cb[i]; mcia[k] = cm[i];
    rgb[k * 3] = cr[i]; rgb[k * 3 + 1] = cg[i]; rgb[k * 3 + 2] = cbl[i];
  });
  return solveOnCpu({ a, b, mcia, rgb }, limits, state.singleDouble, state.largeSuperlattice, mode);
}

/** A representative worst-case cost, used to normalise the heatmap ramp. */
function heatmapScale(table, rows) {
  const m = best2d(table, VIEW_2D.aMax, VIEW_2D.bMax, state.singleDouble,
                   state.largeSuperlattice, 1, maskOf(table, rows));
  return m.length ? Math.max(m[0].cost, 1e-6) * 4 : 1;
}

function maskOf(table, rows) {
  const mask = new Uint8Array(table.rows);
  for (const i of rows) mask[i] = 1;
  return mask;
}

/* ----------------------------------------------------------------- render */

/** @param {{substrates?: boolean}} [opts] */
function refresh(opts = {}) {
  if (opts.substrates) recomputeVisible();
  recomputeMatches();
  renderTable();
  renderReadouts();

  const raster = updateRaster();
  render2d($('figure-2d'), {
    superlattices: data.tables.superlattices_2d,
    rows: visible['2d'],
    films: data.tables.films_2d,
    filmRows: state.filmRows2d,
    largeSuperlattice: state.largeSuperlattice,
    mode: state.mode,
    raster,
    webgl: webglAvailable(),
  }).then(attachClick2d);

  render1d($('figure-1d'), {
    superlattices: data.tables.superlattices_1d,
    rows: visible['1d'],
    films: data.tables.films_1d,
    filmRows: state.filmRows1d,
    largeSuperlattice: state.largeSuperlattice,
    mode: state.mode,
    webgl: webglAvailable(),
  }).then(attachClick1d);
}

function renderTable() {
  const dim = state.dimension;
  const table = superlattices(dim);
  const head = $('matches-head');
  const body = $('matches-body');

  head.innerHTML = dim === 1
    ? '<th>Substrate</th><th>Dimensions</th><th>a mismatch</th><th>MCIA</th>'
    : '<th>Substrate</th><th>Dimensions</th><th>a mismatch</th><th>b mismatch</th><th>MCIA</th>';

  const start = state.page * PAGE_SIZE;
  const page = state.matches.slice(start, start + PAGE_SIZE);
  const pct = (/** @type {number} */ v) => `${(v * 100).toFixed(2)}%`;

  body.innerHTML = page.map((m) => {
    const cells = [
      table.label('substrate', m.index),
      table.label('dimensions', m.index),
      pct(m.aMismatch),
      ...(dim === 1 ? [] : [pct(m.bMismatch)]),
      String(Math.round(m.mcia)),
    ];
    return `<tr>${cells.map((c) => `<td>${c}</td>`).join('')}</tr>`;
  }).join('') || '<tr><td colspan="5">No substrates selected</td></tr>';

  const total = state.matches.length;
  $('page-status').textContent = total
    ? `${start + 1}–${Math.min(start + PAGE_SIZE, total)} of top ${total}`
    : '—';
  /** @type {HTMLButtonElement} */ ($('page-prev')).disabled = state.page === 0;
  /** @type {HTMLButtonElement} */ ($('page-next')).disabled = start + PAGE_SIZE >= total;
}

function renderReadouts() {
  const n = state.substrates.size;
  const total = allSubstrateNames().length;
  const names = [...state.substrates];
  $('substrates-list').textContent = n === 0 ? 'None selected'
    : n === total ? 'All substrates selected'
    : n > 15 ? `${names.slice(0, 15).join(', ')} +${n - 15} more`
    : names.join(', ');

  $('selected-film').textContent = state.selectedFilm;
  const note = $('display-note');
  const messages = [];
  if (!webglAvailable()) {
    messages.push('WebGL is unavailable in this browser, so the plots are drawn as SVG. ' +
                  'They will be slower with many points but are otherwise identical.');
  }
  if (state.mode !== 'Scatter') {
    messages.push(getSolver().supported
      ? `Map solved on the GPU over ${visible['2d'].length.toLocaleString()} superlattices.`
      : 'Map computed on the CPU at reduced resolution (WebGL2 unavailable).');
  }
  note.textContent = messages.join(' ');
}

/* ------------------------------------------------------------ interaction */

function attachClick2d(gd) {
  if (gd._latticeClick) return;
  gd._latticeClick = true;
  gd.on('plotly_click', (ev) => {
    const pt = ev.points && ev.points[0];
    // Films are the only trace carrying customdata; substrate markers are not
    // selectable. The Dash app tested curveNumber === 1, which broke whenever
    // no films were loaded and the trace did not exist.
    if (!pt || pt.data.customdata == null) return;
    const row = pt.data.customdata[pt.pointIndex];
    const films = data.tables.films_2d;
    state.a = films.col('a')[row];
    state.b = films.col('b')[row];
    state.dimension = 2;
    state.selectedFilm =
      `${films.label('formula', row)} ${films.label('plane', row)} ${films.label('crystal_system', row)}`;
    syncInputs();
    refresh();
  });
}

function attachClick1d(gd) {
  if (gd._latticeClick) return;
  gd._latticeClick = true;
  gd.on('plotly_click', (ev) => {
    const pt = ev.points && ev.points[0];
    if (!pt || pt.data.customdata == null) return;
    const row = pt.data.customdata[pt.pointIndex];
    const films = data.tables.films_1d;
    state.a = films.col('a')[row];
    state.dimension = 1;
    state.selectedFilm =
      `${films.label('formula', row)} ${films.label('plane', row)} ${films.label('crystal_system', row)}`;
    syncInputs();
    refresh();
  });
}

function syncInputs() {
  /** @type {HTMLInputElement} */ ($('a-input')).value = String(Number(state.a.toFixed(4)));
  /** @type {HTMLInputElement} */ ($('b-input')).value = String(Number(state.b.toFixed(4)));
  $('b-field').hidden = state.dimension === 1;
  document.querySelectorAll('input[name="dimension"]').forEach((el) => {
    const input = /** @type {HTMLInputElement} */ (el);
    input.checked = Number(input.value) === state.dimension;
  });
}

function wireControls() {
  const on = (id, ev, fn) => $(id).addEventListener(ev, fn);

  on('single-double', 'input', (e) => {
    state.singleDouble = Number(/** @type {HTMLInputElement} */ (e.target).value);
    refresh();
  });
  on('large-super', 'input', (e) => {
    state.largeSuperlattice = Number(/** @type {HTMLInputElement} */ (e.target).value);
    refresh();
  });
  on('a-input', 'input', (e) => {
    const v = Number(/** @type {HTMLInputElement} */ (e.target).value);
    if (Number.isFinite(v) && v > 0) { state.a = v; state.selectedFilm = 'Manual Entry'; refresh(); }
  });
  on('b-input', 'input', (e) => {
    const v = Number(/** @type {HTMLInputElement} */ (e.target).value);
    if (Number.isFinite(v) && v > 0) { state.b = v; state.selectedFilm = 'Manual Entry'; refresh(); }
  });

  document.querySelectorAll('input[name="dimension"]').forEach((el) =>
    el.addEventListener('change', (e) => {
      state.dimension = Number(/** @type {HTMLInputElement} */ (e.target).value);
      $('b-field').hidden = state.dimension === 1;
      refresh();
    }));

  document.querySelectorAll('input[name="mode"]').forEach((el) =>
    el.addEventListener('change', (e) => {
      state.mode = /** @type {HTMLInputElement} */ (e.target).value;
      refresh();
    }));

  on('page-prev', 'click', () => { state.page--; renderTable(); });
  on('page-next', 'click', () => { state.page++; renderTable(); });
  on('download-csv', 'click', downloadCsv);

  // Modals
  const open = (id) => { $(id).hidden = false; };
  const close = (id) => { $(id).hidden = true; };
  on('learn-more', 'click', () => open('methodology-modal'));
  on('open-substrates', 'click', () => open('substrates-modal'));
  on('open-films', 'click', () => { open('films-modal'); ensureFilms(); });
  on('close-substrates', 'click', () => { close('substrates-modal'); refresh({ substrates: true }); });
  on('close-films', 'click', () => close('films-modal'));
  document.querySelectorAll('.modal-close').forEach((b) =>
    b.addEventListener('click', (e) => {
      const backdrop = /** @type {HTMLElement} */ (e.target).closest('.modal-backdrop');
      if (backdrop) /** @type {HTMLElement} */ (backdrop).hidden = true;
    }));
  document.querySelectorAll('.modal-backdrop').forEach((el) =>
    el.addEventListener('click', (e) => { if (e.target === el) /** @type {HTMLElement} */ (el).hidden = true; }));

  on('substrates-all', 'click', () => { selectAllSubstrates(); syncSubstrateBoxes(); });
  on('substrates-none', 'click', () => { state.substrates.clear(); syncSubstrateBoxes(); });
  on('substrate-search', 'input', (e) => {
    const q = /** @type {HTMLInputElement} */ (e.target).value.toLowerCase();
    $('substrate-list').querySelectorAll('label').forEach((l) => {
      /** @type {HTMLElement} */ (l).hidden = !l.textContent.toLowerCase().includes(q);
    });
  });

  on('pull-materials', 'click', pullMaterials);
}

/* ------------------------------------------------------------- substrates */

function allSubstrateNames() {
  const names = new Set();
  for (const dim of ['superlattices_2d', 'superlattices_1d']) {
    for (const n of data.tables[dim].categories('substrate')) names.add(n);
  }
  return [...names].sort();
}

function selectAllSubstrates() {
  state.substrates = new Set(allSubstrateNames());
}

function buildSubstrateList() {
  $('substrate-list').innerHTML = allSubstrateNames().map((n) =>
    `<label><input type="checkbox" value="${n}"> ${n}</label>`).join('');
  $('substrate-list').addEventListener('change', (e) => {
    const box = /** @type {HTMLInputElement} */ (e.target);
    if (box.checked) state.substrates.add(box.value);
    else state.substrates.delete(box.value);
  });
}

function syncSubstrateBoxes() {
  $('substrate-list').querySelectorAll('input').forEach((el) => {
    const box = /** @type {HTMLInputElement} */ (el);
    box.checked = state.substrates.has(box.value);
  });
}

/* ------------------------------------------------------------------ films */

async function ensureFilms() {
  if (data.tables.films_2d) return;
  $('film-status').textContent = 'Loading film catalogue…';
  await Promise.all([data.load('films_2d'), data.load('films_1d')]);
  $('film-status').textContent = 'Film catalogue ready. Set filters, then pull materials.';
}

function buildFilmFilters() {
  $('crystal-systems').innerHTML = Object.keys(PLANES_BY_SYSTEM).map((s) =>
    `<label><input type="checkbox" name="system" value="${s}"> ${s}</label>`).join('');
  $('plane-filters').innerHTML = Object.entries(PLANES_BY_SYSTEM).map(([sys, planes]) =>
    `<div class="plane-row"><strong>${sys}</strong>${planes.map((p) =>
      `<label><input type="checkbox" name="plane" data-system="${sys}" value="${p}" checked> ${p}</label>`
    ).join('')}</div>`).join('');
}

function parseElements(id) {
  return /** @type {HTMLInputElement} */ ($(id)).value
    .split(/[,\s]+/).map((s) => s.trim()).filter(Boolean);
}

async function pullMaterials() {
  await ensureFilms();
  const checked = (name) => [...document.querySelectorAll(`input[name="${name}"]:checked`)]
    .map((el) => /** @type {HTMLInputElement} */ (el).value);

  /** @type {Record<string, string[]>} */
  const planesBySystem = {};
  document.querySelectorAll('input[name="plane"]:checked').forEach((el) => {
    const input = /** @type {HTMLInputElement} */ (el);
    const sys = input.dataset.system || '';
    (planesBySystem[sys] = planesBySystem[sys] || []).push(input.value);
  });

  const numRaw = /** @type {HTMLInputElement} */ ($('num-elements')).value.trim();
  const filter = {
    mustInclude: parseElements('must-include'),
    canInclude: parseElements('can-include'),
    exclude: parseElements('exclude-elements'),
    numElements: numRaw === '' ? null : Number(numRaw),
    crystalSystems: checked('system'),
    pointGroups: parseElements('point-groups'),
    planesBySystem,
    maxRows: MAX_FILMS,
  };
  if (filter.numElements != null && !Number.isFinite(filter.numElements)) {
    $('film-status').textContent = 'Number of elements must be a whole number.';
    return;
  }

  const elements = data.manifest.elements;
  const r2 = filterFilms(data.tables.films_2d, elements, filter);
  const r1 = filterFilms(data.tables.films_1d, elements, filter);

  state.filmRows2d = [];
  state.filmRows1d = [];
  for (let i = 0; i < r2.mask.length; i++) if (r2.mask[i]) state.filmRows2d.push(i);
  for (let i = 0; i < r1.mask.length; i++) if (r1.mask[i]) state.filmRows1d.push(i);

  const total = state.filmRows2d.length + state.filmRows1d.length;
  const msg = total === 0 ? 'Search returned no films' : `${total.toLocaleString()} film planes selected`;
  $('film-status').textContent = msg;
  $('films-list').textContent = msg;
  refresh();
}

/* --------------------------------------------------------------- download */

function downloadCsv() {
  const dim = state.dimension;
  const table = superlattices(dim);
  const header = dim === 1
    ? ['substrate', 'dimensions', 'angle', 'a_mismatch', 'mcia', 'cost']
    : ['substrate', 'dimensions', 'angle', 'a_mismatch', 'b_mismatch', 'mcia', 'cost'];

  const rows = state.matches.map((m) => {
    const base = [
      table.label('substrate', m.index),
      table.label('dimensions', m.index),
      table.label('angle', m.index),
      m.aMismatch.toFixed(6),
    ];
    if (dim !== 1) base.push(m.bMismatch.toFixed(6));
    base.push(m.mcia.toFixed(3), m.cost.toFixed(6));
    return base;
  });

  const esc = (/** @type {string} */ v) => (/[",\n]/.test(v) ? `"${v.replace(/"/g, '""')}"` : v);
  const csv = [header, ...rows].map((r) => r.map(String).map(esc).join(',')).join('\n');

  const film = state.dimension === 1
    ? `a=${state.a.toFixed(4)}`
    : `a=${state.a.toFixed(4)}_b=${state.b.toFixed(4)}`;
  const blob = new Blob([csv], { type: 'text/csv' });
  const link = document.createElement('a');
  link.href = URL.createObjectURL(blob);
  link.download = `lattice-match_${film}.csv`;
  link.click();
  URL.revokeObjectURL(link.href);
}

main();
