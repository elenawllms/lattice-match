// @ts-check
/**
 * Plotly figure construction.
 *
 * Two things differ from the Dash version beyond running client-side:
 *
 *  - `scattergl` rather than `scatter`. The old figures drew ~8,200 SVG path
 *    nodes on first paint, which is what actually froze browsers.
 *  - The Voronoi and heatmap arrive as a `layout.images` entry backed by an
 *    offscreen canvas, instead of a `go.Image` trace carrying a 100x100x3
 *    array as JSON. That is why the gl2d bundle suffices: no image or heatmap
 *    trace module is needed.
 */

/* global Plotly */

// Recessive chrome: hairline grid, muted tick labels, no heavy frame. Values
// are the data-viz reference tokens, so the plots and the page agree.
const INK = { primary: '#0b0b0b', muted: '#898781', grid: '#e1e0d9', axis: '#c3c2b7' };

const AXIS = {
  mirror: false,
  ticks: 'outside',
  ticklen: 4,
  tickcolor: INK.axis,
  showline: true,
  linecolor: INK.axis,
  linewidth: 1,
  gridcolor: INK.grid,
  gridwidth: 1,
  zeroline: false,
  tickfont: { size: 11, color: INK.muted },
  title: { font: { size: 12, color: INK.muted } },
};

const FONT = { family: 'system-ui, -apple-system, "Segoe UI", Roboto, sans-serif' };

export const VIEW_2D = { aMin: 3, aMax: 15, bMin: 3, bMax: 15 };
const INITIAL_RANGE = [4, 10];

const CONFIG = {
  displaylogo: false,
  responsive: true,
  toImageButtonOptions: { filename: 'lattice-match', scale: 2 },
  modeBarButtonsToRemove: ['select2d', 'lasso2d'],
};

/**
 * @param {import('../data/bundle.js').Table} table
 * @param {number[]} rows row indices to draw
 * @param {number} largeSuperlattice
 * @param {boolean} sizeByMcia
 * @param {boolean} webgl false falls back to SVG, which is slow past a few
 *   thousand points but is the only option without WebGL
 */
/**
 * Positions, colours and hover strings for a substrate selection.
 *
 * These depend only on *which* rows are drawn, so they are built once per
 * selection and reused. Only marker size depends on the sliders. Rebuilding all
 * of it on every slider event meant ~7,300 dictionary lookups and string
 * concatenations per frame, which is what made dragging feel heavy.
 *
 * @type {WeakMap<object, Map<string, any>>}
 */
const staticTraceCache = new WeakMap();

function staticTraceData(table, rows, cacheKey) {
  let perTable = staticTraceCache.get(table);
  if (!perTable) staticTraceCache.set(table, (perTable = new Map()));
  const hit = perTable.get(cacheKey);
  if (hit) return hit;

  const a = table.col('a');
  const b = table.columns.b;
  const mcia = table.col('mcia');
  const red = table.col('red');
  const green = table.col('green');
  const blue = table.col('blue');

  const x = new Array(rows.length);
  const y = new Array(rows.length);
  const color = new Array(rows.length);
  const text = new Array(rows.length);

  for (let k = 0; k < rows.length; k++) {
    const i = rows[k];
    x[k] = a[i];
    y[k] = b ? b[i] : 0;
    color[k] = `rgb(${red[i]},${green[i]},${blue[i]})`;
    text[k] = `${table.label('substrate', i)} ${table.label('dimensions', i)}` +
              `<br>MCIA: ${Math.round(mcia[i])} sq. Å` +
              `<br>angle: ${table.label('angle', i)}`;
  }

  const built = { x, y, color, text, mcia };
  // One entry per table: selections change, and holding every past selection
  // would leak. The WeakMap keys on the table so unloaded data can be freed.
  perTable.clear();
  perTable.set(cacheKey, built);
  return built;
}

function substrateTrace(table, rows, largeSuperlattice, sizeByMcia, webgl, cacheKey) {
  const { x, y, color, text, mcia } = staticTraceData(table, rows, cacheKey);

  const size = new Array(rows.length);
  if (sizeByMcia) {
    const exponent = 1 - largeSuperlattice;
    for (let k = 0; k < rows.length; k++) size[k] = 10 * Math.pow(100 / mcia[rows[k]], exponent);
  } else {
    size.fill(7);
  }

  return {
    type: webgl ? 'scattergl' : 'scatter',
    mode: 'markers',
    x,
    y,
    // A 1px surface ring separates overlapping marks without the heavy black
    // outline the original drew on all 7,343 of them.
    marker: { color, size, symbol: 'square',
              line: { width: 1, color: 'rgba(252,252,251,0.9)' } },
    text,
    hovertemplate: '%{text}<br>a: %{x:.4f} Å<br>b: %{y:.4f} Å<extra></extra>',
    showlegend: false,
    name: '',
  };
}

/**
 * @param {import('../data/bundle.js').Table} films
 * @param {number[]} rows
 * @param {boolean} twoD
 * @param {boolean} webgl
 */
function filmTrace(films, rows, twoD, webgl) {
  const a = films.col('a');
  const b = twoD ? films.col('b') : null;
  return {
    type: webgl ? 'scattergl' : 'scatter',
    mode: 'markers',
    x: rows.map((i) => a[i]),
    y: rows.map((i) => (b ? b[i] : 0)),
    marker: { color: INK.primary, size: 6, symbol: 'circle',
              line: { width: 1, color: 'rgba(252,252,251,0.9)' } },
    text: rows.map((i) =>
      `${films.label('formula', i)} ${films.label('plane', i)} ${films.label('crystal_system', i)}`),
    hovertemplate: twoD
      ? '%{text}<br>a: %{x:.4f} Å<br>b: %{y:.4f} Å<extra></extra>'
      : '%{text}<br>a: %{x:.4f} Å<extra></extra>',
    showlegend: false,
    name: '',
    customdata: rows,
  };
}

/**
 * Draw the square/rectangular plot.
 *
 * @param {HTMLElement} el
 * @param {object} opts
 */
export function render2d(el, opts) {
  const { superlattices, rows, films, filmRows, largeSuperlattice, mode, raster } = opts;
  const webgl = opts.webgl !== false;
  const traces = [substrateTrace(superlattices, rows, largeSuperlattice,
                                 mode === 'Scatter', webgl, opts.cacheKey)];
  if (films && filmRows.length) traces.push(filmTrace(films, filmRows, true, webgl));

  /** @type {any} */
  const layout = {
    xaxis: { ...AXIS, title: { text: 'a (Å)' }, range: INITIAL_RANGE, rangemode: 'nonnegative' },
    yaxis: { ...AXIS, title: { text: 'b (Å)' }, range: INITIAL_RANGE, rangemode: 'nonnegative' },
    plot_bgcolor: '#fcfcfb',
    paper_bgcolor: '#fcfcfb',
    font: FONT,
    margin: { l: 56, r: 12, t: 8, b: 48 },
    hovermode: 'closest',
    hoverlabel: { bgcolor: '#ffffff', bordercolor: INK.grid,
                  font: { size: 12, color: INK.primary, ...FONT } },
    uirevision: 'keep',
  };

  if (raster) {
    // Anchored to data coordinates, so it pans and zooms with the scatter.
    layout.images = [{
      source: raster.toDataURL('image/png'),
      xref: 'x', yref: 'y',
      x: VIEW_2D.aMin, y: VIEW_2D.bMax,
      sizex: VIEW_2D.aMax - VIEW_2D.aMin,
      sizey: VIEW_2D.bMax - VIEW_2D.bMin,
      sizing: 'stretch',
      opacity: mode === 'Voronoi' ? 0.45 : 0.75,
      layer: 'below',
    }];
  }

  return Plotly.react(el, traces, layout, CONFIG);
}

/**
 * Draw the triangular strip.
 * @param {HTMLElement} el
 * @param {object} opts
 */
export function render1d(el, opts) {
  const { superlattices, rows, films, filmRows, largeSuperlattice, mode } = opts;
  const webgl = opts.webgl !== false;
  const traces = [substrateTrace(superlattices, rows, largeSuperlattice,
                                 mode === 'Scatter', webgl, opts.cacheKey)];
  if (films && filmRows.length) traces.push(filmTrace(films, filmRows, false, webgl));

  const layout = {
    height: 190,
    xaxis: { ...AXIS, title: { text: 'a (Å)' }, range: INITIAL_RANGE },
    yaxis: {
      ...AXIS, range: [-1, 1], fixedrange: true,
      zeroline: true, zerolinewidth: 1, zerolinecolor: INK.axis, showticklabels: false,
    },
    plot_bgcolor: '#fcfcfb',
    paper_bgcolor: '#fcfcfb',
    font: FONT,
    margin: { l: 56, r: 12, t: 8, b: 44 },
    hovermode: 'closest',
    hoverlabel: { bgcolor: '#ffffff', bordercolor: INK.grid,
                  font: { size: 12, color: INK.primary, ...FONT } },
    uirevision: 'keep',
  };
  return Plotly.react(el, traces, layout, CONFIG);
}
