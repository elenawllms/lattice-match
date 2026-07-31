// @ts-nocheck
/**
 * Differential test: the browser implementation against the Python one.
 *
 * The cost function now exists twice -- once in `src/app.py` / `pipeline/` and
 * once in `web/src/core/`. Two implementations of the same physics is exactly
 * how the old app ended up with a table and a Voronoi map that ranked matches
 * differently. This runs both over the same fixtures and the same shipped
 * binary bundle and requires them to agree.
 *
 * Driven by `pipeline/tests/test_web_differential.py`, which generates the
 * reference values. Runs under JavaScriptCore (`jsc -m`), which ships with
 * macOS, so it needs no Node installation.
 *
 * Usage: jsc -m differential.mjs -- <reference.json> <bundle.b64> <meta.json> <coreDir>
 */

const [refPath, binPath, metaPath, coreDir] = arguments;

const { mismatch, cost1d, cost2d } = await import(`${coreDir}/cost.js`);
const { best2d, best1d } = await import(`${coreDir}/match.js`);

const ref = JSON.parse(readFile(refPath));
const meta = JSON.parse(readFile(metaPath));

let pass = 0;
const failures = [];

const near = (a, b, tol) => Math.abs(a - b) <= tol * Math.max(1, Math.abs(b));

function check(name, got, want, tol = 1e-12) {
  if (near(got, want, tol)) pass++;
  else failures.push(`${name}: js=${got} py=${want}`);
}

// --- scalar functions ------------------------------------------------------
for (const [s, f, want] of ref.mismatch) check(`mismatch(${s},${f})`, mismatch(s, f), want);
for (const c of ref.cost2d) check(`cost2d(${c.slice(0, 5)})`, cost2d(c[0], c[1], c[2], c[3], c[4]), c[5]);
for (const c of ref.cost1d) check(`cost1d(${c.slice(0, 3)})`, cost1d(c[0], c[1], c[2]), c[3]);

// --- ranking over the real catalogue, read from the shipped bundle ---------
const bin = Uint8Array.from(atob(readFile(binPath).trim()), (ch) => ch.charCodeAt(0));
const VIEWS = {
  float32: Float32Array, float64: Float64Array,
  uint8: Uint8Array, uint16: Uint16Array, uint32: Uint32Array,
};
const table = {
  rows: meta.rows,
  col: (name) => {
    const c = meta.columns[name];
    return new VIEWS[c.dtype](bin.buffer, c.offset, c.length);
  },
};

const r = ref.ranking;
const got = best2d(table, r.a, r.b, r.p, r.q, r.top.length);
if (got.length !== r.top.length) {
  failures.push(`ranking length: js=${got.length} py=${r.top.length}`);
} else {
  for (let k = 0; k < r.top.length; k++) {
    if (got[k].index !== r.top[k][0]) {
      failures.push(`rank ${k} index: js=${got[k].index} py=${r.top[k][0]}`);
    } else pass++;
    check(`rank ${k} cost`, got[k].cost, r.top[k][1]);
  }
}

// --- 1d ranking ------------------------------------------------------------
if (ref.ranking1d) {
  const meta1 = JSON.parse(readFile(ref.ranking1d.meta));
  const bin1 = Uint8Array.from(atob(readFile(ref.ranking1d.bin).trim()), (ch) => ch.charCodeAt(0));
  const t1 = {
    rows: meta1.rows,
    col: (name) => {
      const c = meta1.columns[name];
      return new VIEWS[c.dtype](bin1.buffer, c.offset, c.length);
    },
  };
  const g1 = best1d(t1, ref.ranking1d.a, ref.ranking1d.q, ref.ranking1d.top.length);
  for (let k = 0; k < ref.ranking1d.top.length; k++) {
    if (g1[k].index !== ref.ranking1d.top[k][0]) {
      failures.push(`1d rank ${k} index: js=${g1[k].index} py=${ref.ranking1d.top[k][0]}`);
    } else pass++;
    check(`1d rank ${k} cost`, g1[k].cost, ref.ranking1d.top[k][1]);
  }
}

print(JSON.stringify({ pass, failures }));
