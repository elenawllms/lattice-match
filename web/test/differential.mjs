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
 * reference values.
 *
 * Runs under either runtime, so it works locally on macOS via JavaScriptCore
 * (which ships with the OS, needing no Node install) and in CI on Linux via
 * Node. Skipping in CI would defeat the point of having it.
 *
 *   jsc  -m differential.mjs -- <reference.json> <bundle.b64> <meta.json> <coreDir>
 *   node    differential.mjs    <reference.json> <bundle.b64> <meta.json> <coreDir>
 */

const isNode = typeof process !== 'undefined' && process.versions && process.versions.node;

/** @type {(p: string) => string} */
let readText;
/** @type {(s: string) => void} */
let emit;
/** @type {string[]} */
let argv;
/** @type {(dir: string, file: string) => string} */
let moduleSpecifier;

if (isNode) {
  const { readFileSync } = await import('node:fs');
  const { pathToFileURL } = await import('node:url');
  readText = (p) => readFileSync(p, 'utf8');
  emit = (s) => console.log(s);
  argv = process.argv.slice(2);
  moduleSpecifier = (dir, file) => pathToFileURL(`${dir}/${file}`).href;
} else {
  // JavaScriptCore: `readFile` and `print` are globals, and arguments after
  // `--` arrive in the `arguments` global.
  readText = (p) => readFile(p);
  emit = (s) => print(s);
  argv = Array.from(arguments);
  moduleSpecifier = (dir, file) => `${dir}/${file}`;
}

const [refPath, binPath, metaPath, coreDir] = argv;
if (!refPath || !binPath || !metaPath || !coreDir) {
  emit(JSON.stringify({ pass: 0, failures: [`bad arguments: ${JSON.stringify(argv)}`] }));
  throw new Error('differential: missing arguments');
}

const { mismatch, cost1d, cost2d } = await import(moduleSpecifier(coreDir, 'cost.js'));
const { best2d, best1d } = await import(moduleSpecifier(coreDir, 'match.js'));

const ref = JSON.parse(readText(refPath));
const meta = JSON.parse(readText(metaPath));

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
const bin = Uint8Array.from(atob(readText(binPath).trim()), (ch) => ch.charCodeAt(0));
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
  const meta1 = JSON.parse(readText(ref.ranking1d.meta));
  const bin1 = Uint8Array.from(atob(readText(ref.ranking1d.bin).trim()), (ch) => ch.charCodeAt(0));
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

emit(JSON.stringify({ pass, failures }));
