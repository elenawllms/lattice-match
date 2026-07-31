// @ts-check
/**
 * Ranking substrate superlattices against a film.
 *
 * One pass over typed arrays with a bounded insertion, rather than scoring
 * every row into an object and sorting all of them. For 7,343 superlattices
 * this is well under a millisecond, which is what makes the sliders feel live
 * -- the Dash app round-tripped 1.7 MB to the server for the same answer.
 */

import { cost1d, cost2d, mismatch } from './cost.js';

/**
 * @typedef {object} Match
 * @property {number} index row in the superlattice table
 * @property {number} cost
 * @property {number} aMismatch
 * @property {number} [bMismatch] absent for triangular matches
 * @property {number} mcia
 */

/**
 * Best `limit` matches for a rectangular film, ranked by cost.
 *
 * @param {{col: (n: string) => any, rows: number}} table
 * @param {number} filmA
 * @param {number} filmB
 * @param {number} singleDouble
 * @param {number} largeSuperlattice
 * @param {number} limit
 * @param {Uint8Array} [mask] optional per-row filter; rows with 0 are skipped
 * @returns {Match[]}
 */
export function best2d(table, filmA, filmB, singleDouble, largeSuperlattice, limit, mask) {
  const a = table.col('a');
  const b = table.col('b');
  const mcia = table.col('mcia');

  /** @type {Match[]} */
  const heap = [];
  let worst = Infinity;

  for (let i = 0; i < table.rows; i++) {
    if (mask && !mask[i]) continue;
    const am = mismatch(a[i], filmA);
    const bm = mismatch(b[i], filmB);
    const c = cost2d(am, bm, mcia[i], singleDouble, largeSuperlattice);
    // Cheap reject before touching the array at all.
    if (heap.length === limit && c >= worst) continue;

    const entry = { index: i, cost: c, aMismatch: am, bMismatch: bm, mcia: mcia[i] };
    insert(heap, entry, limit);
    worst = heap[heap.length - 1].cost;
  }
  return heap;
}

/**
 * Best `limit` matches for a triangular film.
 *
 * @param {{col: (n: string) => any, rows: number}} table
 * @param {number} filmA
 * @param {number} largeSuperlattice
 * @param {number} limit
 * @param {Uint8Array} [mask]
 * @returns {Match[]}
 */
export function best1d(table, filmA, largeSuperlattice, limit, mask) {
  const a = table.col('a');
  const mcia = table.col('mcia');

  /** @type {Match[]} */
  const heap = [];
  let worst = Infinity;

  for (let i = 0; i < table.rows; i++) {
    if (mask && !mask[i]) continue;
    const am = mismatch(a[i], filmA);
    const c = cost1d(am, mcia[i], largeSuperlattice);
    if (heap.length === limit && c >= worst) continue;

    insert(heap, { index: i, cost: c, aMismatch: am, mcia: mcia[i] }, limit);
    worst = heap[heap.length - 1].cost;
  }
  return heap;
}

/**
 * Insertion into a short sorted array, capped at `limit`.
 *
 * A binary search beats a real heap here: `limit` is ~200, and the early-reject
 * above means very few candidates ever reach this function.
 *
 * The comparison is `<=`, not `<`, so an entry lands *after* any existing
 * entries of equal cost. That makes the sort stable, breaking ties by table
 * order and matching Python's `sorted`. Exact ties are common and not
 * pathological: two substrates with the same lattice constant -- InSb and CdTe
 * are both cubic at 6.48 A -- produce numerically identical superlattices.
 * Without this, the two implementations disagree on their order and the
 * displayed ranking is not reproducible.
 *
 * @param {Match[]} sorted
 * @param {Match} entry
 * @param {number} limit
 */
function insert(sorted, entry, limit) {
  let lo = 0;
  let hi = sorted.length;
  while (lo < hi) {
    const mid = (lo + hi) >> 1;
    if (sorted[mid].cost <= entry.cost) lo = mid + 1;
    else hi = mid;
  }
  sorted.splice(lo, 0, entry);
  if (sorted.length > limit) sorted.pop();
}
