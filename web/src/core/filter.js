// @ts-check
/**
 * Film catalogue filtering.
 *
 * Element membership is a 128-bit mask packed as 4 uint32s per row, so the
 * must/can/exclude tests are a handful of bitwise ops. The Dash app ran six
 * row-wise Python lambdas over 81,000 rows for the same thing, each allocating
 * an intermediate list.
 */

/** uint32 words per element mask; must match `MASK_WORDS` in pipeline/encode.py. */
export const MASK_WORDS = 4;

/**
 * Build a 4-word mask from element symbols.
 * @param {string[]} symbols
 * @param {string[]} elements the full element list from the manifest
 * @returns {Uint32Array}
 */
export function maskOf(symbols, elements) {
  const mask = new Uint32Array(MASK_WORDS);
  for (const sym of symbols) {
    const idx = elements.indexOf(sym);
    if (idx < 0) continue;
    mask[idx >>> 5] |= 1 << (idx & 31);
  }
  return mask;
}

/**
 * @typedef {object} FilmFilter
 * @property {string[]} [mustInclude] every one of these elements must be present
 * @property {string[]} [canInclude] no element outside this set may be present
 * @property {string[]} [exclude] none of these may be present
 * @property {number|null} [numElements]
 * @property {string[]} [crystalSystems]
 * @property {string[]} [pointGroups]
 * @property {Record<string, string[]>} [planesBySystem] allowed planes per crystal system
 * @property {number} [maxRows] stop after this many matches
 */

/**
 * Apply a filter, returning a per-row 0/1 mask.
 *
 * Empty and absent criteria both mean "no constraint". The Dash app tested
 * `is None` here, but a cleared multi-select yields `[]`, so `all(e in [] ...)`
 * was False for every row and searches silently returned nothing.
 *
 * @param {import('../data/bundle.js').Table} films
 * @param {string[]} elements
 * @param {FilmFilter} filter
 * @returns {{mask: Uint8Array, count: number}}
 */
export function filterFilms(films, elements, filter) {
  const rows = films.rows;
  const mask = new Uint8Array(rows);

  const must = maskOf(filter.mustInclude || [], elements);
  const can = maskOf(filter.canInclude || [], elements);
  const excl = maskOf(filter.exclude || [], elements);
  const hasMust = (filter.mustInclude || []).length > 0;
  const hasCan = (filter.canInclude || []).length > 0;
  const hasExcl = (filter.exclude || []).length > 0;

  const em = films.col('element_mask');
  const numElements = films.col('num_elements');
  const systemCodes = films.col('crystal_system');
  const pointGroupCodes = films.col('point_group');
  const planeCodes = films.col('plane');

  const systems = films.categories('crystal_system');
  const pointGroups = films.categories('point_group');
  const planes = films.categories('plane');

  // Precompute per-category allow tables so the row loop only does lookups.
  const systemAllowed = allowTable(systems, filter.crystalSystems);
  const pointGroupAllowed = allowTable(pointGroups, filter.pointGroups);

  // planesBySystem is keyed by crystal system, so build a (system, plane) grid.
  /** @type {Uint8Array|null} */
  let planeGrid = null;
  if (filter.planesBySystem) {
    planeGrid = new Uint8Array(systems.length * planes.length);
    for (let s = 0; s < systems.length; s++) {
      const allowed = filter.planesBySystem[systems[s]];
      for (let p = 0; p < planes.length; p++) {
        // A system with no entry is unconstrained, matching the UI, where an
        // untouched checklist means "all planes".
        planeGrid[s * planes.length + p] = !allowed || allowed.includes(planes[p]) ? 1 : 0;
      }
    }
  }

  const wantNum = filter.numElements == null ? -1 : filter.numElements;
  const limit = filter.maxRows == null ? Infinity : filter.maxRows;
  let count = 0;

  for (let i = 0; i < rows; i++) {
    if (wantNum >= 0 && numElements[i] !== wantNum) continue;
    if (systemAllowed && !systemAllowed[systemCodes[i]]) continue;
    if (pointGroupAllowed && !pointGroupAllowed[pointGroupCodes[i]]) continue;
    if (planeGrid && !planeGrid[systemCodes[i] * planes.length + planeCodes[i]]) continue;

    if (hasMust || hasCan || hasExcl) {
      const base = i * MASK_WORDS;
      let ok = true;
      for (let w = 0; w < MASK_WORDS; w++) {
        const row = em[base + w];
        if (hasMust && (row & must[w]) !== must[w]) { ok = false; break; }
        if (hasExcl && (row & excl[w]) !== 0) { ok = false; break; }
        // "can include" means the row introduces no element outside the set.
        if (hasCan && (row & ~can[w]) !== 0) { ok = false; break; }
      }
      if (!ok) continue;
    }

    mask[i] = 1;
    if (++count >= limit) break;
  }
  return { mask, count };
}

/**
 * @param {string[]} categories
 * @param {string[]|undefined} allowed
 * @returns {Uint8Array|null} null when unconstrained
 */
function allowTable(categories, allowed) {
  if (!allowed || allowed.length === 0) return null;
  const table = new Uint8Array(categories.length);
  for (let i = 0; i < categories.length; i++) {
    table[i] = allowed.includes(categories[i]) ? 1 : 0;
  }
  return table;
}
