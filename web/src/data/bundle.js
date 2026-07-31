// @ts-check
/**
 * Decoder for the binary bundles written by `pipeline/encode.py`.
 *
 * Each table is one flat ArrayBuffer of concatenated typed arrays plus a
 * sidecar JSON describing every column's dtype, byte offset and length. Columns
 * become zero-copy views over that buffer -- no per-row objects are ever
 * allocated, which is the whole point: the Dash app spent most of its time
 * building and serialising dictionaries.
 */

/** @type {Record<string, new (b: ArrayBuffer, o: number, n: number) => any>} */
const VIEWS = {
  float32: Float32Array,
  float64: Float64Array,
  uint8: Uint8Array,
  uint16: Uint16Array,
  uint32: Uint32Array,
  int8: Int8Array,
};

/**
 * @typedef {object} Column
 * @property {string} dtype
 * @property {number} offset
 * @property {number} length
 */

/**
 * @typedef {object} TableMeta
 * @property {number} rows
 * @property {Record<string, Column>} columns
 * @property {Record<string, string[]>} dictionaries
 * @property {string} file
 */

/** A decoded table: typed-array columns plus string dictionaries. */
export class Table {
  /**
   * @param {TableMeta} meta
   * @param {ArrayBuffer} buffer
   */
  constructor(meta, buffer) {
    this.rows = meta.rows;
    this.dictionaries = meta.dictionaries || {};
    /** @type {Record<string, any>} */
    this.columns = {};

    for (const [name, col] of Object.entries(meta.columns)) {
      const View = VIEWS[col.dtype];
      if (!View) throw new Error(`bundle: unknown dtype ${col.dtype} for column ${name}`);
      // Throws if the encoder ever emits a misaligned offset, rather than
      // silently decoding neighbouring bytes as values.
      this.columns[name] = new View(buffer, col.offset, col.length);
    }
  }

  /** @param {string} name */
  col(name) {
    const c = this.columns[name];
    if (!c) throw new Error(`bundle: no column ${name}`);
    return c;
  }

  /** Resolve a dictionary-encoded column to its string value for one row. */
  /** @param {string} name @param {number} row */
  label(name, row) {
    return this.dictionaries[name][this.col(name)[row]];
  }

  /** The full string table for a dictionary column. */
  /** @param {string} name */
  categories(name) {
    return this.dictionaries[name];
  }
}

/**
 * Fetch and decode one table.
 * @param {string} base directory the bundles are served from
 * @param {{meta: string, file: string}} entry
 * @returns {Promise<Table>}
 */
export async function loadTable(base, entry) {
  const [meta, buffer] = await Promise.all([
    fetch(`${base}/${entry.meta}`).then((r) => {
      if (!r.ok) throw new Error(`bundle: ${entry.meta} -> HTTP ${r.status}`);
      return r.json();
    }),
    fetch(`${base}/${entry.file}`).then((r) => {
      if (!r.ok) throw new Error(`bundle: ${entry.file} -> HTTP ${r.status}`);
      return r.arrayBuffer();
    }),
  ]);
  return new Table(meta, buffer);
}

/**
 * Load the manifest and every table marked eager. Films are left alone until
 * something asks for them -- they are ~730 KB gzipped against 24 KB for the
 * superlattices, and the opening view does not use them.
 * @param {string} base
 */
export async function loadEager(base) {
  const manifest = await fetch(`${base}/manifest.json`).then((r) => {
    if (!r.ok) throw new Error(`bundle: manifest.json -> HTTP ${r.status}`);
    return r.json();
  });

  /** @type {Record<string, Table>} */
  const tables = {};
  await Promise.all(
    Object.entries(manifest.tables)
      .filter(([, t]) => /** @type {any} */ (t).eager)
      .map(async ([name, t]) => {
        tables[name] = await loadTable(base, /** @type {any} */ (t));
      })
  );

  const pending = new Map();
  /**
   * Fetch a lazy table, memoised so concurrent callers share one request.
   * @param {string} name
   * @returns {Promise<Table>}
   */
  const load = (name) => {
    if (tables[name]) return Promise.resolve(tables[name]);
    if (!pending.has(name)) {
      const entry = manifest.tables[name];
      if (!entry) throw new Error(`bundle: no table ${name}`);
      pending.set(
        name,
        loadTable(base, entry).then((t) => {
          tables[name] = t;
          return t;
        })
      );
    }
    return pending.get(name);
  };

  return { manifest, tables, load };
}
