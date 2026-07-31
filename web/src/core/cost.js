// @ts-check
/**
 * The matching cost functions.
 *
 * Direct port of `src/app.py`. This file and the Python must agree exactly --
 * `pipeline/tests/` and `web/test.html` run the same fixtures against both.
 * They drifted once before: the table normalised strain by the substrate while
 * the Voronoi normalised by the film and flipped its sign, so the map ranked
 * matches differently from the table beneath it.
 */

/**
 * Fractional strain of a film against a substrate.
 *
 * Substrate first. This is the convention used everywhere, matching the
 * "measure of strain" description in the app's own documentation.
 *
 * @param {number} sub substrate lattice parameter, in angstroms
 * @param {number} film film lattice parameter, in angstroms
 * @returns {number} signed fractional mismatch
 */
export function mismatch(sub, film) {
  return (film - sub) / sub;
}

/**
 * Cost of a square/rectangular match.
 *
 * An l-p quasi-norm over the two axis mismatches, times an MCIA penalty.
 * At p -> 0.05 the quasi-norm blows up unless one axis is nearly perfect, so
 * small p rewards a single very well matched axis; at p = 1 it is a plain sum,
 * rewarding both axes being decent.
 *
 * @param {number} aMismatch
 * @param {number} bMismatch
 * @param {number} mcia minimum coincident interface area, square angstroms
 * @param {number} singleDouble p, in [0.05, 1]
 * @param {number} largeSuperlattice q, in [0.05, 1]; at q = 1 the MCIA penalty vanishes
 */
export function cost2d(aMismatch, bMismatch, mcia, singleDouble, largeSuperlattice) {
  const p = singleDouble;
  return (
    Math.pow(Math.pow(Math.abs(aMismatch), p) + Math.pow(Math.abs(bMismatch), p), 1 / p) *
    Math.pow(mcia, 1 - largeSuperlattice)
  );
}

/**
 * Cost of a triangular match. Only one axis, so no p-norm.
 * @param {number} aMismatch
 * @param {number} mcia
 * @param {number} largeSuperlattice
 */
export function cost1d(aMismatch, mcia, largeSuperlattice) {
  return Math.abs(aMismatch) * Math.pow(mcia, 1 - largeSuperlattice);
}

/**
 * Marker area scaling: smaller MCIA means a better match, so a bigger dot.
 * Mirrors the sizing in the Dash figures.
 * @param {number} mcia
 * @param {number} largeSuperlattice
 */
export function markerSize(mcia, largeSuperlattice) {
  return 10 * Math.pow(100 / mcia, 1 - largeSuperlattice);
}
