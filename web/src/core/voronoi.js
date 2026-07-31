// @ts-check
/**
 * The many-to-many map: which substrate wins at every point in (a, b) space,
 * and how good that best match is.
 *
 * This was the feature that made the old site unusable. It ran `np.vectorize`
 * -- a Python loop over 10,000 grid points, each doing ~8 pandas Series ops and
 * evaluating the cost function twice -- then shipped four 100x100 float grids
 * through the browser and back. It was capped at 1,000 superlattices and still
 * took seconds, and it ran even on the Scatter tab where the result was thrown
 * away.
 *
 * Here every grid cell is a fragment. The superlattice catalogue is uploaded
 * once as a float texture and the shader loops over it, so the whole solve is
 * one draw call at full canvas resolution. The 1,000-superlattice cap is gone.
 *
 * A CPU fallback with identical output covers machines without WebGL2.
 */

import { cost2d, mismatch } from './cost.js';

const VERT = `#version 300 es
in vec2 position;
void main() { gl_Position = vec4(position, 0.0, 1.0); }`;

/**
 * For each fragment, scan the catalogue and keep the lowest cost.
 *
 * The substrate columns live in one RGBA32F texture: r = a, g = b, b = mcia,
 * and the colour is fetched from a second texture at the winning index. Looping
 * in the shader rather than blending across draws keeps it to a single pass.
 */
const FRAG = `#version 300 es
precision highp float;
precision highp int;

uniform sampler2D uParams;   // rgba32f: a, b, mcia, unused
uniform sampler2D uColors;   // rgba8:   substrate colour
uniform int   uCount;        // superlattices in the catalogue
uniform int   uTexWidth;     // texture width; rows wrap at this
uniform vec4  uLimits;       // aMin, aMax, bMin, bMax
uniform vec2  uResolution;
uniform float uSingleDouble; // p
uniform float uLargeSuper;   // q
uniform int   uMode;         // 0 = winner colour (Voronoi), 1 = cost heatmap
uniform float uCostScale;    // normalisation for the heatmap ramp

out vec4 fragColor;

vec3 temps(float t) {
  // Approximation of Plotly's "Temps" ramp, so the heatmap and its colourbar
  // stay visually consistent with the rest of the app.
  vec3 c0 = vec3(0.0090, 0.3400, 0.4160);
  vec3 c1 = vec3(0.4000, 0.7000, 0.6600);
  vec3 c2 = vec3(0.9700, 0.8600, 0.6300);
  vec3 c3 = vec3(0.8900, 0.4700, 0.3200);
  vec3 c4 = vec3(0.6600, 0.1400, 0.1900);
  t = clamp(t, 0.0, 1.0);
  if (t < 0.25) return mix(c0, c1, t / 0.25);
  if (t < 0.50) return mix(c1, c2, (t - 0.25) / 0.25);
  if (t < 0.75) return mix(c2, c3, (t - 0.50) / 0.25);
  return mix(c3, c4, (t - 0.75) / 0.25);
}

void main() {
  vec2 uv = gl_FragCoord.xy / uResolution;
  float filmA = mix(uLimits.x, uLimits.y, uv.x);
  float filmB = mix(uLimits.z, uLimits.w, uv.y);

  float bestCost = 1e30;
  int   bestIdx  = 0;

  for (int i = 0; i < 65536; i++) {
    if (i >= uCount) break;
    ivec2 p = ivec2(i % uTexWidth, i / uTexWidth);
    vec4 s = texelFetch(uParams, p, 0);

    // Substrate first, film second -- the same convention as cost.js.
    float am = (filmA - s.r) / s.r;
    float bm = (filmB - s.g) / s.g;
    float c = pow(pow(abs(am), uSingleDouble) + pow(abs(bm), uSingleDouble),
                  1.0 / uSingleDouble) * pow(s.b, 1.0 - uLargeSuper);
    if (c < bestCost) { bestCost = c; bestIdx = i; }
  }

  if (uMode == 1) {
    fragColor = vec4(temps(bestCost / uCostScale), 1.0);
  } else {
    ivec2 p = ivec2(bestIdx % uTexWidth, bestIdx / uTexWidth);
    fragColor = vec4(texelFetch(uColors, p, 0).rgb, 1.0);
  }
}`;

/** Solves the grid on the GPU, rendering to an offscreen canvas. */
export class VoronoiSolver {
  /** @param {number} size square canvas edge, in pixels */
  constructor(size = 512) {
    this.size = size;
    this.canvas = document.createElement('canvas');
    this.canvas.width = size;
    this.canvas.height = size;
    /** @type {WebGL2RenderingContext|null} */
    this.gl = /** @type {any} */ (
      this.canvas.getContext('webgl2', { preserveDrawingBuffer: true, antialias: false })
    );
    this.ready = false;
    this.count = 0;
    if (this.gl) this._init();
  }

  get supported() {
    return this.gl !== null && this.ready;
  }

  _init() {
    const gl = /** @type {WebGL2RenderingContext} */ (this.gl);
    // RGBA32F textures are an optional WebGL2 feature; without it, fall back.
    if (!gl.getExtension('EXT_color_buffer_float') && !gl.getExtension('OES_texture_float_linear')) {
      // texelFetch on a float texture still works without these in practice,
      // so only bail if the program itself fails to link below.
    }
    const program = linkProgram(gl, VERT, FRAG);
    if (!program) return;
    this.program = program;

    const quad = gl.createBuffer();
    gl.bindBuffer(gl.ARRAY_BUFFER, quad);
    gl.bufferData(gl.ARRAY_BUFFER, new Float32Array([-1, -1, 3, -1, -1, 3]), gl.STATIC_DRAW);
    const loc = gl.getAttribLocation(program, 'position');
    gl.enableVertexAttribArray(loc);
    gl.vertexAttribPointer(loc, 2, gl.FLOAT, false, 0, 0);

    this.paramTex = gl.createTexture();
    this.colorTex = gl.createTexture();
    this.ready = true;
  }

  /**
   * Upload the superlattice catalogue. Only needed when the selection changes,
   * not when a slider moves.
   * @param {Float32Array} a
   * @param {Float32Array} b
   * @param {Float32Array} mcia
   * @param {Uint8Array} rgb interleaved r,g,b per row
   */
  upload(a, b, mcia, rgb) {
    if (!this.supported) return;
    const gl = /** @type {WebGL2RenderingContext} */ (this.gl);
    const n = a.length;
    // Square-ish texture; sampler dimensions are capped, a flat row is not.
    const width = Math.min(2048, Math.max(1, Math.ceil(Math.sqrt(n))));
    const height = Math.ceil(n / width);
    this.count = n;
    this.texWidth = width;

    const params = new Float32Array(width * height * 4);
    const colors = new Uint8Array(width * height * 4);
    for (let i = 0; i < n; i++) {
      params[i * 4] = a[i];
      params[i * 4 + 1] = b[i];
      params[i * 4 + 2] = mcia[i];
      colors[i * 4] = rgb[i * 3];
      colors[i * 4 + 1] = rgb[i * 3 + 1];
      colors[i * 4 + 2] = rgb[i * 3 + 2];
      colors[i * 4 + 3] = 255;
    }
    uploadTexture(gl, this.paramTex, gl.RGBA32F, gl.RGBA, gl.FLOAT, width, height, params);
    uploadTexture(gl, this.colorTex, gl.RGBA8, gl.RGBA, gl.UNSIGNED_BYTE, width, height, colors);
  }

  /**
   * Solve and return the canvas, ready to be used as a Plotly layout image.
   * @param {[number, number, number, number]} limits aMin, aMax, bMin, bMax
   * @param {number} singleDouble
   * @param {number} largeSuperlattice
   * @param {'voronoi'|'heatmap'} mode
   * @param {number} costScale
   */
  render(limits, singleDouble, largeSuperlattice, mode, costScale = 1) {
    if (!this.supported || this.count === 0) return null;
    const gl = /** @type {WebGL2RenderingContext} */ (this.gl);
    const p = this.program;
    gl.useProgram(p);
    gl.viewport(0, 0, this.size, this.size);

    gl.activeTexture(gl.TEXTURE0);
    gl.bindTexture(gl.TEXTURE_2D, this.paramTex);
    gl.uniform1i(gl.getUniformLocation(p, 'uParams'), 0);
    gl.activeTexture(gl.TEXTURE1);
    gl.bindTexture(gl.TEXTURE_2D, this.colorTex);
    gl.uniform1i(gl.getUniformLocation(p, 'uColors'), 1);

    gl.uniform1i(gl.getUniformLocation(p, 'uCount'), this.count);
    gl.uniform1i(gl.getUniformLocation(p, 'uTexWidth'), this.texWidth);
    gl.uniform4f(gl.getUniformLocation(p, 'uLimits'), ...limits);
    gl.uniform2f(gl.getUniformLocation(p, 'uResolution'), this.size, this.size);
    gl.uniform1f(gl.getUniformLocation(p, 'uSingleDouble'), singleDouble);
    gl.uniform1f(gl.getUniformLocation(p, 'uLargeSuper'), largeSuperlattice);
    gl.uniform1i(gl.getUniformLocation(p, 'uMode'), mode === 'heatmap' ? 1 : 0);
    gl.uniform1f(gl.getUniformLocation(p, 'uCostScale'), costScale || 1);

    gl.drawArrays(gl.TRIANGLES, 0, 3);
    return this.canvas;
  }
}

/**
 * CPU equivalent, for machines without WebGL2. Same maths, same output.
 * Kept deliberately simple; it runs at a lower resolution.
 *
 * @param {{a: Float32Array, b: Float32Array, mcia: Float32Array, rgb: Uint8Array}} data
 * @param {[number, number, number, number]} limits
 * @param {number} singleDouble
 * @param {number} largeSuperlattice
 * @param {'voronoi'|'heatmap'} mode
 * @param {number} resolution
 * @returns {HTMLCanvasElement}
 */
export function solveOnCpu(data, limits, singleDouble, largeSuperlattice, mode, resolution = 160) {
  const canvas = document.createElement('canvas');
  canvas.width = canvas.height = resolution;
  const ctx = /** @type {CanvasRenderingContext2D} */ (canvas.getContext('2d'));
  const img = ctx.createImageData(resolution, resolution);
  const [aMin, aMax, bMin, bMax] = limits;
  const n = data.a.length;

  let maxCost = 0;
  const costs = mode === 'heatmap' ? new Float32Array(resolution * resolution) : null;
  const winners = new Int32Array(resolution * resolution);

  for (let y = 0; y < resolution; y++) {
    const filmB = bMin + ((bMax - bMin) * y) / (resolution - 1);
    for (let x = 0; x < resolution; x++) {
      const filmA = aMin + ((aMax - aMin) * x) / (resolution - 1);
      let best = Infinity;
      let bestIdx = 0;
      for (let i = 0; i < n; i++) {
        const c = cost2d(
          mismatch(data.a[i], filmA),
          mismatch(data.b[i], filmB),
          data.mcia[i],
          singleDouble,
          largeSuperlattice
        );
        if (c < best) { best = c; bestIdx = i; }
      }
      const cell = y * resolution + x;
      winners[cell] = bestIdx;
      if (costs) { costs[cell] = best; if (best > maxCost) maxCost = best; }
    }
  }

  for (let cell = 0; cell < resolution * resolution; cell++) {
    // Canvas rows run top-down; the plot's y axis runs bottom-up.
    const x = cell % resolution;
    const y = Math.floor(cell / resolution);
    const out = ((resolution - 1 - y) * resolution + x) * 4;
    if (costs) {
      const [r, g, b] = tempsRamp(costs[cell] / (maxCost || 1));
      img.data[out] = r; img.data[out + 1] = g; img.data[out + 2] = b;
    } else {
      const i = winners[cell];
      img.data[out] = data.rgb[i * 3];
      img.data[out + 1] = data.rgb[i * 3 + 1];
      img.data[out + 2] = data.rgb[i * 3 + 2];
    }
    img.data[out + 3] = 255;
  }
  ctx.putImageData(img, 0, 0);
  return canvas;
}

/** @param {number} t @returns {[number, number, number]} */
function tempsRamp(t) {
  const stops = [
    [2, 87, 106], [102, 179, 168], [247, 219, 161], [227, 120, 82], [168, 36, 48],
  ];
  t = Math.max(0, Math.min(1, t)) * (stops.length - 1);
  const i = Math.min(stops.length - 2, Math.floor(t));
  const f = t - i;
  return /** @type {[number, number, number]} */ (
    stops[i].map((v, k) => Math.round(v + (stops[i + 1][k] - v) * f))
  );
}

/**
 * @param {WebGL2RenderingContext} gl
 * @param {WebGLTexture|null} tex
 */
function uploadTexture(gl, tex, internalFormat, format, type, width, height, data) {
  gl.bindTexture(gl.TEXTURE_2D, tex);
  gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MIN_FILTER, gl.NEAREST);
  gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MAG_FILTER, gl.NEAREST);
  gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_S, gl.CLAMP_TO_EDGE);
  gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_T, gl.CLAMP_TO_EDGE);
  gl.texImage2D(gl.TEXTURE_2D, 0, internalFormat, width, height, 0, format, type, data);
}

/**
 * @param {WebGL2RenderingContext} gl
 * @param {string} vertSrc
 * @param {string} fragSrc
 */
function linkProgram(gl, vertSrc, fragSrc) {
  const compile = (type, src) => {
    const sh = gl.createShader(type);
    if (!sh) return null;
    gl.shaderSource(sh, src);
    gl.compileShader(sh);
    if (!gl.getShaderParameter(sh, gl.COMPILE_STATUS)) {
      console.warn('voronoi: shader failed to compile', gl.getShaderInfoLog(sh));
      return null;
    }
    return sh;
  };
  const vs = compile(gl.VERTEX_SHADER, vertSrc);
  const fs = compile(gl.FRAGMENT_SHADER, fragSrc);
  if (!vs || !fs) return null;
  const p = gl.createProgram();
  if (!p) return null;
  gl.attachShader(p, vs);
  gl.attachShader(p, fs);
  gl.linkProgram(p);
  if (!gl.getProgramParameter(p, gl.LINK_STATUS)) {
    console.warn('voronoi: program failed to link', gl.getProgramInfoLog(p));
    return null;
  }
  return p;
}
