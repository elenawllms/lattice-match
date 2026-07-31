"""Pack the generated tables into compact binary bundles for the web app.

    python -m pipeline.encode --data src/assets/data --out web/data

A bundle is a JSON manifest plus one flat ``.bin`` of concatenated typed
arrays. The manifest records each column's dtype, byte offset and length, so
the browser can create zero-copy typed-array views over a single fetch.

Why not JSON: ``sublattices_2d`` is 1.83 MB as JSON records but only ~100 KB
here, because the numbers travel as float32 rather than decimal text and the
repeated strings become dictionary indices.

Why not Arrow: the schema is small and fixed, and the JS Arrow reader is
~200 KB gzipped -- more than the data it would be decoding.

Element membership is packed as a 128-bit mask (4 x uint32) per film, so the
must/can/exclude filters become two bitwise ops per row in the browser rather
than the six row-wise Python string scans the Dash app ran over 81,000 rows.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

#: Element symbols in atomic-number order; index is the bit position.
ELEMENTS: tuple[str, ...] = (
    "H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne", "Na", "Mg", "Al",
    "Si", "P", "S", "Cl", "Ar", "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe",
    "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr",
    "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn",
    "Sb", "Te", "I", "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm",
    "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf", "Ta", "W",
    "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", "At", "Rn",
    "Fr", "Ra", "Ac", "Th", "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf",
    "Es", "Fm", "Md", "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds",
    "Rg", "Cn", "Nh", "Fl", "Mc", "Lv", "Ts", "Og",
)
ELEMENT_INDEX = {sym: i for i, sym in enumerate(ELEMENTS)}

#: uint32 words per element mask. 118 elements need 4.
MASK_WORDS = 4


class BundleWriter:
    """Accumulates columns into one buffer and emits a manifest describing it."""

    def __init__(self) -> None:
        self._chunks: list[bytes] = []
        self._offset = 0
        self.columns: dict[str, dict] = {}
        self.dictionaries: dict[str, list[str]] = {}

    def _add(self, name: str, array: np.ndarray, dtype: str) -> None:
        if name in self.columns:
            # Silent clobbering here is invisible: the bytes of both columns
            # land in the buffer, but the manifest points at the later one, so
            # the browser decodes the wrong bytes as plausible-looking values.
            raise ValueError(f"duplicate column name {name!r} in bundle")
        data = np.ascontiguousarray(array, dtype=dtype)
        # Typed-array views require the offset to be a multiple of the element
        # size, so pad to keep every column aligned.
        align = data.dtype.itemsize
        if pad := (-self._offset) % align:
            self._chunks.append(b"\0" * pad)
            self._offset += pad
        self.columns[name] = {
            "dtype": dtype, "offset": self._offset, "length": int(data.size)
        }
        raw = data.tobytes()
        self._chunks.append(raw)
        self._offset += len(raw)

    def numeric(self, name: str, values, dtype: str = "float32") -> None:
        self._add(name, np.asarray(values), dtype)

    def dictionary(self, name: str, values: pd.Series) -> None:
        """Store repeated strings as indices into a shared string table."""
        cat = pd.Categorical(values)
        table = [str(c) for c in cat.categories]
        dtype = "uint8" if len(table) < 256 else "uint16" if len(table) < 65536 else "uint32"
        if len(table) >= 2**32:
            raise ValueError(f"{name}: too many distinct values to index")
        self._add(name, cat.codes, dtype)
        self.dictionaries[name] = table

    def element_masks(self, name: str, element_lists: pd.Series) -> None:
        """Pack each row's element set into MASK_WORDS uint32s, row-major."""
        masks = np.zeros((len(element_lists), MASK_WORDS), dtype=np.uint32)
        for row, raw in enumerate(element_lists):
            for sym in str(raw).split(", "):
                idx = ELEMENT_INDEX.get(sym.strip())
                if idx is None:
                    continue
                masks[row, idx // 32] |= np.uint32(1 << (idx % 32))
        self._add(name, masks.ravel(), "uint32")

    def finish(self, rows: int) -> tuple[bytes, dict]:
        return b"".join(self._chunks), {
            "rows": rows,
            "columns": self.columns,
            "dictionaries": self.dictionaries,
        }


def encode_superlattices(df: pd.DataFrame) -> tuple[bytes, dict]:
    w = BundleWriter()
    # float64, not float32. mismatch() is (film - sub) / sub, a difference of
    # two nearby numbers, so float32's ~6e-8 relative error on `a` is amplified
    # by ~5e4 in the result -- around 1e-4 relative on the mismatch itself.
    # That is enough to reorder near-equal matches, which are exactly the ones
    # this tool exists to rank. A differential test against the Python
    # implementation caught two adjacent entries swapping places.
    w.numeric("a", df["a"], "float64")
    if "b" in df.columns:
        w.numeric("b", df["b"], "float64")
    w.numeric("mcia", df["mcia"], "float64")
    # Colours as uint8 RGB; the CSV carries both float channels and a
    # pre-rendered "rgb(...)" string, which is the same information three times.
    # Spelled out rather than lowercased: "B".lower() collides with the "b"
    # lattice parameter and silently overwrote it.
    for column, name in (("R", "red"), ("G", "green"), ("B", "blue")):
        w.numeric(name, np.round(df[column].to_numpy() * 255), "uint8")
    w.dictionary("substrate", df["substrate"])
    w.dictionary("dimensions", df["dimensions"])
    w.dictionary("angle", df["angle"])
    return w.finish(len(df))


def encode_films(df: pd.DataFrame) -> tuple[bytes, dict]:
    w = BundleWriter()
    # float64 for the same reason as the superlattices: these feed the same
    # cancellation-prone mismatch(). This table is lazy-loaded, so the extra
    # bytes do not touch first paint.
    w.numeric("a", df["a"], "float64")
    if "b" in df.columns:
        w.numeric("b", df["b"], "float64")
    w.numeric("num_elements", df["num_elements"], "uint8")
    w.element_masks("element_mask", df["elements"])
    w.dictionary("formula", df["formula"])
    w.dictionary("plane", df["plane"])
    w.dictionary("crystal_system", df["crystal_system"])
    w.dictionary("point_group", df["point_group"])
    return w.finish(len(df))


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data", type=Path, required=True, help="generated CSV directory")
    parser.add_argument("--out", type=Path, required=True, help="bundle output directory")
    args = parser.parse_args(argv)
    args.out.mkdir(parents=True, exist_ok=True)

    # Each table's metadata is a sidecar file rather than one combined
    # manifest. The film formula dictionary alone is ~100 KB gzipped, and the
    # opening view needs only the superlattices -- so films stay unfetched
    # until the user opens the film panel.
    manifest: dict = {"elements": list(ELEMENTS), "maskWords": MASK_WORDS, "tables": {}}
    total_csv = total_bin = 0

    for name, fn, csv_name, eager in [
        ("superlattices_2d", encode_superlattices, "sublattices_2d.csv", True),
        ("superlattices_1d", encode_superlattices, "sublattices_1d.csv", True),
        ("films_2d", encode_films, "stable_films_2d.csv", False),
        ("films_1d", encode_films, "stable_films_1d.csv", False),
    ]:
        csv_path = args.data / csv_name
        df = pd.read_csv(csv_path)
        blob, meta = fn(df)
        meta["file"] = f"{name}.bin"
        (args.out / f"{name}.bin").write_bytes(blob)
        (args.out / f"{name}.meta.json").write_text(json.dumps(meta, separators=(",", ":")))
        manifest["tables"][name] = {
            "meta": f"{name}.meta.json", "file": meta["file"],
            "rows": meta["rows"], "eager": eager,
        }
        total_csv += csv_path.stat().st_size
        total_bin += len(blob)
        print(f"{name:18s} {len(df):7,} rows  {csv_path.stat().st_size / 1e6:6.2f} MB CSV "
              f"-> {len(blob) / 1e6:5.2f} MB bin  {'eager' if eager else 'lazy'}")

    (args.out / "manifest.json").write_text(json.dumps(manifest, separators=(",", ":")))
    print(f"{'total':18s} {'':7s}       {total_csv / 1e6:6.2f} MB      -> {total_bin / 1e6:5.2f} MB")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
