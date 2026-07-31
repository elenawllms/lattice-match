"""Round-trip tests for the binary bundle format.

The browser reads these buffers as zero-copy typed-array views, so an offset or
dtype that is off by one silently yields plausible-looking garbage rather than
an error. These tests decode the bundle the same way the TypeScript loader does
and compare against the source frame.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from pipeline.encode import (
    ELEMENT_INDEX,
    ELEMENTS,
    MASK_WORDS,
    BundleWriter,
    encode_films,
    encode_superlattices,
)

REPO = Path(__file__).resolve().parents[2]
BUNDLES = REPO / "web" / "public" / "data"

NUMPY_DTYPES = {
    "float32": np.float32, "float64": np.float64,
    "uint8": np.uint8, "uint16": np.uint16, "uint32": np.uint32, "int8": np.int8,
}


def read_column(blob: bytes, meta: dict, name: str) -> np.ndarray:
    """Decode one column exactly as the TypeScript loader does."""
    col = meta["columns"][name]
    dtype = NUMPY_DTYPES[col["dtype"]]
    return np.frombuffer(blob, dtype=dtype, count=col["length"], offset=col["offset"])


def decode_dictionary(blob: bytes, meta: dict, name: str) -> list[str]:
    codes = read_column(blob, meta, name)
    table = meta["dictionaries"][name]
    return [table[c] for c in codes]


@pytest.fixture(scope="module")
def superlattice_frame() -> pd.DataFrame:
    return pd.read_csv(REPO / "src" / "assets" / "data" / "sublattices_2d.csv")


def test_superlattice_numerics_round_trip(superlattice_frame):
    blob, meta = encode_superlattices(superlattice_frame)
    for col in ("a", "b", "mcia"):
        got = read_column(blob, meta, col)
        # float32 keeps ~7 significant digits; lattice parameters need ~5.
        assert np.allclose(got, superlattice_frame[col].to_numpy(), rtol=1e-6)


def test_float32_precision_is_adequate_for_lattice_parameters(superlattice_frame):
    """Guard the format choice itself: the error must be far below the ~0.01%
    mismatch differences the tool is asked to distinguish."""
    a = superlattice_frame["a"].to_numpy()
    rel = np.abs(a.astype(np.float32).astype(np.float64) - a) / a
    assert rel.max() < 1e-6


def test_superlattice_strings_round_trip(superlattice_frame):
    blob, meta = encode_superlattices(superlattice_frame)
    for col in ("substrate", "dimensions", "angle"):
        assert decode_dictionary(blob, meta, col) == list(superlattice_frame[col])


def test_colours_round_trip_as_bytes(superlattice_frame):
    blob, meta = encode_superlattices(superlattice_frame)
    for column, name in (("R", "red"), ("G", "green"), ("B", "blue")):
        got = read_column(blob, meta, name)
        expected = np.round(superlattice_frame[column].to_numpy() * 255)
        assert np.array_equal(got, expected.astype(np.uint8))


def test_every_column_is_aligned_for_typed_array_views(superlattice_frame):
    """A TypedArray view throws unless byteOffset is a multiple of BYTES_PER_ELEMENT."""
    blob, meta = encode_superlattices(superlattice_frame)
    for name, col in meta["columns"].items():
        size = np.dtype(NUMPY_DTYPES[col["dtype"]]).itemsize
        assert col["offset"] % size == 0, f"{name} is misaligned for {col['dtype']}"


def test_columns_do_not_overlap(superlattice_frame):
    blob, meta = encode_superlattices(superlattice_frame)
    spans = sorted(
        (c["offset"], c["offset"] + c["length"] * np.dtype(NUMPY_DTYPES[c["dtype"]]).itemsize)
        for c in meta["columns"].values()
    )
    for (_, end), (start, _) in zip(spans, spans[1:]):
        assert start >= end, "columns overlap in the buffer"
    assert spans[-1][1] <= len(blob)


# --- element bitmasks ------------------------------------------------------


def mask_of(blob: bytes, meta: dict, row: int) -> int:
    words = read_column(blob, meta, "element_mask")
    value = 0
    for w in range(MASK_WORDS):
        value |= int(words[row * MASK_WORDS + w]) << (32 * w)
    return value


def test_element_masks_encode_the_right_bits():
    df = pd.DataFrame([
        dict(formula="GaN", elements="Ga, N", num_elements=2, a=3.19, b=5.5,
             plane="(0001)", crystal_system="Hexagonal", point_group="6mm"),
        dict(formula="H2O", elements="H, O", num_elements=2, a=4.0, b=4.0,
             plane="(001)", crystal_system="Cubic", point_group="m-3m"),
    ])
    blob, meta = encode_films(df)
    assert mask_of(blob, meta, 0) == (1 << ELEMENT_INDEX["Ga"]) | (1 << ELEMENT_INDEX["N"])
    assert mask_of(blob, meta, 1) == (1 << ELEMENT_INDEX["H"]) | (1 << ELEMENT_INDEX["O"])


def test_mask_words_cover_every_element():
    assert len(ELEMENTS) <= MASK_WORDS * 32
    assert ELEMENT_INDEX["H"] == 0 and ELEMENT_INDEX["Og"] == len(ELEMENTS) - 1


def test_mask_semantics_match_the_filters():
    """The three film filters, expressed as the bitwise ops the browser runs."""
    df = pd.DataFrame([
        dict(formula="GaN", elements="Ga, N", num_elements=2, a=3.19, b=5.5,
             plane="(0001)", crystal_system="Hexagonal", point_group="6mm"),
    ])
    blob, meta = encode_films(df)
    row = mask_of(blob, meta, 0)
    ga, n, o = (1 << ELEMENT_INDEX[s] for s in ("Ga", "N", "O"))

    assert row & (ga | n) == (ga | n)      # must-include Ga and N -> passes
    assert row & o == 0                    # exclude O            -> passes
    assert row & ~(ga | n | o) == 0        # can-include {Ga,N,O} -> passes
    assert row & ~(ga | o) != 0            # can-include {Ga,O}   -> fails, has N


# --- the shipped bundles ---------------------------------------------------


@pytest.mark.skipif(not (BUNDLES / "manifest.json").exists(), reason="bundles not built")
def test_shipped_bundles_decode_and_match_the_csvs():
    manifest = json.loads((BUNDLES / "manifest.json").read_text())
    pairs = {
        "superlattices_2d": "sublattices_2d.csv",
        "superlattices_1d": "sublattices_1d.csv",
        "films_2d": "stable_films_2d.csv",
        "films_1d": "stable_films_1d.csv",
    }
    for table, csv_name in pairs.items():
        entry = manifest["tables"][table]
        meta = json.loads((BUNDLES / entry["meta"]).read_text())
        blob = (BUNDLES / entry["file"]).read_bytes()
        df = pd.read_csv(REPO / "src" / "assets" / "data" / csv_name)

        assert meta["rows"] == len(df), f"{table}: bundle has {meta['rows']} rows, CSV has {len(df)}"
        assert np.allclose(read_column(blob, meta, "a"), df["a"].to_numpy(), rtol=1e-6)
        if "b" in meta["columns"]:
            assert np.allclose(read_column(blob, meta, "b"), df["b"].to_numpy(), rtol=1e-6)


@pytest.mark.skipif(not (BUNDLES / "manifest.json").exists(), reason="bundles not built")
def test_only_superlattices_are_eager():
    """Film data must not load on first paint; it is ~730 KB gzipped."""
    manifest = json.loads((BUNDLES / "manifest.json").read_text())
    eager = {n for n, t in manifest["tables"].items() if t["eager"]}
    assert eager == {"superlattices_2d", "superlattices_1d"}


def test_dictionary_index_width_scales_with_cardinality():
    small = pd.Series(["a", "b", "a"])
    w = BundleWriter()
    w.dictionary("x", small)
    assert w.columns["x"]["dtype"] == "uint8"

    big = pd.Series([f"v{i}" for i in range(300)])
    w2 = BundleWriter()
    w2.dictionary("y", big)
    assert w2.columns["y"]["dtype"] == "uint16"


def test_duplicate_column_names_are_rejected():
    """The blue channel was lowercased to "b" and silently overwrote the b
    lattice parameter in the manifest, so the browser would have decoded
    colour bytes as lattice constants."""
    w = BundleWriter()
    w.numeric("b", np.array([1.0, 2.0]))
    with pytest.raises(ValueError, match="duplicate column"):
        w.numeric("b", np.array([3, 4]), "uint8")
