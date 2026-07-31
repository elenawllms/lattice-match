"""Run the browser cost/ranking code against the Python implementation.

The physics now exists twice: once in Python (`pipeline/`, `src/app.py`) and
once in JavaScript (`web/src/core/`). Two implementations of the same formulas
is precisely how the old app ended up with a table and a Voronoi map that
ranked matches differently, so they are pinned to each other here.

Executed with JavaScriptCore, which ships with macOS at a fixed path, so this
needs no Node installation. Skipped where `jsc` is absent.
"""

from __future__ import annotations

import base64
import json
import shutil
import subprocess
import sys
from pathlib import Path

import pandas as pd
import pytest

REPO = Path(__file__).resolve().parents[2]
WEB = REPO / "web"
JSC = Path(
    "/System/Library/Frameworks/JavaScriptCore.framework/Versions/A/Helpers/jsc"
)


def js_runtime() -> list[str] | None:
    """Whichever JS runtime is available: JavaScriptCore locally, Node in CI."""
    if JSC.exists():
        return [str(JSC), "-m"]
    node = shutil.which("node")
    return [node] if node else None


pytestmark = pytest.mark.skipif(
    js_runtime() is None or not (WEB / "data" / "superlattices_2d.bin").exists(),
    reason="needs a JS runtime (JavaScriptCore or Node) and built web bundles",
)

sys.path.insert(0, str(REPO / "src"))


def _reference(tmp_path: Path) -> Path:
    """Compute reference values with the Python implementation."""
    import app as A  # noqa: PLC0415 - imported late so the skip applies first

    cases = [
        (0.1, 0.05, 50.0, 0.5, 0.5),
        (0.0, 0.0, 10.0, 0.5, 0.5),
        (-0.02, 0.3, 120.0, 0.05, 1.0),
        (0.25, -0.25, 7.5, 1.0, 0.05),
        (1e-9, 1e-9, 200.0, 0.35, 0.75),
    ]
    ref = {
        "mismatch": [
            [s, f, A.mismatch(s, f)]
            for s, f in [(5, 5.5), (4.763, 13.003), (3.19, 3.0), (12.0, 12.0)]
        ],
        "cost2d": [[*c, float(A.costFunction2d(*c))] for c in cases],
        "cost1d": [
            [am, mc, q, float(A.costFunction1d(am, mc, q))]
            for am, mc, q in [(0.1, 50.0, 0.5), (0.0, 10.0, 0.5), (-0.2, 120.0, 0.05)]
        ],
    }

    # Rank the real catalogue against sapphire M-plane's own parameters.
    sl = pd.read_csv(REPO / "src" / "assets" / "data" / "sublattices_2d.csv")
    a, b, p, q, n = 4.763, 13.003, 0.5, 0.5, 40
    scored = sorted(
        (
            (
                i,
                float(
                    A.costFunction2d(
                        A.mismatch(row.a, a), A.mismatch(row.b, b), row.mcia, p, q
                    )
                ),
            )
            for i, row in enumerate(sl.itertuples())
        ),
        key=lambda t: t[1],
    )[:n]
    ref["ranking"] = {"a": a, "b": b, "p": p, "q": q, "top": [list(t) for t in scored]}

    # And the triangular table.
    sl1 = pd.read_csv(REPO / "src" / "assets" / "data" / "sublattices_1d.csv")
    a1, q1 = 3.19, 0.5
    scored1 = sorted(
        ((i, float(A.costFunction1d(A.mismatch(row.a, a1), row.mcia, q1)))
         for i, row in enumerate(sl1.itertuples())),
        key=lambda t: t[1],
    )[:20]
    b64_1d = tmp_path / "superlattices_1d.b64"
    b64_1d.write_text(
        base64.b64encode((WEB / "data" / "superlattices_1d.bin").read_bytes()).decode()
    )
    ref["ranking1d"] = {
        "a": a1, "q": q1, "top": [list(t) for t in scored1],
        "meta": str(WEB / "data" / "superlattices_1d.meta.json"),
        "bin": str(b64_1d),
    }

    path = tmp_path / "reference.json"
    path.write_text(json.dumps(ref))
    return path


def test_browser_matches_python(tmp_path):
    ref_path = _reference(tmp_path)

    b64 = tmp_path / "superlattices_2d.b64"
    b64.write_text(
        base64.b64encode((WEB / "data" / "superlattices_2d.bin").read_bytes()).decode()
    )

    runtime = js_runtime()
    assert runtime is not None
    script = [str(WEB / "test" / "differential.mjs")]
    # jsc needs a bare `--` to separate script arguments; node does not.
    separator = ["--"] if runtime[0] == str(JSC) else []
    args = [
        str(ref_path), str(b64),
        str(WEB / "data" / "superlattices_2d.meta.json"),
        str(WEB / "src" / "core"),
    ]
    result = subprocess.run(
        runtime + script + separator + args,
        capture_output=True, text=True, timeout=180,
    )
    assert result.returncode == 0, f"jsc failed:\n{result.stderr}"

    payload = json.loads(result.stdout.strip().splitlines()[-1])
    assert not payload["failures"], "JS and Python disagree:\n  " + "\n  ".join(
        payload["failures"]
    )
    # Guard against the harness silently checking nothing.
    assert payload["pass"] > 100, f"only {payload['pass']} assertions ran"
