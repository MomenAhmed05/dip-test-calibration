"""F_M from the published Balabdaoui & Filali table (Table 2, x = 0.33 .. 2.54, step 0.01), strict range.

F_M = cdf of M = sup(Brownian bridge - its concave majorant) (Balabdaoui & Pitman 2011).
Replaces the guessed extension 0.393*(x/0.9)**8 in dt2_farm.FM below x = 0.90 (DT3c-RESULT §2).

  python code/fm_table.py --build   # parse the local paper text once -> deposits/dt4/bf_table_full.json
  python code/fm_table.py --test    # regression tests (table shape, known values, vs fm_check.json MC)
"""
import json
import re
import sys
from pathlib import Path

import numpy as np

W = Path(__file__).resolve().parents[1]
SRC = W.parent / "prior" / "priorart" / "step2" / "dl_balabdaoui_filali_2011" / "bf.txt"
TAB = W / "deposits" / "dt4" / "bf_table_full.json"
X_LO, X_HI = 0.33, 2.54
TAIL_LO, TAIL_HI = 9.25e-12, 1.44e-5  # F_M(0.33) and 1 - F_M(2.54) from the same table


def build():
    """Parse the x/F_M pairs of Table 2 from the paper text (numbers only; nothing is executed)."""
    pairs = {}
    for line in SRC.read_text(encoding="utf-8", errors="replace").splitlines():
        toks = line.split()
        if len(toks) >= 2 and len(toks) % 2 == 0 and all(re.fullmatch(r"\d\.\d\d", toks[j]) for j in range(0, len(toks), 2)):
            for j in range(0, len(toks), 2):
                x = float(toks[j])
                if X_LO - 1e-9 <= x <= X_HI + 1e-9:
                    pairs.setdefault(round(x, 2), float(toks[j + 1]))
    xs = sorted(pairs)
    TAB.parent.mkdir(parents=True, exist_ok=True)
    TAB.write_text(json.dumps({"source": "Balabdaoui & Filali, Table 2 (bf.txt)", "x": xs,
                               "F_M": [pairs[x] for x in xs]}, indent=0))
    print("built", TAB, len(xs), "points", xs[0], "..", xs[-1])


_T = None


def _table():
    global _T
    if _T is None:
        try:
            d = json.loads(TAB.read_text())
        except (OSError, json.JSONDecodeError) as e:
            raise RuntimeError(f"cannot read {TAB}: {e}; run fm_table.py --build first") from e
        _T = (np.array(d["x"]), np.array(d["F_M"]))
    return _T


def FM_pub(x, clip_tails=False):
    """F_M(x) by linear interpolation in the published table.

    Outside [0.33, 2.54] raise, unless clip_tails=True: then return 0 below (error <= 9.25e-12)
    and 1 above (error <= 1.44e-5). Never extrapolates by a guessed formula.
    """
    xs, fv = _table()
    if x < X_LO or x > X_HI:
        if not clip_tails:
            raise ValueError(f"F_M requested at x={x:.4f}, outside the published table [{X_LO}, {X_HI}]")
        return 0.0 if x < X_LO else 1.0
    return float(np.interp(x, xs, fv))


def test():
    xs, fv = _table()
    assert len(xs) == 222 and abs(xs[0] - X_LO) < 1e-12 and abs(xs[-1] - X_HI) < 1e-12, "table extent"
    assert np.allclose(np.diff(xs), 0.01), "uniform grid"
    assert np.all(np.diff(fv) > 0), "F_M must increase"
    for x, v in [(0.80, 0.220643921928), (1.00, 0.560961512572), (1.79, 0.993859705663)]:
        assert abs(FM_pub(x) - v) < 1e-12, x
    # agreement with the old 0.90-1.79 transcription (dt1b_confirm.bf_table)
    sys.path.insert(0, str(Path(__file__).parent))
    from dt1b_confirm import bf_table
    ox, ov = bf_table()
    assert max(abs(FM_pub(float(x)) - v) for x, v in zip(ox, ov, strict=True)) < 1e-12, "old transcription"
    # independent Monte Carlo check (fm_check.py, N = 1e7)
    try:
        mc = json.loads((W / "deposits" / "dt4" / "fm_check.json").read_text())
    except (OSError, json.JSONDecodeError) as e:
        raise RuntimeError(f"cannot read fm_check.json: {e}") from e
    d = max(abs(FM_pub(r["x"]) - r["F_M_est"]) for r in mc["rows"])
    assert d < mc["hoeffding_eps"], f"MC disagreement {d}"
    for bad in (0.2, 3.0):
        try:
            FM_pub(bad)
            raise AssertionError("no error outside range")
        except ValueError:
            pass
    print(f"fm_table tests PASS: 222 points, max |table - MC| = {d:.2e} (< Hoeffding {mc['hoeffding_eps']:.1e})")


if __name__ == "__main__":
    if "--build" in sys.argv:
        build()
    if "--test" in sys.argv:
        test()
