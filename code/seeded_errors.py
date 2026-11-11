"""Step 3 item 5: seeded errors (one catch rule per error). DEV tool check.

Each validator is a check that the real pipeline runs (or must run) before a number is read.
Each seeded error is a deliberate bug. A bug counts as CAUGHT iff at least one validator fails on it,
and the clean pipeline passes every validator. Writes deposits/dt1b/seeded_errors.json.
"""

import json
import sys
from pathlib import Path

import diptest
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from dt0c_screen import g_ramp, g_witness  # noqa: E402, I001

W = Path(__file__).resolve().parents[1]
REF = W / "deposits" / "dt1b"
SEED = 2026110604  # DEV root for this tool check (declared in CV-6)


# ---------------- pipeline pieces (clean versions) ----------------
def dip_clean(x):
    x = np.sort(np.asarray(x, float))
    return 0.0 if x[0] == x[-1] else float(diptest.dipstat(x, allow_zero=False, sort_x=False))  # type: ignore[arg-type]


def qhat_clean(ref_sorted, a):
    return ref_sorted[int(np.ceil((1 - a) * len(ref_sorted))) - 1]


# ---------------- validators ----------------
def v_fixtures(dip):
    """Convention fixtures (CV-1 / DT1a §4)."""
    ok = abs(dip([0, 1, 2, 3]) - 0.125) < 1e-12
    ok &= abs(dip([0, 0, 1, 1, 2]) - 0.2) < 1e-12
    ok &= dip([2, 2, 2, 2, 2]) == 0.0
    ok &= abs(dip([0, 1]) - 0.25) < 1e-12
    return bool(ok)


def v_quantile(qhat):
    """The point quantile equals numpy's inverted-CDF quantile on a known array."""
    r = np.arange(1, 1001) / 1000.0
    return bool(qhat(r, 0.05) == np.quantile(r, 0.95, method="inverted_cdf")
                and qhat(r, 0.01) == np.quantile(r, 0.99, method="inverted_cdf"))


def v_anchor(dip, ref_for_n, n=1000, reps=4000):
    """Uniform rejection at the reference's 5% point must be in a 99.9% band around 5%."""
    rng = np.random.Generator(np.random.PCG64(np.random.SeedSequence(SEED, spawn_key=(1, n))))
    ref = ref_for_n(n)
    q = qhat_clean(ref, 0.05)
    d = np.array([dip(rng.random(n)) for _ in range(reps)])
    r = np.mean(d >= q)
    se = np.sqrt(0.05 * 0.95 / reps)
    return bool(abs(r - 0.05) <= 3.3 * se)


def v_law(gen, w, L):
    """Law check: the top mass, and that the shoulder density is flat (not sloped) and below the top."""
    rng = np.random.Generator(np.random.PCG64(np.random.SeedSequence(SEED, spawn_key=(2,))))
    x = gen(rng, (1, 400_000))[0]
    top = np.mean(x >= L)
    sh = x[(x > 0) & (x < L)]
    h, _ = np.histogram(sh, bins=10, range=(0, L))
    slope = np.polyfit(np.arange(10), h / h.mean(), 1)[0]
    in_range = np.all((x >= -1e-12) & (x <= 1 + 1e-12))
    return bool(abs(top - w) < 0.005 and abs(slope) < 0.01 and in_range)


def v_roots(roots_used, declared_conf, banned_dev):
    """Provenance: confirmation roots are declared, distinct, and not on the DEV list."""
    return bool(set(roots_used) <= set(declared_conf) and len(set(roots_used)) == len(roots_used)
                and not set(roots_used) & set(banned_dev))


DECLARED = [2026110601, 2026110602, 2026110603]
BANNED = [20261103, 7, 20261107, 2026110530, 2026110540, 2026110541, 2026110550, 2026110542, 2026110560]


def ref_clean(n):
    return np.load(REF / f"ref_n{n}.npy")


def run_all(dip=dip_clean, qhat=qhat_clean, ref_for_n=ref_clean, gen=None, w=0.1, L=0.995,
            roots=(2026110601, 2026110602)):
    gen = gen or g_witness(w=0.1)
    return {"fixtures": v_fixtures(dip), "quantile": v_quantile(qhat), "anchor": v_anchor(dip, ref_for_n),
            "law": v_law(gen, w, L), "roots": v_roots(list(roots), DECLARED, BANNED)}


def main():
    out: dict = {"clean": run_all()}
    # seeded errors
    errs = {
        "E1 allow_zero=True": {"dip": lambda x: float(diptest.dipstat(np.sort(np.asarray(x, float)), allow_zero=True))},  # type: ignore[arg-type]
        "E2 off-by-one quantile index": {"qhat": lambda r, a: r[int(np.ceil((1 - a) * len(r)))]},
        "E3 ties dropped before the dip": {"dip": lambda x: dip_clean(np.unique(x)) if len(np.unique(x)) > 1 else 0.0},
        "E4 reference at the wrong n (n=2000 used for n=1000)": {"ref_for_n": lambda n: ref_clean(2000)},
        "E5 top mass 0.01 instead of 0.1": {"gen": g_witness(w=0.01)},
        "E6 shoulder slope sign flipped (beta=-1: decreasing shoulder)": {"gen": g_witness(w=0.1, beta=-1.0)},
        "E7 noise scale x10 (sd 0.05 added)": {"gen": g_witness(w=0.1, noise_sd=0.05)},
        "E8 confirmation root reused from a DEV run": {"roots": (2026110540, 2026110602)},
        "E9 same root for reference and draws": {"roots": (2026110602, 2026110602)},
        "E10 ramp law passed where W is expected (wrong top region)": {"gen": g_ramp()},
    }
    out["seeded"] = {}
    for name, kw in errs.items():
        r = run_all(**kw)
        out["seeded"][name] = {"validators": r, "caught": not all(r.values()),
                               "caught_by": [k for k, v in r.items() if not v]}
    caught = sum(v["caught"] for v in out["seeded"].values())
    out["catch_rate"] = f"{caught}/{len(errs)}"
    out["clean_all_pass"] = all(out["clean"].values())
    (REF / "seeded_errors.json").write_text(json.dumps(out, indent=1))
    print(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
