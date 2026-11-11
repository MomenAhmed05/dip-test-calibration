"""Step 14 (DT13) seeded errors (DT1b pattern). DEV tool check, root 2026110747.

Each validator is a check that the DT13 pipeline runs before a number is read. Each seeded error is a
deliberate bug. A bug counts as CAUGHT iff at least one validator fails on it; the clean pipeline must pass all.
Validators:
  fixtures  - CV-1 dip fixtures (allow_zero, ties, scale);
  quantile  - inverted-CDF order statistic;
  exact_n4  - MC (40,000 draws, DEV) of W(0.1, 0.005) and W(0.4, 0.005) at n = 4 at the EXACT cutoff matches the
              exact Part X value within 4.5 SE at all three alpha;
  anchor    - uniform at its own reference quantile (n = 30, 20,000 draws) within 4.5 SE of alpha;
  gate      - unimodality gate on the law list (density nondecreasing then nonincreasing; N-1);
  law       - top mass of the frozen witness within 0.005 of 0.1;
  roots     - confirmation roots declared, distinct, not DEV;
  direction - the enclosure uses counts at X_(u) for the lower and X_(l) for the upper limit;
  strict    - rejection is D >= c (CV-4), not D > c;
  provenance- every size lower bound comes from a CONFIRM-root record, not from the DEV search.
Usage: python dt13_seeded.py -> deposits/dt13/seeded_errors.json
"""
import json
import sys
from pathlib import Path

import diptest
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from dt0c_screen import cp_lower, cp_upper, g_ramp, g_uniform, g_witness  # noqa: E402, I001

W = Path(__file__).resolve().parents[1]
OUT = W / "deposits" / "dt13"
SEED = 2026110747
DECLARED = [2026110750, 2026110751, 2026110752, 2026110753, 2026110754, 2026110755, 2026110756]
BANNED_DEV = [2026110740, 2026110741, 2026110742, 2026110743, 2026110744, 2026110745, 2026110746, 2026110747,
              2026110602, 2026110711, 2026110722, 2026110734]


def dip_clean(x):
    x = np.sort(np.asarray(x, float))
    return 0.0 if x[0] == x[-1] else float(diptest.dipstat(x, allow_zero=False, sort_x=False))  # type: ignore[arg-type]


def qhat_clean(r, a):
    return r[int(np.ceil((1 - a) * len(r))) - 1]


def reject_clean(d, c):
    return d >= c


def enclosure_clean(k_Xl, k_q, k_Xu, N, conf):
    return cp_lower(k_Xu, N, conf), cp_upper(k_Xl, N, conf)


def q4_exact(a):
    r = (np.sqrt(1 + 8 / a) - 3) / 2
    return r / (4 * (r + 1))


EXACT4 = {}   # filled from deposits/dt13/n4_exact.json


def v_fixtures(dip):
    ok = abs(dip([0, 1, 2, 3]) - 0.125) < 1e-12
    ok &= abs(dip([0, 0, 1, 1, 2]) - 0.2) < 1e-12
    ok &= dip([2, 2, 2, 2, 2]) == 0.0
    ok &= abs(dip([0, 1]) - 0.25) < 1e-12
    return bool(ok)


def v_quantile(qhat):
    r = np.arange(1, 1001) / 1000.0
    return bool(all(qhat(r, a) == np.quantile(r, 1 - a, method="inverted_cdf") for a in (0.01, 0.05, 0.10)))


def v_exact_n4(dip, cutoff4, reject):
    rng = np.random.Generator(np.random.PCG64(np.random.SeedSequence(SEED, spawn_key=(1,))))
    for name, w in (("W(0.1,0.005)", 0.1), ("W(0.4,0.005)", 0.4)):
        x = g_witness(w=w, eps=0.005)(rng, (40_000, 4))
        d = np.array([dip(row) for row in x])
        for a in (0.01, 0.05, 0.10):
            p = EXACT4[(name, a)]
            m = float(np.mean(reject(d, cutoff4(a))))
            if abs(m - p) > 4.5 * np.sqrt(p * (1 - p) / len(d)):
                return False
    return True


def v_anchor(dip, ref_for_n, reject):
    rng = np.random.Generator(np.random.PCG64(np.random.SeedSequence(SEED, spawn_key=(2,))))
    ref = ref_for_n(30)
    d = np.array([dip(rng.random(30)) for _ in range(20_000)])
    for a in (0.05, 0.10):
        m = float(np.mean(reject(d, qhat_clean(ref, a))))
        if abs(m - a) > 4.5 * np.sqrt(a * (1 - a) / len(d)):
            return False
    return True


def density_unimodal(f, grid=np.linspace(1e-6, 1 - 1e-6, 20001)):
    v = f(grid)
    dv = np.diff(v)
    i = int(np.argmax(v))
    return bool(np.all(dv[:i] >= -1e-9) and np.all(dv[i:] <= 1e-9))


def v_gate(laws):
    return all(density_unimodal(f) for f in laws)


def v_law(gen):
    rng = np.random.Generator(np.random.PCG64(np.random.SeedSequence(SEED, spawn_key=(3,))))
    x = gen(rng, (1, 400_000))[0]
    return bool(abs(np.mean(x >= 0.995) - 0.1) < 0.005)


def v_roots(roots):
    return bool(set(roots) <= set(DECLARED) and len(set(roots)) == len(roots) and not set(roots) & set(BANNED_DEV))


def v_direction(enc):
    """Counts satisfy k_Xl >= k_q >= k_Xu (X_(l) <= q <= X_(u)); lower must use k_Xu, upper k_Xl."""
    lo, up = enc(540, 520, 500, 10_000, 0.999)
    return bool(abs(lo - cp_lower(500, 10_000, 0.999)) < 1e-15 and abs(up - cp_upper(540, 10_000, 0.999)) < 1e-15)


def v_strict(reject):
    return bool(reject(0.2, 0.2))


def v_provenance(lrows):
    return all(r["root"] in DECLARED for r in lrows)


def two_step_density(w, eps):
    L = 1 - eps
    return lambda x: np.where(x < L, (1 - w) / L, w / eps)


def ramp_density(w=0.1, eps_r=0.01):
    Lp = 1 - eps_r
    d0 = (1 - w) / Lp
    ctop = 2 * w / eps_r - d0
    return lambda x: np.where(x < Lp, d0, d0 + (ctop - d0) * (x - Lp) / eps_r)


LAWS_CLEAN = [two_step_density(0.1, 0.005), ramp_density(), two_step_density(0.4, 0.005), two_step_density(0.2, 1e-4)]


def run_all(dip=dip_clean, qhat=qhat_clean, reject=reject_clean, enc=enclosure_clean, cutoff4=q4_exact,
            ref_for_n=None, laws=None, gen=None, roots=(2026110750, 2026110751, 2026110752), lrows=None):
    rng = np.random.Generator(np.random.PCG64(np.random.SeedSequence(SEED, spawn_key=(4,))))
    ref30 = np.sort([dip_clean(rng.random(30)) for _ in range(100_000)])
    ref_for_n = ref_for_n or (lambda n: ref30)
    laws = LAWS_CLEAN if laws is None else laws
    gen = gen or g_witness(w=0.1, eps=0.005)
    lrows = lrows if lrows is not None else [{"root": 2026110752}, {"root": 2026110751}]
    return {"fixtures": v_fixtures(dip), "quantile": v_quantile(qhat), "exact_n4": v_exact_n4(dip, cutoff4, reject),
            "anchor": v_anchor(dip, ref_for_n, reject), "gate": v_gate(laws), "law": v_law(gen),
            "roots": v_roots(list(roots)), "direction": v_direction(enc), "strict": v_strict(reject),
            "provenance": v_provenance(lrows)}


def main():
    ex = json.loads((OUT / "n4_exact.json").read_text())
    for e in ex["confirm_laws"]:
        EXACT4[(e["law"], e["alpha"])] = e["P"][0]
    out: dict = {"clean": run_all()}
    rng = np.random.Generator(np.random.PCG64(np.random.SeedSequence(SEED, spawn_key=(5,))))
    ref50 = np.sort([dip_clean(rng.random(50)) for _ in range(100_000)])
    errs = {
        "E1 allow_zero=True (diptest default)": {"dip": lambda x: float(diptest.dipstat(np.sort(np.asarray(x, float)), allow_zero=True))},  # type: ignore[arg-type]
        "E2 dip returned x2 (excess-mass convention)": {"dip": lambda x: 2 * dip_clean(x)},
        "E3 critical value from the wrong n (n=50 reference used at n=30)": {"ref_for_n": lambda n: ref50},
        "E4 n=4 cutoff from the wrong level (alpha/2)": {"cutoff4": lambda a: q4_exact(a / 2)},
        "E5 rejection D > c instead of D >= c": {"reject": lambda d, c: np.greater(d, c)},
        "E6 bimodal law slipped into the null family (U-shaped two-step)": {"laws": LAWS_CLEAN + [lambda x: np.where((x < 0.1) | (x > 0.9), 4.0, 0.2)]},
        "E7 in-search (DEV) maximum reported as a size lower bound": {"lrows": [{"root": 2026110743}, {"root": 2026110751}]},
        "E8 DEV root reused for confirmation": {"roots": (2026110743, 2026110751, 2026110752)},
        "E9 flipped enclosure (lower from X_(l) counts, upper from X_(u))": {"enc": lambda kl, kq, ku, N, c: (cp_lower(kl, N, c), cp_upper(ku, N, c))},
        "E10 off-by-one quantile index": {"qhat": lambda r, a: r[int(np.ceil((1 - a) * len(r)))]},
        "E11 witness top mass 0.01 instead of 0.1": {"gen": g_witness(w=0.01, eps=0.005)},
        "E12 same root for reference and draws": {"roots": (2026110750, 2026110750)},
        "E13 ramp generator where W(0.1) is expected": {"gen": g_ramp()},
        "E14 uniform generator used as the witness": {"gen": g_uniform},
    }
    out["seeded"] = {}
    for name, kw in errs.items():
        r = run_all(**kw)
        out["seeded"][name] = {"validators": r, "caught": not all(r.values()), "caught_by": [k for k, v in r.items() if not v]}
        print(name, "->", "CAUGHT by " + ",".join(out["seeded"][name]["caught_by"]) if out["seeded"][name]["caught"] else "MISSED", flush=True)
    caught = sum(v["caught"] for v in out["seeded"].values())
    out["catch_rate"] = f"{caught}/{len(errs)}"
    out["clean_all_pass"] = all(out["clean"].values())
    print("clean:", out["clean"], "| catch rate", out["catch_rate"])
    try:
        (OUT / "seeded_errors.json").write_text(json.dumps(out, indent=1))
    except OSError as err:
        sys.exit(f"cannot write: {err}")


if __name__ == "__main__":
    main()
