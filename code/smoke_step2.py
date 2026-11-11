"""Step 2 remaining smokes (DT0-b). DEV tier. Writes deposits/dt0b/smoke_step2.json.

1. diptest (D_alg, CV-1) uniform percentiles vs H&H 1985 Table 1 at n = 50, 100, 200.
   n = 50: new reference, root 2026110541 (CV-6), 200,000 dips. n = 100, 200: the DT0-c reference (1e6 each).
2. Interval / multiprecision libraries: mpmath.iv, python-flint arb (if installed).
3. PyTorch CPU: time one adversary-style training step at n = 200.
"""

import json
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from dt0c_screen import dips_of  # noqa: E402, I001

W = Path(__file__).resolve().parents[1]
OUT = W / "deposits" / "dt0b"
P = (0.01, 0.05, 0.10, 0.50, 0.90, 0.95, 0.99, 0.995, 0.999)
HH = {  # H&H 1985 Table 1 (OCR, checked against the page): 9,999 dips, max SE .001
    50: (.0312, .0352, .0378, .0489, .0645, .0702, .0804, .0842, .0926),
    100: (.0228, .0256, .0274, .0355, .0471, .0510, .0586, .0619, .0687),
    200: (.0165, .0185, .0197, .0255, .0341, .0370, .0429, .0449, .0496),
}


def table1():
    out = {}
    rng = np.random.Generator(np.random.PCG64(np.random.SeedSequence(2026110541, spawn_key=(0, 50))))
    ref50 = np.concatenate([dips_of(rng.random((2000, 50))) for _ in range(100)])
    refs = {50: ref50,
            100: np.load(W / "deposits" / "dt0c" / "ref_n100.npy"),
            200: np.load(W / "deposits" / "dt0c" / "ref_n200.npy")}
    for n, ref in refs.items():
        ours = np.quantile(ref, P, method="inverted_cdf")
        diff = ours - np.array(HH[n])
        out[n] = {"R": len(ref), "ours": [round(float(x), 5) for x in ours], "hh": HH[n],
                  "max_abs_diff": float(np.max(np.abs(diff))),
                  "within_0.002": bool(np.all(np.abs(diff) <= 0.002))}
    return out


def intervals():
    res = {}
    try:
        import mpmath
        from mpmath import iv
        iv.dps = 30
        x = iv.mpf(1) / 3
        y = x * 3
        res["mpmath"] = {"version": mpmath.__version__, "third_times_3": str(y),
                         "contains_1": bool(1 in y), "sqrt2_interval": str(iv.sqrt(2))}
    except Exception as e:  # noqa: BLE001
        res["mpmath"] = {"error": repr(e)}
    try:
        import flint
        flint.ctx.prec = 100
        a = flint.arb(1) / 3 * 3
        res["flint"] = {"version": flint.__version__, "third_times_3": str(a), "contains_1": bool(a.contains(1))}
    except Exception as e:  # noqa: BLE001
        res["flint"] = {"error": repr(e)}
    return res


def torch_step():
    try:
        import torch
    except Exception as e:  # noqa: BLE001
        return {"error": repr(e)}
    torch.manual_seed(2026110541)
    n, batch, k = 200, 64, 64
    # generator: maps noise to a sorted monotone sample via cumulative softplus increments (a unimodal-ish
    # parametric law is NOT enforced here; this only times the forward/backward pass of a small MLP + sort).
    net = torch.nn.Sequential(torch.nn.Linear(k, 256), torch.nn.ReLU(), torch.nn.Linear(256, n))
    opt = torch.optim.Adam(net.parameters(), lr=1e-3)
    times = []
    for _ in range(20):
        t0 = time.perf_counter()
        z = torch.randn(batch, k)
        x = torch.cumsum(torch.nn.functional.softplus(net(z)), dim=1)
        x = (x - x[:, :1]) / (x[:, -1:] - x[:, :1] + 1e-9)
        ecdf = torch.arange(1, n + 1) / n
        loss = -(ecdf - x).abs().amax(dim=1).mean()  # a sup-distance surrogate
        opt.zero_grad()
        loss.backward()
        opt.step()
        times.append(time.perf_counter() - t0)
    return {"torch": torch.__version__, "cuda": torch.cuda.is_available(), "threads": torch.get_num_threads(),
            "n": n, "batch": batch, "median_step_ms": 1000 * float(np.median(times[5:]))}


def main():
    t0 = time.time()
    res = {"tier": "DEV", "table1": table1(), "intervals": intervals(), "torch": torch_step()}
    res["seconds"] = time.time() - t0
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "smoke_step2.json").write_text(json.dumps(res, indent=1))
    print(json.dumps(res, indent=1))


if __name__ == "__main__":
    main()
