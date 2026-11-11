"""O5-ML secondary, descriptive (DT5-ML-DESIGN v0.2 §1 "Real activations"): dip p-values of real ReLU activations
under the uniform table and under LF, on sklearn `digits` (1,797 images, 64 features).

A 1-hidden-layer MLP (64 ReLU units, sklearn MLPClassifier, random_state 0, max_iter 500) is trained on all data.
For each hidden unit: (a) all activations (atom at 0 + positive part: outside Theorem 4 because of the atom);
(b) the positive part only (continuous); (c) per class, positive part. For each sample: root dip, p_uniform
(diptest table) and p_LF = 1 - F_M(2 sqrt(n) dip). Reported: counts of p < 0.05 under each calibration, and the
fraction of exact zeros. There is no ground truth, so no "spurious" claim is made (design).
Output: deposits/o5_activations.json. Seed: fixed (random_state 0); no hypothesis test is pre-registered.
"""
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from fm_table import FM_pub  # noqa: E402


def pvals(x):
    import diptest
    x = np.sort(np.asarray(x, float))
    n = x.size
    if n < 10 or x[0] == x[-1]:
        return None
    res = diptest.diptest(x)
    d, pu = float(res[0]), float(res[1])  # type: ignore[arg-type]
    plf = float(1.0 - FM_pub(2.0 * np.sqrt(n) * d, clip_tails=True))
    return {"n": int(n), "dip": d, "p_uniform": pu, "p_LF": plf}


def main():
    from sklearn.datasets import load_digits
    from sklearn.neural_network import MLPClassifier
    X, y = load_digits(return_X_y=True)
    X = X / 16.0
    clf = MLPClassifier(hidden_layer_sizes=(64,), activation="relu", random_state=0, max_iter=500).fit(X, y)
    H = np.maximum(X @ clf.coefs_[0] + clf.intercepts_[0], 0.0)
    rows = []
    for u in range(H.shape[1]):
        h = H[:, u]
        rec = {"unit": u, "frac_zero": float(np.mean(h == 0)), "all": pvals(h), "pos": pvals(h[h > 0]),
               "pos_by_class": [pvals(h[(h > 0) & (y == c)]) for c in range(10)]}
        rows.append(rec)

    def count(key):
        cells = [r[key] for r in rows if r[key] is not None]
        return {"cells": len(cells), "uniform_p<.05": sum(c["p_uniform"] < 0.05 for c in cells),
                "LF_p<.05": sum(c["p_LF"] < 0.05 for c in cells)}
    cls = [c for r in rows for c in r["pos_by_class"] if c is not None]
    summary = {"train_acc": float(clf.score(X, y)), "all": count("all"), "pos": count("pos"),
               "pos_by_class": {"cells": len(cls), "uniform_p<.05": sum(c["p_uniform"] < 0.05 for c in cls),
                                "LF_p<.05": sum(c["p_LF"] < 0.05 for c in cls)},
               "dead_units": sum(r["frac_zero"] == 1.0 for r in rows),
               "median_frac_zero": float(np.median([r["frac_zero"] for r in rows]))}
    out = Path(__file__).resolve().parents[1] / "deposits" / "o5_activations.json"
    out.write_text(json.dumps({"summary": summary, "rows": rows}, indent=1))
    print(json.dumps(summary, indent=1))


if __name__ == "__main__":
    main()
