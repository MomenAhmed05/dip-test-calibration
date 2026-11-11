"""Print compact tables from deposits/dt3c/confirm.json (Step 5 CONFIRM DT3c). Read-only.
Run: python code/dt3c_report.py
"""

import json
from pathlib import Path

P = Path(__file__).resolve().parents[1] / "deposits" / "dt3c" / "confirm.json"


def main():
    try:
        d = json.loads(P.read_text())
    except (OSError, json.JSONDecodeError) as e:
        print(f"cannot read {P}: {e} (has DT3c finished?)")
        return
    print("tier", d["tier"], "roots", d["roots"], "ref", d["ref_reps"], "law", d["law_reps"],
          "seconds", round(d.get("seconds", 0)))
    print("sqrt(n) q_n:", {n: {a: round(v, 4) for a, v in q.items()} for n, q in d["q_sqrtn"].items()})
    print("sup bound 1-F_M(2q_n):", {n: {a: round(100 * v, 2) for a, v in q.items()} for n, q in d["sup_bound"].items()})
    print()
    print(f"{'law':28s} {'n':>6s} | {'obs5':>6s} {'BF5':>6s} {'d5':>6s} | {'obs1':>6s} {'BF1':>6s} {'d1':>6s} | "
          f"{'maxtail':>7s} | {'q99|R|':>7s} {'P|R|>.01':>8s} {'P(R<0)':>7s}")
    for c in d["cells"]:
        a5, a1 = c["a0.05"], c["a0.01"]
        mt = max(abs(v["diff_pts"]) for v in c["tails"].values())
        print(f"{c['name'][:28]:28s} {c['n']:6d} | {100 * a5['obs']:6.2f} {100 * a5['r_BF']:6.2f} {a5['diff_pts']:+6.2f} | "
              f"{100 * a1['obs']:6.2f} {100 * a1['r_BF']:6.2f} {a1['diff_pts']:+6.2f} | {mt:7.2f} | "
              f"{c['q99_absR']:7.4f} {c['P(|R|>0.01)']:8.4f} {c['P(R<-1e-9)']:7.4f}")
    print()
    for k, v in d["readings"].items():
        print(k, "->", v["verdict"])


if __name__ == "__main__":
    main()
