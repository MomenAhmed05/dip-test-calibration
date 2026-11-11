"""Step 16 rebuild helper (DT15-A): re-run the deterministic / seeded parts and write a deterministic record.

1. NE' Arb check (deterministic) on the frozen H_geom file, tag 'rebuild'; compare with `neprime_arb_confirm.json`
   field by field (timing field 'sec' and 'tag' excluded).
2. H_geom Monte Carlo re-run from its CONFIRM root (2026110848, N = 4e7, ~4 min); compare counts and bounds.
3. Paired-MC decision: re-run `dt16_pair_cert_summary.py confirm 2026110847 1e-3` on the frozen shards; keep stdout.
Writes deposits/dt16/rebuild_dt16.json (no timings; byte-stable). Temporary 'rebuild' outputs are deleted.
Usage (from the project root): python code/dt16_rebuild.py
"""
import json
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
DEP = os.path.join(HERE, "..", "deposits", "dt16")


def load(path):
    try:
        with open(path) as fh:
            return json.load(fh)
    except (OSError, ValueError) as err:
        sys.exit(f"cannot read {path}: {err}")


def run(args):
    r = subprocess.run([sys.executable, *args], cwd=HERE, capture_output=True, text=True, check=False)
    if r.returncode != 0:
        sys.exit(f"{args[0]} failed: {r.stderr[-400:]}")
    return r.stdout


def main():
    hg = os.path.join(DEP, "hgeom_confirm_r2026110848.json")
    run(["dt16_neprime_arb.py", hg, "rebuild"])
    a, b = load(os.path.join(DEP, "neprime_arb_confirm.json")), load(os.path.join(DEP, "neprime_arb_rebuild.json"))
    ne_match = all(a[k] == b[k] for k in a if k not in ("sec", "tag"))
    run(["dt16_hgeom_cert.py", "40000000", "2026110848", "rebuild", "1e-4"])
    c, d = load(hg), load(os.path.join(DEP, "hgeom_rebuild_r2026110848.json"))
    hg_match = all(c[k] == d[k] for k in ("counts", "bound", "cp_upper", "rs", "N", "K", "delta"))
    summary = run(["dt16_pair_cert_summary.py", "confirm", "2026110847", "1e-3"])
    for f in ("neprime_arb_rebuild.json", "hgeom_rebuild_r2026110848.json"):
        os.remove(os.path.join(DEP, f))
    rec = {"neprime_match_confirm": ne_match,
           "neprime_reach": [lv["reach"] for lv in a["levels"]],
           "neprime_analytic": [lv["analytic"] for lv in a["levels"]],
           "hgeom_match_confirm": hg_match,
           "pair_decisions": [line for line in summary.splitlines() if line.startswith(("COMPLETE", "DECISION"))],
           "pair_summary_stdout": summary}
    out = os.path.join(DEP, "rebuild_dt16.json")
    try:
        with open(out, "w") as fh:
            json.dump(rec, fh, indent=1, sort_keys=True)
    except OSError as err:
        sys.exit(f"cannot write {out}: {err}")
    print(json.dumps({k: rec[k] for k in ("neprime_match_confirm", "hgeom_match_confirm", "pair_decisions")}, indent=1))
    return 0 if ne_match and hg_match else 1


if __name__ == "__main__":
    sys.exit(main())
