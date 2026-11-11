"""DT9 orchestrator (DT9-PREREG): runs the frozen stages and applies the frozen stop rule.

Stage l = 1, 2, 3: every shard s = 0..10 runs, one after the other, n = 800 then n = 1600 (only the n not yet
stopped) up to the look-l chunk count (dt9_block_cert.py resumes: existing chunks are skipped). After the
stage, dt9_cert_analysis.look() is evaluated; an n stops once 0.7315 is certified at alpha = 0.05.
Usage: python dt9_orchestrate.py root
"""
import os
import subprocess
import sys
import threading
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from dt9_cert_analysis import LOOKS, NSH, look  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
LOGDIR = os.path.join(HERE, "..", "deposits", "dt9")


def log(msg, root):
    line = f"{time.strftime('%Y-%m-%d %H:%M:%S', time.gmtime())}Z {msg}"
    print(line, flush=True)
    try:
        with open(os.path.join(LOGDIR, f"orchestrate_r{root}.log"), "a") as fh:
            fh.write(line + "\n")
    except OSError as err:
        sys.exit(f"cannot write log: {err}")


def shard_chain(s, jobs, root, errs):
    for n, L in jobs:
        cmd = [sys.executable, os.path.join(HERE, "dt9_block_cert.py"), str(n), str(L), str(s), str(NSH),
               str(root), "cert", "500"]
        try:
            with open(os.path.join(LOGDIR, f"cert_n{n}_r{root}_s{s}.log"), "a") as fh:
                rc = subprocess.run(cmd, stdout=fh, stderr=subprocess.STDOUT, check=False).returncode
        except OSError as err:
            errs.append((s, n, str(err)))
            return
        if rc != 0:
            errs.append((s, n, rc))
            return


def main():
    root = int(sys.argv[1])
    stopped = dict.fromkeys(LOOKS, False)
    for li in range(3):
        jobs = [(n, LOOKS[n][li]) for n in LOOKS if not stopped[n]]
        if not jobs:
            break
        log(f"stage {li + 1} start jobs={jobs}", root)
        errs = []
        th = [threading.Thread(target=shard_chain, args=(s, jobs, root, errs)) for s in range(NSH)]
        for t in th:
            t.start()
        for t in th:
            t.join()
        if errs:
            log(f"stage {li + 1} ERRORS {errs}; stop (resume by rerunning this script with the same root)", root)
            sys.exit(1)
        for n, _ in jobs:
            r = look(root, n, li)
            if r is None:
                log(f"n={n} look {li + 1}: chunks missing; stop", root)
                sys.exit(1)
            a5, a1 = r["alpha"]["0.05"], r["alpha"]["0.01"]
            log(f"n={n} look {li + 1} N={r['N']} 5% cert={a5['certified']} K={a5['trail'][0]['K']} b={a5['b']} "
                f"| 1% cert={a1['certified']} b={a1['b']}", root)
            if a5["certified"] == 0.7315:
                stopped[n] = True
    log(f"done stopped={stopped}", root)


if __name__ == "__main__":
    main()
