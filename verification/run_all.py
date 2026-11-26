"""Run every verification script; logs in logs/, summary in logs/SUMMARY.txt."""
import os
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPTS = ["check_series.py", "check_wz.py", "check_limits.py",
           "check_completion.py", "check_campbell.py"]


def main():
    logdir = os.path.join(HERE, "logs")
    os.makedirs(logdir, exist_ok=True)
    summary, bad = [], 0
    for s in SCRIPTS:
        t0 = time.time()
        r = subprocess.run([sys.executable, s], cwd=HERE, capture_output=True, text=True)
        with open(os.path.join(logdir, s[:-3] + ".log"), "w") as fh:
            fh.write(r.stdout + ("\n--- stderr ---\n" + r.stderr if r.stderr else ""))
        last = [ln for ln in r.stdout.splitlines() if ln.startswith("--")]
        line = (last[-1] if last else f"-- {s}: crashed") + f"   [{time.time() - t0:.1f}s]"
        if r.returncode:
            bad += 1
            line += "   <-- FAILED"
        print(line, flush=True)
        summary.append(line)
    summary.append("ALL PASSED" if not bad else f"{bad} SCRIPT(S) FAILED")
    open(os.path.join(logdir, "SUMMARY.txt"), "w").write("\n".join(summary) + "\n")
    print(summary[-1])
    sys.exit(1 if bad else 0)


if __name__ == "__main__":
    main()
