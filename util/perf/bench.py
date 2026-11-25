#!/usr/bin/env python3
"""Benchmark openMalaria builds: CPU time (median of N reps, interleaved), peak RSS, malloc count.

usage: bench.py <repo> <reps> <label>=<openMalaria> [<label>=<openMalaria> ...]
Reps are interleaved across binaries (A B A B ...) so drift affects all equally.
Metric is user+sys CPU seconds (BENCH_WALL=1 for wall time); BENCH_SCENARIOS="A B" overrides the scenario list.
Correctness: golden.py runs every regression scenario and compares outputs byte for byte.
"""
import os, re, shutil, statistics, subprocess, sys, tempfile, time

SCENARIOS = os.environ["BENCH_SCENARIOS"].split() if "BENCH_SCENARIOS" in os.environ else ["Genotypes", "Penny", "MSAT_HRP2", "TBV", "VecFullTest", "NonHumanHosts", "Molineaux"]
CPU = os.environ.get("BENCH_WALL") is None  # default metric: user+sys CPU seconds (less scheduler noise than wall)
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from golden import schema_of

def run(binary, repo, name, count_mallocs=False):
    testdir = os.path.join(repo, "test")
    xml = os.path.join(testdir, "scenario%s.xml" % name)
    sim = tempfile.mkdtemp()
    xsd = schema_of(xml)
    src = [os.path.join(d, xsd) for d in (os.path.join(repo, "schema"), os.path.join(os.path.dirname(binary), "schema"))]
    shutil.copy(next(p for p in src if os.path.isfile(p)), sim)
    env = dict(os.environ)
    if count_mallocs:
        env["LD_PRELOAD"] = os.path.join(HERE, "malloccount.so")
    t = time.perf_counter()
    p = subprocess.Popen([binary, "--resource-path", testdir, "--scenario", xml], cwd=sim, env=env,
                         stdout=subprocess.DEVNULL, stderr=subprocess.PIPE)
    _, status, ru = os.wait4(p.pid, 0)
    dt = time.perf_counter() - t
    if CPU: dt = ru.ru_utime + ru.ru_stime
    err = p.stderr.read().decode()
    shutil.rmtree(sim)
    assert status == 0, "%s %s failed: %s" % (binary, name, err[-500:])
    m = re.search(r"MALLOC_COUNT (\d+)", err)
    return dt, ru.ru_maxrss, int(m.group(1)) if m else None

def main():
    repo, reps = os.path.abspath(sys.argv[1]), int(sys.argv[2])
    bins = [a.split("=", 1) for a in sys.argv[3:]]
    times = {(l, s): [] for l, _ in bins for s in SCENARIOS}
    rss, mallocs = {}, {}
    for l, b in bins:
        for s in SCENARIOS:
            _, rss[l, s], mallocs[l, s] = run(b, repo, s, count_mallocs=True)
    for _ in range(reps):
        for s in SCENARIOS:
            for l, b in bins:
                times[l, s].append(run(b, repo, s)[0])
    base = bins[0][0]
    print("%-14s %-14s %9s %7s %8s %8s %14s" % ("scenario", "build", "median_s", "spread", "vs_base", "rss_MB", "mallocs"))
    tot = {l: 0.0 for l, _ in bins}
    for s in SCENARIOS:
        for l, _ in bins:
            t = times[l, s]; med = statistics.median(t); tot[l] += med
            b = statistics.median(times[base, s])
            print("%-14s %-14s %9.3f %6.1f%% %+7.1f%% %8.1f %14s" % (s, l, med, 100 * (max(t) - min(t)) / med,
                  100 * (med / b - 1), rss[l, s] / 1024, mallocs[l, s]))
    print()
    for l, _ in bins:
        print("TOTAL %-14s %8.2fs  %+6.1f%% vs %s" % (l, tot[l], 100 * (tot[l] / tot[base] - 1), base))

if __name__ == "__main__":
    main()
