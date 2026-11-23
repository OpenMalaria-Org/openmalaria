#!/usr/bin/env python3
"""Run every regression scenario with a binary and store all output files.

usage: golden.py run <openMalaria> <repo> <outdir> [-j N]
       golden.py cmp <dirA> <dirB>      (byte-for-byte; exit 1 on any difference)
"""
import os, re, shutil, subprocess, sys, tempfile, time, filecmp
from concurrent.futures import ThreadPoolExecutor

def scenario_names(repo):
    txt = open(os.path.join(repo, "test/CMakeLists.txt")).read()
    block = re.search(r"set \(OM_BOXTEST_NAMES(.*?)\n\)", txt, re.S).group(1)
    return [n for l in block.splitlines() for n in l.split("#")[0].split()]

def schema_of(xml):
    head = open(xml).read(4000)
    m = re.search(r'noNamespaceSchemaLocation="([^"]+)"', head) or re.search(r'schemaLocation="\S+ ([^"]+)"', head)
    return m.group(1)

def run_one(binary, repo, outdir, name):
    # scenario_current.xsd is generated next to the binary; others live in repo/schema
    schemas = [os.path.join(repo, "schema"), os.path.join(os.path.dirname(binary), "schema")]
    testdir = os.path.join(repo, "test")
    xml = os.path.join(testdir, "scenario%s.xml" % name)
    sim = tempfile.mkdtemp(prefix=name + "-", dir=outdir)
    xsd = schema_of(xml)
    shutil.copy(next(os.path.join(d, xsd) for d in schemas if os.path.isfile(os.path.join(d, xsd))), sim)
    t = time.time()
    r = subprocess.run([binary, "--resource-path", testdir, "--scenario", xml],
                       cwd=sim, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE)
    dt = time.time() - t
    dest = os.path.join(outdir, name)
    os.makedirs(dest)
    for f in os.listdir(sim):
        if not f.endswith(".xsd"):
            shutil.move(os.path.join(sim, f), dest)
    shutil.rmtree(sim)
    open(os.path.join(dest, "_exit"), "w").write("%d\n" % r.returncode)
    return name, r.returncode, dt

def run(binary, repo, outdir, jobs):
    os.makedirs(outdir)
    names = scenario_names(repo)
    with ThreadPoolExecutor(jobs) as ex:
        res = list(ex.map(lambda n: run_one(os.path.abspath(binary), repo, outdir, n), names))
    for n, rc, dt in sorted(res, key=lambda r: -r[2]):
        print("%-32s rc=%d %7.2fs" % (n, rc, dt))

def cmp(a, b):
    bad = []
    for name in sorted(os.listdir(a)):
        fa, fb = set(os.listdir(os.path.join(a, name))), set(os.listdir(os.path.join(b, name))) if os.path.isdir(os.path.join(b, name)) else set()
        if fa != fb:
            bad.append("%s: file sets differ %s vs %s" % (name, sorted(fa), sorted(fb)))
            continue
        for f in sorted(fa):
            if not filecmp.cmp(os.path.join(a, name, f), os.path.join(b, name, f), shallow=False):
                bad.append("%s/%s differs" % (name, f))
    print("\n".join(bad) if bad else "IDENTICAL: all %d scenarios" % len(os.listdir(a)))
    return 1 if bad else 0

if __name__ == "__main__":
    if sys.argv[1] == "run":
        jobs = int(sys.argv[sys.argv.index("-j") + 1]) if "-j" in sys.argv else 1
        run(sys.argv[2], os.path.abspath(sys.argv[3]), sys.argv[4], jobs)
    else:
        sys.exit(cmp(sys.argv[2], sys.argv[3]))
