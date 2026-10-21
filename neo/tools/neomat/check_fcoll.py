"""Relative error of NEO's F, FB (real*8, real*16) vs mpmath, for one lambda."""
import subprocess, sys
from mpmath import mpf
from fref import F, FB

m0, lam = int(sys.argv[1]), sys.argv[2]
dump = lambda exe: {(int(t[0]), int(t[1])): (mpf(t[2]), mpf(t[3])) for t in
        (l.split() for l in subprocess.run([f"./{exe}", str(m0), lam],
         capture_output=True, text=True).stdout.splitlines())}
d8, d16 = dump("fcoll8"), dump("fcoll16")
worst = {}
for (m, n) in d8:
    if (m + n) % 2 == 0 or m + n < 0 or m < -m0+2 or n < -m0+2: continue  # NEO fills only m+n odd
    # NEO uses F only with n >= 0 and FB only with m >= 0
    for k, ok, fn in ((0, n >= 0, F), (1, m >= 0, FB)):
        if not ok: continue
        ref = fn(m, n, lam)
        for tag, d in (("8", d8), ("16", d16)):
            e = float(abs(d[(m, n)][k]/ref - 1))
            key = ("F", "FB")[k] + tag
            if e > worst.get(key, (0,))[0]: worst[key] = (e, m, n)
for k in sorted(worst): print(f"  {k:5s} max rel err {worst[k][0]:.1e} at (m,n)={worst[k][1:]}")
