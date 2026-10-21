"""ENERGY_basis_ints arrays: max |error| normalized by sqrt(h_i h_j),
h_n = Gamma(n+alpha+1)/(2 n!).  Reference: old code (HEAD) in real*16.

usage: python3 check_basis.py n_xi laguerre_method n_energy [n_energy ...]
"""
import subprocess, sys
from math import lgamma, exp, sqrt

nx, meth = int(sys.argv[1]), int(sys.argv[2])
REF = "old16"
if sys.argv[3] == "--new":    # reference: new code in real*16
    REF = "new16"; del sys.argv[3]
alpha = [0.5] + [1.5]*nx if meth == 1 else [i+0.5 for i in range(nx+1)]
h = lambda n, a: 0.5*exp(lgamma(n+a+1)-lgamma(n+1))


def load(exe, ne):
    R = {}
    out = subprocess.run([f"./{exe}", str(ne), str(nx), str(meth)],
                         capture_output=True, text=True).stdout
    for line in out.splitlines():
        t = line.split()
        tag = t[0]
        if tag in ('v', 'e0'):
            ix, ie, je = map(int, t[1:4]); vals = t[4:]; kx = 0
        else:
            ix, ie, je, kx = map(int, t[1:5]); vals = t[5:]
        names = {'v': ['e0v', 'e05v', 'e1v', 'e105v', 'e2v'], 'e0': ['e0'],
                 'k2': ['e05', 'en05', 'e05de'], 'e1': ['e1']}[tag]
        for nm, v in zip(names, vals):
            R[(nm, ix, ie, je, kx)] = float(v)
    return R


def jx_of(nm, ix, kx):
    if nm in ('e05', 'en05', 'e05de'): return ix+2*kx-3
    if nm == 'e1': return ix+2*kx-4
    return ix


def err(A, B):
    worst = {}
    for k, b in B.items():
        nm, ix, ie, je, kx = k
        if nm.endswith('v'):
            sc = sqrt(h(ie, alpha[ix]))
        else:
            jx = jx_of(nm, ix, kx)
            if jx < 0 or jx > nx: continue
            sc = sqrt(h(ie, alpha[ix])*h(je, alpha[jx]))
        e = abs(A[k]-b)/sc
        worst[nm] = max(worst.get(nm, 0.0), e)
    return worst


names = ['e0', 'e05', 'en05', 'e05de', 'e1', 'e0v', 'e05v', 'e1v', 'e105v', 'e2v']
print(f"n_xi={nx} method={meth}: max normalized error vs {REF}")
print(f"{'n_e':>4} {'code':>6} " + " ".join(f"{n:>8}" for n in names))
for ne in map(int, sys.argv[3:]):
    ref = load(REF, ne)
    for exe in (("old8", "new8") if REF == "old16" else ("new8",)):
        w = err(load(exe, ne), ref)
        print(f"{ne:4d} {exe:>6} " + " ".join(f"{w[n]:8.1e}" for n in names))
