"""emat_coll_test or emat_coll_field (COL=1,2 via env), 3 species (D, e, W): max over (is,js,ix) blocks of
max|error|/max|block|, with entries scaled by 1/sqrt(h_i h_j).

usage: python3 check_coll.py n_xi laguerre_method collision_model ref n_energy...
       ref = old16 or new16
"""
import os, subprocess, sys
COL = int(os.environ.get('COL', '1'))
UNC = os.environ.get('UNC', '0')
from math import lgamma, exp, sqrt

nx, meth, model, ref = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]), sys.argv[4]
alpha = {1: [0.5] + [1.5]*nx, 2: [i+0.5 for i in range(nx+1)],
         3: [0.5]*(nx+1), 4: [1.5]*(nx+1)}[meth]
h = lambda n, a: 0.5*exp(lgamma(n+a+1)-lgamma(n+1))
names = {1: 'D', 2: 'e', 3: 'W'}


def load(exe, ne):
    out = subprocess.run([f"./c{exe}", str(ne), str(nx), str(meth), str(model), UNC],
                         capture_output=True, text=True).stdout
    R = {}
    for line in out.splitlines():
        t = line.split()
        is_, js, ix, ie, je = map(int, t[:5])
        R[(is_, js, ix, ie, je)] = float(t[4+COL])/sqrt(h(ie, alpha[ix])*h(je, alpha[ix]))
    return R


def worst(A, B, ne):
    res = {}
    for is_ in (1, 2, 3):
        for js in (1, 2, 3):
            w = 0.0
            for ix in range(nx+1):
                blk = [(ix, i, j) for i in range(ne+1) for j in range(ne+1)]
                sc = max(abs(B[(is_, js) + k]) for k in blk)
                if sc == 0: continue
                w = max(w, max(abs(A[(is_, js) + k]-B[(is_, js) + k]) for k in blk)/sc)
            res[names[is_]+names[js]] = w
    return res


pairs = [a+b for a in 'DeW' for b in 'DeW']
print(f"{['test','field'][COL-1]} model={model} method={meth} n_xi={nx} uncoupled={UNC}: vs {ref}")
print(f"{'n_e':>4} {'code':>6} " + " ".join(f"{p:>7}" for p in pairs))
for ne in map(int, sys.argv[5:]):
    R = load(ref, ne)
    for exe in (("old8", "new8") if ref == "old16" else ("new8",)):
        w = worst(load(exe, ne), R, ne)
        print(f"{ne:4d} {exe:>6} " + " ".join(f"{w[p]:7.0e}" for p in pairs))
