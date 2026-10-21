"""NEO collision_model=4 matrices for one species pair: real*8 and real*16
(drv8, drv16, using neo_compute_fcoll) vs mpmath (same formulas, F from fref).

usage: python3 check_matrix.py n_energy laguerre_method n_xi m_a/m_b T_a/T_b
"""
import os, pickle, subprocess, sys
from mpmath import mp, mpf, gamma, sqrt, pi
from fref import F, FB

ne, meth, nx = map(int, sys.argv[1:4])
mr, tr = mpf(sys.argv[4]), mpf(sys.argv[5])
lam = tr/mr

cache_file = f"fcache_{float(lam):.6e}.pkl"
cache = pickle.load(open(cache_file, "rb")) if os.path.exists(cache_file) else {}


def fc(m, n, bar=False):
    k = (m, n, bar)
    if k not in cache:
        cache[k] = (FB if bar else F)(m, n, lam)
    return cache[k]


g2 = lambda n: gamma(mpf(n)/2)
e_lag = [1] + [3]*nx if meth == 1 else [2*i+1 for i in range(nx+1)]
xb = [0] + [1]*nx if meth == 1 else list(range(nx+1))


def ref(l):
    E = [[mpf(0)]*(ne+1) for _ in range(ne+1)]
    T = [[mpf(0)]*(ne+1) for _ in range(ne+1)]
    Fd = [[mpf(0)]*(ne+1) for _ in range(ne+1)]
    b = xb[l]
    for i in range(ne+1):
        for j in range(ne+1):
            for k in range(i+1):
                for m in range(j+1):
                    z = (-1)**(k+m)*(g2(2+2*i+e_lag[l])/g2(2+2*(i-k))/g2(2+2*k+e_lag[l])) \
                        * (g2(2+2*j+e_lag[l])/g2(2+2*(j-m))/g2(2+2*m+e_lag[l])) \
                        / (g2(2+2*k)*g2(2+2*m))
                    p = 2*(k+m) + 2*b
                    E[i][j] += z*g2(p+3)/2
                    T[i][j] += -sqrt(lam/pi)*l*(l+1)*z*(fc(p-1, 0) - fc(p-3, 2))
                    T[i][j] += 4*sqrt(lam/pi)*(2*k+b)*z*((1-tr)*fc(p-1, 2) - (2*m+b)*fc(p-3, 2)/2)
                    r1 = mr*(1+lam)**(-mpf(p+3)/2)*g2(p+3)
                    xa, ya = 2*k+b, 2*m+b
                    r2 = -2/(l+mpf(1)/2)*(mr-l*(1-mr))*fc(xa-l+1, ya+l+2)
                    r3 = -2/(l+mpf(1)/2)*(1+l*(1-mr))*fc(xa+l+2, ya-l+1, True)
                    w45 = -l*(l-1)/(l*l-mpf(1)/4)
                    w67 = (l+1)*(l+2)/((l+mpf(3)/2)*(l+mpf(1)/2))
                    r45 = w45*(fc(xa-l+3, ya+l+2) + fc(xa+l+2, ya-l+3, True))
                    r67 = w67*(fc(xa-l+1, ya+l+4) + fc(xa+l+4, ya-l+1, True))
                    Fd[i][j] += 2/sqrt(pi)*lam**(mpf(3)/2)*lam**(mpf(2*m+b)/2)*z*(r1+r2+r3+r45+r67)
    return E, T, Fd


def run(exe):
    R = {}
    out = subprocess.run([f"./{exe}", str(ne), str(meth), str(nx), sys.argv[4], sys.argv[5]],
                         capture_output=True, text=True).stdout
    for line in out.splitlines():
        t = line.split()
        R[tuple(map(int, t[:3]))] = [mpf(v) for v in t[3:]]
    return R


R8, R16 = run("drv8"), run("drv16")
print(f"lambda={float(lam):.4g}  n_energy={ne}  method={meth}")
print("   l | |dC_test|/|C_test|  real*8  real*16 | |dC_field|/|C_field|  real*8  real*16 | |C_field|/|C_test|")
for l in range(nx+1):
    E, T, Fd = ref(l)
    d = [1/sqrt(E[i][i]) for i in range(ne+1)]
    def err(R, col, A):
        sc = max(abs(A[i][j]*d[i]*d[j]) for i in range(ne+1) for j in range(ne+1))
        e = max(abs((R[(l, i, j)][col]-A[i][j])*d[i]*d[j]) for i in range(ne+1) for j in range(ne+1))
        return float(e/sc) if sc else float(e), sc
    t8, st = err(R8, 1, T); t16, _ = err(R16, 1, T)
    f8, sf = err(R8, 2, Fd); f16, _ = err(R16, 2, Fd)
    print(f"  {l:2d} |                    {t8:7.1e}  {t16:7.1e} |"
          f"                       {f8:7.1e}  {f16:7.1e} | {float(sf/st) if st else float('nan'):8.1e}")
pickle.dump(cache, open(cache_file, "wb"))
