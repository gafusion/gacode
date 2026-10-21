"""Independent mpmath reference for NEO's F(m,n,lambda), FB(m,n,lambda):

  F  = int_0^inf dx x^m e^{-x^2} int_0^x   dy y^n e^{-lambda y^2}
  FB = int_0^inf dx x^m e^{-x^2} int_x^inf dy y^n e^{-lambda y^2}

The outer integral is always over the narrower Maxwellian, in its own
units, so the integrand has no internal layer.
"""
from mpmath import mp, mpf, gammainc, quad, inf, exp

mp.dps = 50
PTS = [0, 1, 3, 6, inf]


def _I(a, z, upper):
    return gammainc(a, z) if upper else gammainc(a, 0, z)


def _ref(m, n, lam, upper_y):
    lam = mpf(lam)
    if lam <= 1:
        # outer x (narrow), inner y in closed form
        a = mpf(n+1)/2
        g = lambda x: x**m*exp(-x*x)*_I(a, lam*x*x, upper_y)/(2*lam**a)
        return quad(g, PTS, maxdegree=10)
    # outer u = sqrt(lam) y (narrow), inner x in closed form
    a = mpf(m+1)/2
    g = lambda u: u**n*exp(-u*u)*_I(a, u*u/lam, not upper_y)/2
    return quad(g, PTS, maxdegree=10)/lam**(mpf(n+1)/2)


def F(m, n, lam):
    return _ref(m, n, lam, upper_y=False)


def FB(m, n, lam):
    return _ref(m, n, lam, upper_y=True)
