#!/usr/bin/env python
# Plot one moment, [p,e,m]flux(theta) in friction form, with one panel per species in a
# single window, from out.neo.flux_theta (3 columns: pflux, eflux, mflux). The FSA values
# from out.neo.transport (and +gv from out.neo.transport_gv) are quoted on each panel.
# With norm=1, each flux is divided by its GB value and by n_a/n_e (per particle).
# Usage: neo_plot_flux_th.py <simdir> <ftype> [species: all or 1,3,...] [moment: pflux|eflux|mflux] [norm: 0|1]
import os
import sys

import numpy as np
import matplotlib
if sys.argv[2] != 'screen':
    matplotlib.use('Agg')
import matplotlib.pyplot as plt

simdir, ftype = sys.argv[1], sys.argv[2]
spec = sys.argv[3] if len(sys.argv) > 3 else 'all'
moment = sys.argv[4] if len(sys.argv) > 4 else 'pflux'
norm = (len(sys.argv) > 5 and sys.argv[5] in ('1', 'norm'))
if moment not in ('pflux', 'eflux', 'mflux'):
    print('ERROR: moment must be pflux, eflux or mflux (got %s)' % moment)
    sys.exit(1)

fname = None
for name in ('out.neo.flux_theta',):
    if os.path.isfile(simdir + '/' + name):
        fname = name
        break
if fname is None:
    print('out.neo.flux_theta not found; nothing to plot.')
    sys.exit(0)

# Grid: n_species, n_energy, n_xi, n_theta, theta(1:n_theta), n_radial, r(:)
tok = open(simdir + '/out.neo.grid').read().split()
n_species = int(tok[0])
n_theta = int(tok[3])
theta = np.array(tok[4:4 + n_theta], dtype=float)
n_radial = int(tok[4 + n_theta])

ir = 0  # radial index to plot (0-based), relevant only if n_radial > 1

# Species-major within each radius. Columns (friction form):
# pflux_th eflux_th mflux_th
raw = np.loadtxt(simdir + '/' + fname)
ncol = 3
data = raw.reshape(n_radial, n_species, n_theta, ncol)[ir]

# FSA reference values, one row per radius.
# out.neo.transport: 5 leading columns then 8 per species (pflux, eflux, mflux, ...)
# out.neo.transport_gv: r, then 3 per species (pflux_gv, eflux_gv, mflux_gv)
# The friction-form average equals the DKE flux + gyroviscous flux.
tr = np.atleast_2d(np.loadtxt(simdir + '/out.neo.transport'))[ir]
gv = np.atleast_2d(np.loadtxt(simdir + '/out.neo.transport_gv'))[ir]
fsa = [tr[5::8], tr[6::8], tr[7::8]]
fsa_gv = [gv[1::3], gv[2::3], gv[3::3]]
names = ['pflux', 'eflux', 'mflux']

# out.neo.diagnostic_geo: '#' header lines, then theta, v_drift_x, gradpar_Bmag,
# Bmag, w_theta, R (n_theta each), then R(theta=0) and dR(theta=0)/dr
geo = np.loadtxt(simdir + '/out.neo.diagnostic_geo', comments='#')
w_theta = geo[4 * n_theta:5 * n_theta]

# Total flux = sum_theta w_theta * flux(theta); compare with the NEO output
total = np.einsum('t,stk->sk', w_theta, data)

# --- Optional normalization: GB units and per-particle (divide by n_a/n_e) ---
if norm:
    equil = np.atleast_2d(np.loadtxt(simdir + '/out.neo.equil'))[ir]
    rho_s = equil[3]
    dens_a = equil[7::5][:n_species]
    temp_a = equil[8::5][:n_species]
    zsp = np.atleast_1d(np.loadtxt(simdir + '/out.neo.species'))[1::2]
    if np.any(zsp < 0):
        # kinetic electrons: take n_e, T_e directly
        ie_ = int(np.argmax(zsp < 0))
        n_e, T_e = dens_a[ie_], temp_a[ie_]
    else:
        # adiabatic electrons: infer n_e, T_e from the GB-normalized fluxes in
        # out.neo.transport_flux (pgb = n_e rho^2 T_e^1.5, egb = n_e rho^2 T_e^2.5)
        rows = [l.split() for l in open(simdir + '/out.neo.transport_flux') if l.strip() and l[0] != '#']
        gbv = np.array(rows[2 * ir * n_species:(2 * ir + 1) * n_species], dtype=float)  # dke block (gv block follows)
        js = int(np.argmax(np.abs(tr[5::8])))
        pgb_ = tr[5::8][js] / gbv[js, 1]
        egb_ = tr[6::8][js] / gbv[js, 2]
        T_e = egb_ / pgb_
        n_e = pgb_ / (rho_s**2 * T_e**1.5)
    gb = [n_e * rho_s**2 * T_e**p for p in (1.5, 2.5, 2.0)]   # pflux, eflux, mflux
    div = np.ones((n_species, ncol))
    for ks in range(ncol):
        div[:, ks] = gb[ks] * dens_a / n_e
    data = data / div[:, None, :]
    total = total / div
    fsa = [fsa[ks] / div[:, ks] for ks in range(ncol)]
    fsa_gv = [fsa_gv[ks] / div[:, ks] for ks in range(ncol)]
    nlab = ' / (GB n_a/n_e)'  # shown in the title
else:
    nlab = ''

species = range(1, n_species + 1) if spec == 'all' else [int(s) for s in spec.split(',')]

k = names.index(moment)

# One panel per species, all in one window
fig, axs = plt.subplots(len(species), 1, figsize=(7, 3.0 * len(species) + 0.8),
                        sharex=True, squeeze=False)
for ax, is_ in zip(axs[:, 0], species):
    i = is_ - 1
    ax.plot(theta / np.pi, data[i, :, k], '-o', ms=3)
    ax.set_ylabel('species %d' % is_)
    ax.grid(alpha=0.3)
    txt = 'sum(w*%s_th) = %.6e\nNEO %s (+gv) = %.6e (%.6e)' % (
        moment, total[i, k], moment, fsa[k][i], fsa[k][i] + fsa_gv[k][i])
    ax.text(0.03, 0.97, txt, transform=ax.transAxes, va='top', family='monospace', fontsize=8,
            bbox=dict(boxstyle='round', fc='white', ec='gray', alpha=0.85))
axs[-1, 0].set_xlabel(r'$\theta/\pi$')
fig.suptitle(moment + nlab.replace(' / ', '/') if norm else moment)
fig.tight_layout()
if ftype != 'screen':
    fig.savefig('%s/out.neo.flux_theta.%s.%s' % (simdir, moment, ftype), dpi=150)
else:
    plt.show()
