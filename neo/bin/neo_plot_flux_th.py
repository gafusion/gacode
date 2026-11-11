#!/usr/bin/env python
# Plot one moment, [p,e,m]flux(theta) or [p,e,m]flux_gv(theta), with one panel per
# species in a single window, from out.neo.flux_theta; the FSA value from
# out.neo.transport (or out.neo.transport_gv) is shown on each panel. For pflux,
# the friction-form pflux(theta) (7th column) is overlaid if present.
# With norm=1, each species' flux is divided by its GB value (as in out.neo.transport_flux)
# and by n_a/n_e, i.e. plotted per particle: Gamma_a/(Gamma_GB n_a/n_e).
# Usage: neo_plot_flux_th.py <simdir> <ftype> [species: all or 1,3,...] [moment: pflux|eflux|mflux|pflux_gv|eflux_gv|mflux_gv] [norm: 0|1]
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
if moment not in ('pflux', 'eflux', 'mflux', 'pflux_gv', 'eflux_gv', 'mflux_gv'):
    print('ERROR: moment must be pflux, eflux, mflux, pflux_gv, eflux_gv or mflux_gv (got %s)' % moment)
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

# Species-major within each radius. Columns:
# pflux_th pflux_gv_th eflux_th eflux_gv_th mflux_th mflux_gv_th [pflux_fric_th]
# (the 7th column, friction-form pflux(theta), is absent in older output)
raw = np.loadtxt(simdir + '/' + fname)
ncol = raw.size // (n_radial * n_species * n_theta)
data = raw.reshape(n_radial, n_species, n_theta, ncol)[ir]

# FSA values, one row per radius.
# out.neo.transport: 5 leading columns then 8 per species (pflux, eflux, mflux, ...)
# out.neo.transport_gv: r, then 3 per species (pflux_gv, eflux_gv, mflux_gv)
tr = np.atleast_2d(np.loadtxt(simdir + '/out.neo.transport'))[ir]
gv = np.atleast_2d(np.loadtxt(simdir + '/out.neo.transport_gv'))[ir]
fsa = [tr[5::8], gv[1::3], tr[6::8], gv[2::3], tr[7::8], gv[3::3]]
names = ['pflux', 'pflux_gv', 'eflux', 'eflux_gv', 'mflux', 'mflux_gv']
files = ['out.neo.transport', 'out.neo.transport_gv'] * 3

# out.neo.diagnostic_geo: '#' header lines, then theta, v_drift_x, gradpar_Bmag,
# Bmag, w_theta, R (n_theta each), then R(theta=0) and dR(theta=0)/dr
geo = np.loadtxt(simdir + '/out.neo.diagnostic_geo', comments='#')
w_theta = geo[4 * n_theta:5 * n_theta]

# Total flux = sum_theta w_theta * flux(theta); compare with the NEO output
total = np.einsum('t,stk->sk', w_theta, data)  # includes fric column if present

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
    # divisor for each species and column of data / fsa
    div = np.ones((n_species, ncol))
    for ks in range(ncol):
        mm = 0 if ks == 6 else ks // 2
        div[:, ks] = gb[mm] * dens_a / n_e
    data = data / div[:, None, :]
    total = total / div
    fsa = [fsa[ks] / div[:, ks] for ks in range(6)]
    nlab = ' / (GB n_a/n_e)'  # shown in the title
else:
    nlab = ''

species = range(1, n_species + 1) if spec == 'all' else [int(s) for s in spec.split(',')]

#print('%-3s %-9s %16s %16s %10s' % ('is', 'flux', 'sum(w*flux_th)', 'NEO output', 'rel diff'))
#for is_ in species:
#    for k in range(6):
#        ref = fsa[k][is_ - 1]
#        print('%-3d %-9s %16.8e %16.8e %10.2e' % (
#            is_, names[k], total[is_ - 1, k], ref,
#            abs(total[is_ - 1, k] - ref) / max(abs(ref), 1e-300)))

allm = ['pflux', 'pflux_gv', 'eflux', 'eflux_gv', 'mflux', 'mflux_gv']
k = allm.index(moment)           # column of this moment in data / total / fsa
is_gv = moment.endswith('_gv')
show_fric = (k == 0 and ncol > 6)

# One panel per species, all in one window
fig, axs = plt.subplots(len(species), 1, figsize=(7, 3.0 * len(species) + 0.8),
                        sharex=True, squeeze=False)
for ax, is_ in zip(axs[:, 0], species):
    i = is_ - 1
    ax.plot(theta / np.pi, data[i, :, k], '-o', ms=3,
            label=moment + '_th' + (' (drift form)' if show_fric else ''))
    if show_fric:
        ax.plot(theta / np.pi, data[i, :, 6], '--s', ms=3, label='pflux_th (friction form)')
        ax.legend(loc='lower left', bbox_to_anchor=(0.0, 1.0), ncol=2, fontsize=8, frameon=False)
    ax.set_ylabel('species %d' % is_)
    ax.grid(alpha=0.3)
    txt = 'sum(w*%s_th) = %.6e\n%s: %s = %.6e' % (names[k], total[i, k], files[k], names[k], fsa[k][i])
    if show_fric:
        txt += '\nsum(w*pflux_fric_th) = %.6e' % total[i, 6]
    ax.text(0.03, 0.97, txt, transform=ax.transAxes, va='top', family='monospace', fontsize=8,
            bbox=dict(boxstyle='round', fc='white', ec='gray', alpha=0.85))
axs[-1, 0].set_xlabel(r'$\theta/\pi$')
fig.suptitle(moment + nlab.replace(' / ', '/') if norm else moment)
fig.tight_layout()
if ftype != 'screen':
    fig.savefig('%s/out.neo.flux_theta.%s.%s' % (simdir, moment, ftype), dpi=150)
else:
    plt.show()
