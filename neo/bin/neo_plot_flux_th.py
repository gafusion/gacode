#!/usr/bin/env python
# Plot [p,e,m]flux(theta) and [p,e,m]flux_gv(theta) from out.neo.flux_theta
# and print the FSA values from out.neo.transport and out.neo.transport_gv
# Usage: neo_plot_pflux_th.py <simdir> <ftype> [species: all or 1,3,...]
import os
import sys

import numpy as np
import matplotlib
if sys.argv[2] != 'screen':
    matplotlib.use('Agg')
import matplotlib.pyplot as plt

simdir, ftype = sys.argv[1], sys.argv[2]
spec = sys.argv[3] if len(sys.argv) > 3 else 'all'

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
# pflux_th pflux_gv_th eflux_th eflux_gv_th mflux_th mflux_gv_th
data = np.loadtxt(simdir + '/' + fname).reshape(n_radial, n_species, n_theta, 6)[ir]

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
total = np.einsum('t,stk->sk', w_theta, data)

species = range(1, n_species + 1) if spec == 'all' else [int(s) for s in spec.split(',')]

#print('%-3s %-9s %16s %16s %10s' % ('is', 'flux', 'sum(w*flux_th)', 'NEO output', 'rel diff'))
#for is_ in species:
#    for k in range(6):
#        ref = fsa[k][is_ - 1]
#        print('%-3d %-9s %16.8e %16.8e %10.2e' % (
#            is_, names[k], total[is_ - 1, k], ref,
#            abs(total[is_ - 1, k] - ref) / max(abs(ref), 1e-300)))

for is_ in species:
    i = is_ - 1
    fig, axs = plt.subplots(3, 2, figsize=(11, 9), sharex=True)
    for k in range(6):
        ax = axs[k // 2][k % 2]
        ax.plot(theta, data[i, :, k], '-o', ms=3)
        ax.set_ylabel(names[k] + '_th')
        ax.grid(alpha=0.3)
        ax.text(0.03, 0.97, 'sum(w*%s_th) = %.6e\n%s: %s = %.6e'
                % (names[k], total[i, k], files[k], names[k], fsa[k][i]),
                transform=ax.transAxes, va='top', family='monospace', fontsize=9,
                bbox=dict(boxstyle='round', fc='white', ec='gray'))
    axs[2][0].set_xlabel(r'$\theta$')
    axs[2][1].set_xlabel(r'$\theta$')
    fig.suptitle('species %d' % is_)
    fig.tight_layout()
    if ftype != 'screen':
        fig.savefig('%s/out.neo.flux_theta.s%d.%s' % (simdir, is_, ftype), dpi=150)

if ftype == 'screen':
    plt.show()
