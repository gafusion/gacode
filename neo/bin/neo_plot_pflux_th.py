#!/usr/bin/env python
# Plot pflux_th and pflux_gv_th vs. theta (out.neo.pflux_th) and print the
# values from out.neo.transport and out.neo.transport_gv on the plot.
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

if not os.path.isfile(simdir + '/out.neo.pflux_th'):
    print('out.neo.pflux_th not found; nothing to plot.')
    sys.exit(0)

# Grid: n_species, n_energy, n_xi, n_theta, theta(1:n_theta), n_radial, r(:)
tok = open(simdir + '/out.neo.grid').read().split()
n_species = int(tok[0])
n_theta = int(tok[3])
theta = np.array(tok[4:4 + n_theta], dtype=float)
n_radial = int(tok[4 + n_theta])

ir = 0  # radial index to plot (0-based), relevant only if n_radial > 1

# Species-major within each radius: all theta for species 1, then species 2, ...
data = np.loadtxt(simdir + '/out.neo.pflux_th').reshape(n_radial, n_species, n_theta, 2)
pflux_th = data[ir, :, :, 0]
pflux_gv_th = data[ir, :, :, 1]

# One row per radius. out.neo.transport: 5 leading columns (r, d_phi_sqavg, jpar,
# vtor_0order_th0, uparB_0order), then 8 per species starting with pflux.
# out.neo.transport_gv: r, then 3 per species (pflux_gv, eflux_gv, mflux_gv).
tr = np.atleast_2d(np.loadtxt(simdir + '/out.neo.transport'))[ir]
gv = np.atleast_2d(np.loadtxt(simdir + '/out.neo.transport_gv'))[ir]
pflux = tr[5::8][:n_species]
pflux_gv = gv[1::3][:n_species]

species = range(1, n_species + 1) if spec == 'all' else [int(s) for s in spec.split(',')]

fig, axs = plt.subplots(len(species), 2, figsize=(11, 3.4 * len(species)),
                        squeeze=False, sharex=True)
for row, is_ in enumerate(species):
    i = is_ - 1
    for col, (y, name, fsa_file, fsa_name, fsa) in enumerate([
            (pflux_th[i], 'pflux_th', 'out.neo.transport', 'pflux', pflux[i]),
            (pflux_gv_th[i], 'pflux_gv_th', 'out.neo.transport_gv', 'pflux_gv', pflux_gv[i])]):
        ax = axs[row][col]
        ax.plot(theta, y, '-o', ms=3)
        ax.set_title('species %d' % is_)
        ax.set_ylabel(name)
        ax.set_xlabel(r'$\theta$')
        ax.grid(alpha=0.3)
        ax.text(0.03, 0.97, '%s\n%s = %.6e' % (fsa_file, fsa_name, fsa),
                transform=ax.transAxes, va='top', family='monospace', fontsize=9,
                bbox=dict(boxstyle='round', fc='white', ec='gray'))

fig.tight_layout()
if ftype == 'screen':
    plt.show()
else:
    fig.savefig(simdir + '/out.neo.pflux_th.' + ftype, dpi=150)
