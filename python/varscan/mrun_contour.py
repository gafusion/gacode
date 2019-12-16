import sys
import os
import numpy as np
import matplotlib.cm as cm
from gacodeplotdefs import *
from pyrats.tgyro.data import TGYROData

simdir  = sys.argv[1]
index   = int(sys.argv[2])-1
imgfile = sys.argv[3]
code    = sys.argv[4]
norm    = int(sys.argv[5])

if len(sys.argv) > 6:
    x = sys.argv[6]
    var1 = x.split(',')[0]
    var2 = x.split(',')[1] 
else:
    var1 = "x"
    var2 = "y"

nz = 128

title = ['\Gamma_e','Q_e','\Gamma_i','Q_i']

if code == 'sum':
    xm = np.loadtxt(simdir+'/out.neo.mrun_x')
    ym = np.loadtxt(simdir+'/out.neo.mrun_y')
    zm = np.loadtxt(simdir+'/out.neo.mrun_z')+np.loadtxt(simdir+'/out.tglf.mrun_z')
else:
    xm = np.loadtxt(simdir+'/out.'+code+'.mrun_x')
    ym = np.loadtxt(simdir+'/out.'+code+'.mrun_y')
    zm = np.loadtxt(simdir+'/out.'+code+'.mrun_z')

if norm > 0:
    x = TGYROData('./')
    if index == 0 or index == 2:
        gb = x.data['Gamma_GB'][0][norm]
    else:
        gb = x.data['Q_GB'][0][norm]
    print('INFO: (mrun_contour) Using flux units from r/a='+str(x.data['r/a'][0][norm]))
    unit = ['10^{19}/m^2/s','MW/m^2','10^{19}/m^2/s','MW/m^2']
else:
    gb = 1.0
    print('INFO: (mrun_contour) Using GyroBohm units')
    unit = ['\Gamma_\mathrm{GB}','Q_\mathrm{GB}','\Gamma_\mathrm{GB}','Q_\mathrm{GB}']

x = xm
y = ym

nx = len(x)
ny = len(y)

z = (zm[:,index]*gb).reshape((ny,nx))

zmin = np.amin(z)
zmax = np.amax(z)
 
dz = (zmax-zmin)/nz
levels = np.arange(zmin,zmax+dz,dz)

fig = plt.figure(figsize=(6,6))
#fig.subplots_adjust(left=0.15, right=0.97, top=1.05, bottom=0.0)

ax = fig.add_subplot(111)
ax.set_title(r'$'+title[index]+'\; ['+unit[index]+']$') 

f = ax.contourf(x,y,z,levels,cmap=cm.jet,origin='lower')
plt.colorbar(f)
ax.set_xlabel(r'$\mathrm{'+var1+'}$')
ax.set_ylabel(r'$\mathrm{'+var2+'}$')
ax.set_xlim([np.amin(x),np.amax(x)])
ax.set_ylim([np.amin(y),np.amax(y)])

if imgfile == 'screen':
    plt.show()
else:
    plt.savefig(imgfile)
    print("INFO: (mrun_contour) Wrote plot to "+imgfile+".")

#-------------------------------------------------------------------

