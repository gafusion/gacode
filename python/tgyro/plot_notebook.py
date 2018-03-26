#-------------------------------------------------------------
# notebook.py
#
# PURPOSE:
#  Notebook plotter to see tgyro results.
#-------------------------------------------------------------

import wx
import matplotlib
import string
import sys
import re
import numpy as np
matplotlib.use('WXAgg')
from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg as FigureCanvas
import matplotlib.pyplot as plt
from matplotlib import rc
from gacodefuncs import *
from tgyro.data import tgyrodata

rc('text',usetex=True)
rc('font',size=18)

simdir = sys.argv[1]
units  = int(sys.argv[2])
ext    = sys.argv[3]

sim = tgyrodata(simdir)

print 'Number of ions  : ',sim.n_ion
print 'Number of radii : ',sim.n_r
print 'Evolution eqns  : ',sim.n_evolve
print 'Completed iter  : ',sim.n_iterations

n   = sim.n_iterations
x   = sim.data['r/a'][0]
ggb = sim.data['Gamma_GB'][n]
qgb = sim.data['Q_GB'][n]
pgb = sim.data['Pi_GB'][n]

n_ion = sim.n_ion

def plot_gen(ax,tag):

    ax.grid(which="major",ls="-",alpha=0.4)
    ax.grid(which="minor",ls=":",alpha=0.4)
    ax.set_xlabel('$r/a$')
 
    tot=r'$\mathbf{total}$'
    tar=r'$\mathbf{target}$'

    if tag == 'eflux_e_target':
        if units == 0:
            ax.plot(x,sim.data['eflux_e_tot'][n],label=tot)
            ax.plot(x,sim.data['eflux_e_target'][n],label=tar)
            ax.set_ylabel('$Q_e/Q_{GB}$',color='k')
        else:
            ax.plot(x,sim.data['eflux_e_tot'][n]*qgb,label=tot)
            ax.plot(x,sim.data['eflux_e_target'][n]*qgb,label=tar)
            ax.set_ylabel('$Q_e [MW/m^2]$',color='k')
    elif tag == 'eflux_i_target':
        if units == 0:
            ax.plot(x,sim.data['eflux_i_tot'][n],label=tot)
            ax.plot(x,sim.data['eflux_i_target'][n],label=tar)
            ax.set_ylabel('$Q_i/Q_{GB}$',color='k')
        else:
            ax.plot(x,sim.data['eflux_i_tot'][n]*qgb,label=tot)
            ax.plot(x,sim.data['eflux_i_target'][n]*qgb,label=tar)
            ax.set_ylabel('$Q_i~[MW/m^2]$',color='k')
    elif tag == 'pflux_e_target':
        if units == 0:
            ax.plot(x,sim.data['pflux_e_tot'][n],label=tot)
            ax.plot(x,sim.data['pflux_e_target'][n],label=tar)
            ax.set_ylabel('$\Gamma_e/\Gamma_{GB}$',color='k')
        else:
            ax.plot(x,sim.data['pflux_e_tot'][n]*ggb,label=tot)
            ax.plot(x,sim.data['pflux_e_target'][n]*ggb,label=tar)
            ax.set_ylabel('$\Gamma_e~[10^{19}/m^2/s]$',color='k')
            
    for i in range(n_ion):
        pstr = 'pflux_i'+str(i+1)
        if tag == pstr+'_target':
            if units == 0:
                ax.plot(x,sim.data[pstr+'_tot'][n],label=tot)
                ax.plot(x,sim.data[pstr+'_target'][n],label=tar)
                ax.set_ylabel('$\Gamma_{i'+str(i+1)+'}/\Gamma_{GB}$',color='k')
            else:
                ax.plot(x,sim.data[pstr+'_tot'][n]*ggb,label=tot)
                ax.plot(x,sim.data[pstr+'_target'][n]*ggb,label=tar)
                ax.set_ylabel('$\Gamma_{i'+str(i+1)+'}~[10^{19}/m^2/s]$',color='k')
            break

    # Smooth curves

    init=r'$\mathbf{initial}$'
    fin=r'$\mathbf{final}$'

    if tag == 'te':
        xf,pf = smooth_pro(sim.data['r/a'][0],sim.data['a/Lte'][0],sim.data['te'][0],64)
        ax.plot(xf,pf,color='black',label=init)
        xf,pf = smooth_pro(sim.data['r/a'][0],sim.data['a/Lte'][n],sim.data['te'][n],64)
        ax.plot(xf,pf,color='magenta',label=fin)
        ax.set_ylabel(r'$\mathrm{T_e~[keV]}$')
        # Dots
        ax.plot(sim.data['r/a'][0],sim.data['te'][0],'o',color='k')
        ax.plot(sim.data['r/a'][0],sim.data['te'][n],'o',color='k')
    elif tag == 'ti':
        xf,pf = smooth_pro(sim.data['r/a'][0],sim.data['a/Lti1'][0],sim.data['ti1'][0],64)
        ax.plot(xf,pf,color='black',label=init)
        xf,pf = smooth_pro(sim.data['r/a'][0],sim.data['a/Lti1'][n],sim.data['ti1'][n],64)
        ax.plot(xf,pf,color='magenta',label=fin)
        ax.set_ylabel(r'$\mathrm{T_i~[keV]}$')
        # Dots
        ax.plot(sim.data['r/a'][0],sim.data['ti1'][0],'o',color='k')
        ax.plot(sim.data['r/a'][0],sim.data['ti1'][n],'o',color='k')
    elif tag == 'ne':
        xf,pf = smooth_pro(sim.data['r/a'][0],sim.data['a/Lne'][0],sim.data['ne'][0],64)
        ax.plot(xf,pf,color='black',label=init)
        xf,pf = smooth_pro(sim.data['r/a'][0],sim.data['a/Lne'][n],sim.data['ne'][n],64)
        ax.plot(xf,pf,color='magenta',label=fin)
        ax.set_ylabel(r'$\mathrm{n_e~[10^19/m^3]}$')
        # Dots
        ax.plot(sim.data['r/a'][0],sim.data['ne'][0],'o',color='k')
        ax.plot(sim.data['r/a'][0],sim.data['ne'][n],'o',color='k')

    for i in range(n_ion):
        if tag == 'ni'+str(i+1):
            xf,pf = smooth_pro(sim.data['r/a'][0],sim.data['a/L'+tag][0],sim.data[tag][0],64)
            ax.plot(xf,pf,color='black',label=init)
            xf,pf = smooth_pro(sim.data['r/a'][0],sim.data['a/L'+tag][n],sim.data[tag][n],64)
            ax.plot(xf,pf,color='magenta',label=fin)
            ax.set_ylabel(r'$\mathrm{n_i~[10^19/m^3]}$')
            # Dots
            ax.plot(sim.data['r/a'][0],sim.data[tag][0],'o',color='k')
            ax.plot(sim.data['r/a'][0],sim.data[tag][n],'o',color='k')
            break
    
    if 'flux' in tag:
        ax.legend(loc=2)
    else:
        ax.legend(loc=1)

#-------------------------------------------------------------------------------------

class TabPanel(wx.Panel):
    def __init__(self, parent):
        wx.Panel.__init__(self, parent=parent)
 
        self.figure = plt.Figure()
        self.figure.subplots_adjust(left=0.07,right=0.95)
        self.ax = self.figure.add_subplot(111)
        self.canvas = FigureCanvas(self, -1, self.figure)
        self.sizer = wx.BoxSizer(wx.VERTICAL)
        self.sizer.Add(self.canvas, 1, wx.LEFT | wx.TOP | wx.GROW)
        self.SetSizer(self.sizer)
        self.Fit()

    def draw(self,tag):
        plot_gen(self.ax,tag)
        
#-------------------------------------------------------------------------------------

class DemoFrame(wx.Frame):
  def __init__(self):
        """Constructor"""        
        wx.Frame.__init__(self, None, wx.ID_ANY, 
                          "TGYRO plotting notebook",
                          size=(1100,600))
        panel = wx.Panel(self)
 
        notebook = wx.Notebook(panel)

        tab = TabPanel(notebook)
        tab.draw('eflux_e_target')
        notebook.AddPage(tab,'*eflux_e')

        tab = TabPanel(notebook)
        tab.draw('eflux_i_target')
        notebook.AddPage(tab,'*eflux_i')

        tab = TabPanel(notebook)
        tab.draw('pflux_e_target')
        notebook.AddPage(tab,'*pflux_e')

        for i in range(n_ion):
            tab = TabPanel(notebook)
            tab.draw('pflux_i'+str(i+1)+'_target')
            notebook.AddPage(tab,'*pflux_i'+str(i+1))
           
        tab = TabPanel(notebook)
        tab.draw('te')
        notebook.AddPage(tab,'Te')

        tab = TabPanel(notebook)
        tab.draw('ti')
        notebook.AddPage(tab,'Ti')

        tab = TabPanel(notebook)
        tab.draw('ne')
        notebook.AddPage(tab,'ne')

        for i in range(n_ion):
            tab = TabPanel(notebook)
            tab.draw('ni'+str(i+1))
            notebook.AddPage(tab,'ni'+str(i+1))

        sizer = wx.BoxSizer(wx.VERTICAL)
        sizer.Add(notebook, 1, wx.ALL|wx.EXPAND, 5)
        panel.SetSizer(sizer)
        self.Layout()
 
        self.Show()
 
#-------------------------------------------------------------------------------------

if __name__ == "__main__":

  if ext == 'screen':

    # On-screen wxpython notebook
      
    app = wx.App(False)
    frame = DemoFrame()
    app.MainLoop()
    
  else:

    # Generate plots

    list=['eflux_e_target','eflux_i_target','pflux_e_target']
    for i in range(n_ion):
        list.append('pflux_i'+str(i+1)+'_target')
    list=list+['te','ti','ne']
    for i in range(n_ion):
        list.append('ni'+str(i+1)+'_target')

    for x in list:
        figure = plt.figure(figsize=(9,6))
        figure.subplots_adjust(left=0.12,right=0.95,bottom=0.16)
        ax = figure.add_subplot(111)
        plot_gen(ax,x)
        pfile = 'out.'+x+'.'+ext
        plt.savefig(pfile)
        print 'INFO: (notebook.py) Wrote '+pfile
