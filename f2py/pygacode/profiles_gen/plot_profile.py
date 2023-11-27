#-------------------------------------------------------------
# plot_profile.py
#
# PURPOSE:
#  Notebook plotter to see input.gacode profiles
#-------------------------------------------------------------

import os
import wx
import sys
import matplotlib
import numpy as np
matplotlib.use('WXAgg')
from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg as FigureCanvas
import matplotlib.pyplot as plt
from ..gacodefuncs import *
from pygacode import expro
from pygacode import gapystr_get as gapystr

matplotlib.rc('text',usetex=True)
matplotlib.rc('font',size=18)

simdir = './'
wdir = os.path.realpath(simdir)
infile = sys.argv[1]
rvar = sys.argv[2]
rmin = sys.argv[3]
rmax = sys.argv[4]
ext = sys.argv[5]
loc = int(sys.argv[6])
dot = int(sys.argv[7])
therm = int(sys.argv[8])

m1 = 0 ; m2 = 0

alls = not bool(therm)

expro.expro_read(infile,0)

def plotit(ax,x,y,ystr):

   global m1,m2,dot,rvar
   
   ax.plot(x[m1:m2],y[m1:m2],label=r'$'+ystr+'$')
   if dot:
      ax.plot(x[m1:m2],y[m1:m2],'o',color='k',alpha=0.3,ms=4)     
   return

def plot_select(ax,tag):

   global m1,m2,dot,rvar

   if rvar == 'r':
      x = expro.expro_rmin
      xlabel = r'$r/a$'
   elif rvar == 'rho':
      x = expro.expro_rho
      xlabel = r'$\rho$'
   elif rvar == 'pflux':
      x = abs(expro.expro_polflux)
      xlabel = r'$\psi_\mathrm{pol}$'

   x = x/max(x)
   n = expro.expro_n_ion
   
   # normalization
   a = expro.expro_rmin[-1]
   csa = expro.expro_cs/a
   if csa[0] == 0.0:
      # csa may be zero if profile=null
      csa[:] = 1.0

   sname = gapystr(expro.expro_name)
   stype = gapystr(expro.expro_type)

   # Set roman font
   for p in range(n):
      sname[p] = '\mathrm{'+sname[p]+'}'
   
   # Set x-range
   m1 = 0 ; m2 = len(x)
   if rmax != 'auto':
      m2 = np.argmin(np.abs(x-float(rmax)))
   if rmin != 'auto':
      m1 = np.argmin(np.abs(x-float(rmin)))
      
   ax.set_xlim([x[m1],x[m2-1]])
   ax.set_xlabel(xlabel)

   m=m2
   if tag == 'Bunit':
      y = expro.expro_bunit ; ystr = 'B_\mathrm{unit}' ; plotit(ax,x,y,ystr)

   if tag == 'gam_e':
      y = expro.expro_gamma_e ; ystr = '(a/c_s) \gamma_E' ; plotit(ax,x,y/csa,ystr)

   if tag == 'gam_p':
      y = expro.expro_gamma_p ; ystr = '(a/c_s) \gamma_p' ; plotit(ax,x,y/csa,ystr) 

   if tag == 'Mach':
      y = expro.expro_mach ; ystr = r'M' ; plotit(ax,x,y,ystr)

   if tag == 'r':
      y = expro.expro_rho ; ystr = '\\rho' ; plotit(ax,x,y,ystr)
      y = expro.expro_polflux ; ystr = '\psi' ; plotit(ax,x,y/y[-1],ystr)

   if tag == 'R':
      y = expro.expro_rmaj ; ystr = r'R_0' ; plotit(ax,x,y,ystr)
      y = expro.expro_drmaj ; ystr = r'dR_0/dr' ; plotit(ax,x,y,ystr)

   if tag == 'Z':
      y = expro.expro_zmag ; ystr = r'Z_0' ; plotit(ax,x,y,ystr)
      y = expro.expro_dzmag ; ystr = r'dZ_0/dr' ; plotit(ax,x,y,ystr)

   if tag == 'kappa':
      y = expro.expro_kappa ; ystr = r'\kappa' ; plotit(ax,x,y,ystr)
      y = expro.expro_skappa ; ystr = r's_\kappa' ; plotit(ax,x,y,ystr)

   if tag == 'n':
      # ne
      y = expro.expro_ne ; ystr = 'n_e' ; plotit(ax,x,y,ystr)
      # ni
      for p in range(n):
         if alls or stype[p] == '[therm]':
            y = expro.expro_ni[p,:] ; ystr = 'n_i~['+sname[p]+']' ; plotit(ax,x,y,ystr)

   if tag == 'Ln':
      # Lne
      y = a*expro.expro_dlnnedr ; ystr = 'a/L_{ne}' ; plotit(ax,x,y,ystr)
      # Lni
      for p in range(n):
         if alls or stype[p] == '[therm]':
            y = a*expro.expro_dlnnidr[p,:] ; ystr = 'a/L_{ni}~['+sname[p]+']' ; plotit(ax,x,y,ystr)

   if tag == 'sn':
      # sne
      y = a*expro.expro_sdlnnedr ; ystr = 's_{ne}' ; plotit(ax,x,y,ystr)
      # sni
      for p in range(n):
         if alls or stype[p] == '[therm]':
            y = a*expro.expro_sdlnnidr[p,:] ; ystr = 's_{ni}~['+sname[p]+']' ; plotit(ax,x,y,ystr)

   if tag == 'T':
      y = expro.expro_te ; ystr = 'T_e' ; plotit(ax,x,y,ystr)
      for p in range(n):
         if alls or stype[p] == '[therm]':
            y = expro.expro_ti[p,:] ; ystr = 'T_i~['+sname[p]+']' ; plotit(ax,x,y,ystr)

   if tag == 'LT':
      y = a*expro.expro_dlntedr ; ystr = 'a/L_{Te}' ; plotit(ax,x,y,ystr)
      for p in range(n):
         if alls or stype[p] == '[therm]':
            y = a*expro.expro_dlntidr[p,:] ; ystr = 'a/L_{Ti}~['+sname[p]+']' ; plotit(ax,x,y,ystr)

   if tag == 'sT':
      y = a*expro.expro_sdlntedr ; ystr = 's_{Te}' ; plotit(ax,x,y,ystr)
      for p in range(n):
         if alls or stype[p] == '[therm]':
            y = a*expro.expro_sdlntidr[p,:] ; ystr = 's_{Ti}~['+sname[p]+']' ; plotit(ax,x,y,ystr)

   if tag == 'j':
       y = expro.expro_jbs ; ystr = 'J_\mathrm{bs}' ; plotit(ax,x,y,ystr)
       y = expro.expro_johm ; ystr = 'J_\mathrm{ohm}' ; plotit(ax,x,y,ystr)
       y = expro.expro_jbstor ; ystr = 'J_\mathrm{bs,tor}' ; plotit(ax,x,y,ystr)

   if tag == 'sigma':
       y = expro.expro_sigmapar ; ystr = '\sigma_\parallel' ; plotit(ax,x,y,ystr)

   if tag == 'sin':
       y = np.arcsin(expro.expro_delta) ; ystr = 's_1 = \sin^{-1}\delta' ; plotit(ax,x,y,ystr)
       y = -expro.expro_zeta         ; ystr = 's_2 = -\zeta' ; plotit(ax,x,y,ystr)
       y = expro.expro_shape_sin3    ; ystr = 's_3' ; plotit(ax,x,y,ystr)
       y = expro.expro_shape_sin4    ; ystr = 's_4' ; plotit(ax,x,y,ystr)
       y = expro.expro_shape_sin5    ; ystr = 's_5' ; plotit(ax,x,y,ystr)
       y = expro.expro_shape_sin6    ; ystr = 's_6' ; plotit(ax,x,y,ystr)
         
   if tag == 'cos':
       y = expro.expro_shape_cos0 ; ystr = 'c_0' ; plotit(ax,x,y,ystr)
       y = expro.expro_shape_cos1 ; ystr = 'c_1' ; plotit(ax,x,y,ystr)
       y = expro.expro_shape_cos2 ; ystr = 'c_2' ; plotit(ax,x,y,ystr)
       y = expro.expro_shape_cos3 ; ystr = 'c_3' ; plotit(ax,x,y,ystr)
       y = expro.expro_shape_cos4 ; ystr = 'c_4' ; plotit(ax,x,y,ystr)
       y = expro.expro_shape_cos5 ; ystr = 'c_5' ; plotit(ax,x,y,ystr)
       y = expro.expro_shape_cos6 ; ystr = 'c_6' ; plotit(ax,x,y,ystr)

   if tag == 'ssin':
       y = expro.expro_sdelta       ; ystr = 's_1 = \sin^{-1}\delta' ; plotit(ax,x,y,ystr)
       y = -expro.expro_szeta       ; ystr = 's_2 = -\zeta' ; plotit(ax,x,y,ystr)
       y = expro.expro_shape_ssin3  ; ystr = 's_3' ; plotit(ax,x,y,ystr)
       y = expro.expro_shape_ssin4  ; ystr = 's_4' ; plotit(ax,x,y,ystr)
       y = expro.expro_shape_ssin5  ; ystr = 's_5' ; plotit(ax,x,y,ystr)
       y = expro.expro_shape_ssin6  ; ystr = 's_6' ; plotit(ax,x,y,ystr)

   if tag == 'scos':
       y = expro.expro_shape_scos0 ; ystr = 'c_0' ; plotit(ax,x,y,ystr)
       y = expro.expro_shape_scos1 ; ystr = 'c_1' ; plotit(ax,x,y,ystr)
       y = expro.expro_shape_scos2 ; ystr = 'c_2' ; plotit(ax,x,y,ystr)
       y = expro.expro_shape_scos3 ; ystr = 'c_3' ; plotit(ax,x,y,ystr)
       y = expro.expro_shape_scos4 ; ystr = 'c_4' ; plotit(ax,x,y,ystr)
       y = expro.expro_shape_scos5 ; ystr = 'c_5' ; plotit(ax,x,y,ystr)
       y = expro.expro_shape_scos6 ; ystr = 'c_6' ; plotit(ax,x,y,ystr)

   if tag == 'q':
       y = expro.expro_q ; ystr = 'q'
       if y[0] < 0.0:
          y = -y
          ystr = '-'+ystr
       plotit(ax,x,y,ystr)
       y = expro.expro_s ; ystr = 's' ; plotit(ax,x,y,ystr)

   if tag == 'nu':
       y = expro.expro_nuee ; ystr = '(a/c_s) \\nu_{ee}' ; plotit(ax,x,y/csa,ystr)

   if tag == 'rhos':
       y = expro.expro_rhos ; ystr = '\\rho_s/a' ; plotit(ax,x,y/a,ystr)

   if tag == 'V':
       y = expro.expro_vol ; ystr = 'V/a^3' ; plotit(ax,x,y/a**3,ystr)
       y = expro.expro_volp ; ystr = 'V^\prime/a^2' ; plotit(ax,x,y/a**2,ystr)
       y = expro.expro_surf ; ystr = 'S/a^2' ; plotit(ax,x,y/a**2,ystr)

   if tag == 'Pe':
      y = expro.expro_qohme+expro.expro_qbeame+expro.expro_qrfe+expro.expro_qione
      yt = y
      ystr = 'q_{e,\mathrm{aux}}~\mathrm{[MW/m^3]}' ; plotit(ax,x,y,ystr)
      y = expro.expro_qfuse ; ystr = 'q_{\\alpha e}~\mathrm{[MW/m^3]}' ; plotit(ax,x,y,ystr)
      yt = yt+y
      y = -expro.expro_qbrem ; ystr = '-q_\mathrm{brem}~\mathrm{[MW/m^3]}' ; plotit(ax,x,y,ystr)
      yt = yt+y
      y = -expro.expro_qsync ; ystr = '-q_\mathrm{sync}~\mathrm{[MW/m^3]}' ; plotit(ax,x,y,ystr)
      yt = yt+y
      y = -expro.expro_qline ; ystr = '-q_\mathrm{line}~\mathrm{[MW/m^3]}' ; plotit(ax,x,y,ystr)
      yt = yt+y
      y = -expro.expro_qei ; ystr = '-q_{ei}~\mathrm{[MW/m^3]}' ; plotit(ax,x,y,ystr)
      yt = yt+y
      ystr = 'q_{e}^\mathrm{TOT}~\mathrm{[MW/m^3]}' ; plotit(ax,x,yt,ystr)
      
   if tag == 'Pi':
      y = expro.expro_qcxi+expro.expro_qbeami+expro.expro_qrfi+expro.expro_qioni
      yt = y
      ystr = 'q_{i,\mathrm{aux}}~\mathrm{[MW/m^3]}' ; plotit(ax,x,y,ystr)
      y = expro.expro_qfusi ; ystr = 'q_{\\alpha i}~\mathrm{[MW/m^3]}' ; plotit(ax,x,y,ystr)
      yt = yt+y
      y = expro.expro_qei ; ystr = 'q_{ei}~\mathrm{[MW/m^3]}' ; plotit(ax,x,y,ystr)
      yt = yt+y
      ystr = 'q_{i}^\mathrm{TOT}~\mathrm{[MW/m^3]}' ; plotit(ax,x,yt,ystr)

   ax.legend(loc=loc,ncol=2)
       

#-------------------------------------------------------------------------------------

class TabPanel(wx.Panel):
    def __init__(self, parent):
        wx.Panel.__init__(self, parent=parent)
 
        self.figure = plt.Figure()
        self.figure.subplots_adjust(left=0.07,right=0.95)
        self.ax = self.figure.add_subplot(111)
        self.ax.grid(which="both",ls=":")
        self.canvas = FigureCanvas(self, -1, self.figure)
        self.sizer = wx.BoxSizer(wx.VERTICAL)
        self.sizer.Add(self.canvas, 1, wx.LEFT | wx.TOP | wx.GROW)
        self.SetSizer(self.sizer)
        self.Fit()
        
    def draw(self,tag):
        plot_select(self.ax,tag)
  
#-------------------------------------------------------------------------------------

class DemoFrame(wx.Frame):
   
  def __init__(self):

        wx.Frame.__init__(self, None, wx.ID_ANY, 
                          'input.gacode plotting notebook -- '+wdir,
                          size=(1100,800))

        panel = wx.Panel(self)
 
        notebook = wx.Notebook(panel)

        mytabs = ['r','R','Z','kappa','sin','cos','scos','ssin','V','q','Bunit',
                  'n','Ln','sn','T','LT','sT','gam_e','gam_p','Mach','rhos','nu',
                  'j','sigma','Pe','Pi']

        for x in mytabs:
           tab = TabPanel(notebook)
           tab.draw(x)
           notebook.AddPage(tab,x)

        sizer = wx.BoxSizer(wx.VERTICAL)
        sizer.Add(notebook, 1, wx.ALL|wx.EXPAND, 5)
        panel.SetSizer(sizer)
        self.Layout() 
        self.Show()
 

if __name__ == "__main__":

   # On-screen wxpython notebook
      
   app = wx.App(False)
   frame = DemoFrame()
   app.MainLoop()
