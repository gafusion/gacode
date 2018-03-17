import wx
import matplotlib
import numpy as np
matplotlib.use('WXAgg')
from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg as FigureCanvas
from gacodeplotdefs import *

data = np.loadtxt('out.locpargen.geo')
data2 = np.loadtxt('out.locpargen.geo.2')

class TabPanel(wx.Panel):
    def __init__(self, parent):
        wx.Panel.__init__(self, parent=parent)
 
        self.figure = plt.Figure()
        self.ax = self.figure.add_subplot(111)
        self.canvas = FigureCanvas(self, -1, self.figure)
        self.sizer = wx.BoxSizer(wx.VERTICAL)
        self.sizer.Add(self.canvas, 1, wx.LEFT | wx.TOP | wx.GROW)
        self.SetSizer(self.sizer)
        self.Fit()

    def draw(self,i,ylabel):

        self.ax.grid(which="majorminor",ls=":")
        self.ax.grid(which="major",ls=":")
        self.ax.set_xlabel(r'$\theta/\pi$',fontsize=GFONTSIZE)
        self.ax.set_ylabel(r'$'+ylabel+'$',color='k',fontsize=GFONTSIZE)
        self.ax.plot(data[:,0]/np.pi,data[:,i])
        self.ax.plot(data2[:,0]/np.pi,data2[:,i],color='m')
        self.ax.set_xlim([-1,1])
        if ylabel == '\\Theta':
            self.ax.plot(data[:,0]/np.pi,data[:,0],'--k')
           
            
class DemoFrame(wx.Frame):
     def __init__(self):
        """Constructor"""        
        wx.Frame.__init__(self, None, wx.ID_ANY, 
                          "Geometry functions",
                          size=(1000,800))
        panel = wx.Panel(self)
 
        notebook = wx.Notebook(panel)

        tab = TabPanel(notebook)
        tab.draw(1,'B')
        notebook.AddPage(tab,"b")
 
        tab = TabPanel(notebook)
        tab.draw(2,'\\partial B/\\partial \\theta')
        notebook.AddPage(tab,"dbdt")
 
        tab = TabPanel(notebook)
        tab.draw(3,'\\partial^2 B / \\partial \\theta^2')
        notebook.AddPage(tab,"dbdt2")

        tab = TabPanel(notebook)
        tab.draw(4,'gsin')
        notebook.AddPage(tab,"gsin")

        tab = TabPanel(notebook)
        tab.draw(5,'B_p')
        notebook.AddPage(tab,"bp")

        tab = TabPanel(notebook)
        tab.draw(6,'R')
        notebook.AddPage(tab,"bigr")

        tab = TabPanel(notebook)
        tab.draw(7,'| \\nabla R |')
        notebook.AddPage(tab,"grad_r")

        tab = TabPanel(notebook)
        tab.draw(8,'J_r')
        notebook.AddPage(tab,"jac_r")

        tab = TabPanel(notebook)
        tab.draw(9,'\\Theta')
        notebook.AddPage(tab,"captheta")
 
        tab = TabPanel(notebook)
        tab.draw(10,'\\partial \ell / \\partial \\Theta')
        notebook.AddPage(tab,"l_t")
 
        sizer = wx.BoxSizer(wx.VERTICAL)
        sizer.Add(notebook, 1, wx.ALL|wx.EXPAND, 5)
        panel.SetSizer(sizer)
        self.Layout()
 
        self.Show()
 
if __name__ == "__main__":
    app = wx.App(False)
    frame = DemoFrame()
    app.MainLoop()
