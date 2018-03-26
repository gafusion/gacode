import data
import sys
import numpy as np
from gacodefuncs import *
from cgyro.data import cgyrodata

class cgyrodata_dump(data.cgyrodata):

   def dump_flux(self,fc=0,fig=None):

      self.getflux()

      ns = self.n_species
      t = self.t

      ys = np.sum(self.ky_flux,axis=(2,3))
      for moment in ['n','e','v']:
         if moment == 'n':
            ntag = 'Density~flux'
            mtag = '\Gamma'
            ttag = 'G'
            ftag = 'flux_n'
            y = ys[:,0,:]
         elif moment == 'e':
            ntag = 'Energy~flux'
            mtag = 'Q'
            ttag = 'Q'
            ftag = 'flux_e'
            y = ys[:,1,:]
         elif moment == 'v':
            ntag = 'Momentum~flux'
            mtag = '\Pi'
            ttag = 'Pi'
            ftag = 'flux_v'
            y = ys[:,2,:]

         data  = np.column_stack((self.t,y[0,:]))
         head  = '(cs/a) t     '+ttag+'_1/'+ttag+'_GB'
         fname = 'out.cgyro.'+ftag
         for ispec in range(1,ns,1):
            head = head+'       '+ttag+'_'+str(ispec+1)+'/'+ttag+'_GB'
            data = np.column_stack((data,y[ispec,:]))
         np.savetxt(fname,data,fmt='%.8e',header=head)
         print 'INFO: (dump_flux.py) Created '+fname

   def dump_ky_flux(self,w=0.5,field=0,moment='e',fc=0,fig=None):

      self.getflux()
      ns = self.n_species
      t  = self.t

      ky  = self.ky
      ave = np.zeros((self.n_n,ns))

      field_tag = '\mathrm{Total}'

      if fc == 0:
         ys = np.sum(self.ky_flux,axis=(2))
      else:
         ys = self.ky_flux[:,:,field,:,:]
         if field == 0:
            field_tag = '\phi'
         elif field == 1:
            field_tag = 'A_\parallel'
         else:
            field_tag = 'B_\parallel'

      if moment == 'n':
         ntag = 'Density~flux'
         mtag = '\Gamma'
         ttag = 'G'
         ftag = 'flux_n'
         y = ys[:,0,:,:]
      elif moment == 'e':
         ntag = 'Energy~flux'
         mtag = 'Q'
         ttag = 'Q'
         ftag = 'flux_e'
         y = ys[:,1,:,:]
      elif moment == 'v':
         ntag = 'Momentum~flux'
         mtag = '\Pi'
         ttag = 'Pi'
         ftag = 'flux_v'
         y = ys[:,2,:,:]
      else:
         print 'ERROR (plot_ky_flux.py) Invalid moment.'
         sys.exit()

      # Determine tmin
      imin=iwindow(t,w)

      fname = 'out.cgyro.ky_flux.'+ftag
      
      arr = np.zeros([len(ky),ns+1])
      arr[:,0] = ky
      stag = '# (k_y rho_s'
      for ispec in range(ns):
         for j in range(self.n_n):
            ave = average(y[ispec,j,:],self.t,w)
            arr[j,ispec+1] = ave
         stag = stag+' , s'+str(ispec)
            
      fid = open(fname,'w')
      fid.write('# Moment  : '+mtag+'\n')
      fid.write('# Time    : '+str(self.t[imin])+' < (c_s/a) t < '+str(self.t[-1])+'\n')
      fid.write(stag+')\n')
      np.savetxt(fid,arr,fmt='%.5e')
      fid.close()

      print 'INFO: (plot_ky_flux) Created '+fname
