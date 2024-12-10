from . import data
import sys
import numpy as np
from ..gacodefuncs import *
from .data import cgyrodata

class cgyrodata_dump(data.cgyrodata):

   def dump_flux(self,xin):

      moment = xin['moment']
      
      self.getflux()

      ns = self.n_species
      t = self.t

      ys = np.sum(self.ky_flux,axis=(2,3))
      for moment in ['n','e','v']:
         if moment == 'n':
            ntag = 'Density~flux'
            mtag = r'\Gamma'
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
            mtag = r'\Pi'
            ttag = 'Pi'
            ftag = 'flux_v'
            y = ys[:,2,:]
         else:
            raise ValueError('(dump_flux.py) Invalid moment.')

         data  = np.column_stack((self.t,y[0,:]))
         head  = '(cs/a) t     '+ttag+'_1/'+ttag+'_GB'
         fname = 'out.cgyro.'+ftag
         for ispec in range(1,ns,1):
            head = head+'       '+ttag+'_'+str(ispec+1)+'/'+ttag+'_GB'
            data = np.column_stack((data,y[ispec,:]))
         np.savetxt(fname,data,fmt='%.8e',header=head)
         print('INFO: (dump_flux) Created '+fname)

   def dump_ky_flux(self,xin):

      w      = xin['w']
      moment = xin['moment']
      field  = xin['field']
      fc     = xin['fc']
      
      self.getflux()
      ns = self.n_species
      t  = self.t

      ky  = self.ky
      ave = np.zeros((self.n_n,ns))

      field_tag = r'\mathrm{Total}'

      if fc == 0:
         ys = np.sum(self.ky_flux,axis=2)
      else:
         ys = self.ky_flux[:,:,field,:,:]
         if field == 0:
            field_tag = r'\phi'
         elif field == 1:
            field_tag = r'A_\parallel'
         else:
            field_tag = r'B_\parallel'

      if moment == 'e' or moment == 'phi':
         ntag = 'Energy~flux'
         mtag = 'Q'
         ttag = 'Q'
         ftag = 'flux_e'
         y = ys[:,1,:,:]
      elif moment == 'n':
         ntag = 'Density~flux'
         mtag = r'\Gamma'
         ttag = 'G'
         ftag = 'flux_n'
         y = ys[:,0,:,:]
      elif moment == 'v':
         ntag = 'Momentum~flux'
         mtag = r'\Pi'
         ttag = 'Pi'
         ftag = 'flux_v'
         y = ys[:,2,:,:]
      else:
         raise ValueError('(dump_ky_flux.py) Invalid moment.')

      # Determine tmin
      imin,imax=time_index(t,w)

      fname = 'out.cgyro.ky_flux.'+ftag
      
      arr = np.zeros([ns+1,len(ky)])
      arr[0,:] = ky
      arr[1:,:] = time_average(y,self.t,imin,imax)
      stag = '# (k_y rho_s'
      for ispec in range(ns):
         stag = stag+' , s'+str(ispec)

      arr = np.transpose(arr)
      with open(fname,'w') as fid:
         fid.write('# Moment  : '+mtag+'\n')
         fid.write('# Time    : '+str(self.t[imin])+' < (c_s/a) t < '+str(self.t[imax])+'\n')
         fid.write(stag+')\n')
         np.savetxt(fid,arr,fmt='%.5e')

      print('INFO: (dump_ky_flux) Created '+fname)
