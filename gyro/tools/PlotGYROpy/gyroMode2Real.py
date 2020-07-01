#!/usr/bin/env python
# file processed by 2to3
from __future__ import print_function, absolute_import
from builtins import map, filter, range

import numpy
import optparse
import tables
import sys

def dataAtTorAngle(modeData,n0,n_n,nu,omega,phi):
     """ Given a dictionary containing the mode data, and the 
      transform variables (n0,n_n,omega,nu) show the data at a
      particular toroidal angle

     """
     firstkey=list(modeData.keys())[0]
     # Shape will 
     modeShape=modeData[firstkey].shape
     nuShape=modeData[firstkey].shape
     
     for key in list(modeData.keys()):
       if (n0==0):
           istart=2
           for iphi in range(1,nphi+1):
             real_buff[:,:,:,iphi]=real(buffn[:,:,:,1])
       else:
         istart=1
         real_buff[:,:,:,:]=0.
         #Get alpha coordinate on either the coarse or fine mesh.
         # Include doppler shift here
         if (iscoarse):
           alpha_loc=alpha_phi[:,:,iphi]+omega_exp*t_current
         else:
           alpha_loc=alpha_phi_fine[:,:,iphi]+omega_exp*t_current
         for im in range(istart,n_n+1):
            nn=n0+(im-1)*d_n
            for ikin in range(1,n3+1):
              real_buff[:,:,ikin,iphi]=real_buff[:,:,ikin,iphi]\
                    +2.*real(buffn[:,:,ikin,im]*exp(-c_i*nn*alpha_loc[:,:]))

def main():
    """ SEK: Once the pretty print can be figured out, then we can look
    at completing the standalone executable
    """
    parser = optparse.OptionParser(usage="%prog [options] inputFile")
    parser.add_option('-i', '--input', dest='input',
                      help='Name of input file if not specified as argument.',
                      default='')

    options, args = parser.parse_args()

    # Process arguments
    if len(args) > 1:
      parser.print_usage()
      return
    elif len(args)==1:
      inputFile=args[0]
    else:
      if options.input == '':
        print("Must specify an input file")
        return
      else:
        inputFile=options.input
    
    return

if __name__ == "__main__":
        main()
