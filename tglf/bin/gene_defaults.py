#-------------------------------------------------------------------
# gene_defaults.py
#
# PURPOSE
# define defaults for inputs to GENE that are needed for TGLF
#
# AUTHORS:
# Gary Staebler
#--------------------------------------------------------------------
import sys
from gacodeinput import *
gene = SimpleInput()
# note the namelist parser f90nml converts everything to lower case
def set_defaults() :
# box namelist
     gene.add('nky0','16')
     gene.add('kymin','0.05')
     gene.add('n_spec','2')
# general namelist
     gene.add('nonlinear','True')
     gene.add('collision_op','pitch-angle')
     gene.add('coll','0.0')
     gene.add('beta','0.0')
     gene.add('debye2','0.0')
     gene.add('bpar','T')
     gene.add('zeff','1.0')
#  external_contr namelist
     gene.add('exbrate','0.0')
     gene.add('pfsrate','0.0')
# geometry namelist
     gene.add('magn_geometry','miller')
     gene.add('shat','1.0')
     gene.add('q0','2.0')
     gene.add('major_r','3.0')
     gene.add('minor_r','1.0')
     gene.add('major_z','0.0')
     gene.add('trpeps','0.167')
     gene.add('amhd','0.0')
     gene.add('kappa','1.0')
     gene.add('delta','0.0')
     gene.add('s_kappa','0.0')
     gene.add('s_delta','0.0')
     gene.add('drr','0.0')
     gene.add('drz','0.0')
     gene.add('zeta','0.0')
     gene.add('s_zeta','0.0')
     gene.add('dpdx_term','full-drift')
     gene.add('dpdx_pm','0.0')     
     return gene
