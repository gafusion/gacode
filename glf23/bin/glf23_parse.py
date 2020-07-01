#!/usr/bin/env python
# file processed by 2to3
from __future__ import print_function, absolute_import
from builtins import map, filter, range
from gacodeinput import *
import sys

x = SimpleInput()

x.set_extension('.gen')

# GLF23 input parameters
x.add('USE_TRANSPORT_MODEL','.true.')
x.add('USE_ADIABATIC_ELECTRONS','.false.')
x.add('ALPHA_P_MULT','1.0')
x.add('ALPHA_QUENCH_MULT','1.0')
x.add('VERSION','2')
x.add('LPRINT','0')

x.add('NS','2')
x.add('MASS_1','2.723e-4')
x.add('MASS_2','1.0')
x.add('MASS_3','0.0')
x.add('ZS_1','-1.0')      
x.add('ZS_2','1.0') 
x.add('ZS_3','0.0')

x.add('KY','0.3')

x.add('RLNS_1','1.0')
x.add('RLNS_2','1.0') 
x.add('RLNS_3','0.0')
x.add('RLTS_1','3.0')
x.add('RLTS_2','3.0')
x.add('RLTS_3','0.0')
x.add('VPAR_SHEAR_1','0.0')
x.add('VPAR_SHEAR_2','0.0')
x.add('VPAR_SHEAR_3','0.0')
x.add('VEXB_SHEAR','0.0')

x.add('TAUS_1','1.0')
x.add('TAUS_2','1.0')
x.add('TAUS_3','0.0')
x.add('AS_1','1.0')
x.add('AS_2','1.0')
x.add('AS_3','0.0')
x.add('BETAE','0.0')
x.add('XNUE','0.0')

x.add('RMIN_SA','0.5')
x.add('RMAJ_SA','3.0')
x.add('Q_SA','2.0')
x.add('SHAT_SA','1.0')
x.add('ALPHA_SA','0.0')
x.add('XWELL_SA','0.0')
x.add('THETA0_SA','0.0')

# Perform the parsing
x.read_input('input.glf23')

x.printmsg()        

sys.exit(x.error)

