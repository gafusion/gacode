# file processed by 2to3
from __future__ import print_function, absolute_import
from builtins import map, filter, range
from gacodeinput import *
import sys

x = SimpleInput()

x.set_extension('.gen')

# PROFILES_3d input parameters
x.add('PSI_NORM','0.45')
x.add('NTHETA_FOURIER','6')
x.add('A_METERS','0.6')
x.add('IMPURITY_FLAG','0')
x.add('Z_IMP','6')
x.add('RZTEST_MODEL','0')
x.add('NP_MASK0','1')
x.add('NP_MASK1','1')
x.add('NP_MASK2','1')
x.add('NP_MASK3','1')
x.add('NP_MASK4','0')
x.add('NP_MASK5','0')
x.add('NP_MASK6','0')
x.add('NP_MASK7','0')
x.add('NP_MASK8','0')
x.add('NP_MASK9','0')
x.add('NP_MASK10','0')
x.add('NP_MASK11','0')
x.add('NP_MASK12','0')
x.add('NP_MASK13','0')
x.add('NP_MASK14','0')
x.add('NP_MASK15','0')
x.add('NP_MASK16','0')

# Perform the parsing
x.read_input('input.profiles_3d')

x.printmsg()

sys.exit(x.error)
