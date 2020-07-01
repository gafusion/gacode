# file processed by 2to3
from __future__ import print_function, absolute_import
from builtins import map, filter, range
from gacodeinput import *
import sys

x = SimpleInput()

x.set_extension('.gen')

# VGEN input parameters
x.add('ER_METHOD','1')
x.add('VEL_METHOD','2')
x.add('ERSPECIES_INDX','2')
x.add('EPAR_FLAG','0')
x.add('NTHETA_MIN','17')
x.add('NTHETA_MAX','39')
x.add('NN_FLAG','0')
x.add('TYPE','DC')

# Perform the parsing
x.read_input('input.vgen')

x.printmsg()        

sys.exit(x.error)


