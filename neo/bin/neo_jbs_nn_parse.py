from gacodeinput import *
import sys

x = SimpleInput()

x.set_extension('.gen')

# JBSNNGEN input parameters
x.add('IMPURITIES_SPEC','1')
x.add('IMPURITIES_SPEC','2')
x.add('PLASMA_SPEC','1')
x.add('PLASMA_SPEC','2')
x.add('NTHETA_MIN','17')
x.add('NTHETA_MAX','39')
x.add('TYPE','DC')

# Perform the parsing
x.read_input('input.jbsnn')

x.printmsg()        

sys.exit(x.error)


