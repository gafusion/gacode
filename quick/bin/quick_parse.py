from gacodeinput import *
import sys

x = SimpleInput()

x.set_extension('.gen')

# QUICK input parameters
x.add('N_ION','2')
x.add('Z1','1')
x.add('Z2','6')
x.add('Z3','1')

# Perform the parsing
x.read_input('input.quick')

x.printmsg()        

sys.exit(x.error)


