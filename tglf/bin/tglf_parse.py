from gacodeinput import *
import sys
import tglf_defaults

# fill x with the TGLF defaults
x = tglf_defaults.set_defaults()

x.set_extension('.gen')

# Perform the parsing
x.read_input('input.tglf')

x.printmsg()        

sys.exit(x.error)

