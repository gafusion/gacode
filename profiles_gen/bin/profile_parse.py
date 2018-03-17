from gacodeinput import *
import sys

inputfile = sys.argv[1]

x = ProfileInput()

x.set_extension('.gen')

# input.profiles scalar parameters

x.add('SHOT','0')
x.add('N_ION','1')
x.add('N_EXP','41')
x.add('BT_EXP','-2.0')
x.add('IP_EXP','1.0')
x.add('RVBV','1.0')
x.add('ARHO_EXP','0.3')
 
# Perform the parsing
x.read_profile(inputfile)

x.printmsg()        

sys.exit(x.error)

