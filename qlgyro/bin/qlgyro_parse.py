
from gacodeinput import *
import sys

x = SimpleInput()

x.set_extension('.gen')

# QLGYRO input parameters
x.add('N_PARALLEL', '1')
x.add('N_RUNS', '1')
x.add('GAMMA_E', '0.0')
x.add('CODE', '0')
x.add('KYGRID_MODEL', '1')
x.add('KY', '0.3')
x.add('NKY', '12')
x.add('XNU_MODEL', '3')
x.add('AUTO_BOX_SIZE', '0')
x.add('KX_MAX_BOX', '1000000.0')
x.add('SAT_RULE', '1')
x.add('N_PX0', '1')
x.add('PX0GRID_MODEL', '1')
x.add('RESTART_MODE', '0')

# Perform the parsing
x.read_input('input.qlgyro')

x.printmsg()

if x.data_dict['CODE'] == '0':
    print('GYRO')
else:
    print('CGYRO')

sys.exit(x.error)


