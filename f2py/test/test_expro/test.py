from pygacode import expro
import numpy as np
import os, sys

filename = '../data/input.gacode'

for k in range(3):
    print('READ: ' + filename)
    expro.expro_read(filename)

    print(expro.expro_name)
    print(expro.expro_type)

    print('nexp  : ' + str(expro.expro_n_exp))
    print('bunit : ' + str(expro.expro_bunit))

    import tempfile

    tmpdir = tempfile.gettempdir() + os.sep
    filename = tmpdir + 'input.gacode'
    print('WRITE: ' + filename)
    expro.expro_write(filename)

    print('')
