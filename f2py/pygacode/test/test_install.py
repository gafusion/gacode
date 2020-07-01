import os
import sys
import numpy
import tempfile
from pygacode import *
from pygacode import expro

error_message = 'pygacode installation not successful'

filename = os.path.split(__file__)[0] + '/input.gacode'

for k in range(3):
    # read the file
    print('READ: ' + filename)
    expro.expro_read(filename)

    # test that expro_name is a numpy array
    assert isinstance(expro.expro_name, numpy.ndarray), error_message
    # it is one dimensional
    assert len(expro.expro_name.shape) == 1, error_message
    # it has 20 elements
    assert expro.expro_name.shape[0] == 20, error_message
    # the strings that are full are 10 elements long
    assert len(expro.expro_name[0]) == 10, error_message
    # they are bytes (true both in Python 2 and 3)
    assert isinstance(expro.expro_name[0], bytes), error_message

    # test that expro_list returns some fields
    assert len(expro.expro_list) > 0, error_message

    # test gapystr_get returns list
    assert isinstance(gapystr_get(expro.expro_list), list), error_message
    # test gapystr_get returns list of strings
    assert all([isinstance(x, str) for x in gapystr_get(expro.expro_list)]), error_message

    # print some content
    print(gapystr_get(expro.expro_list))
    print(gapystr_get(expro.expro_name))
    print(gapystr_get(expro.expro_type))
    print('nexp  : ' + str(expro.expro_n_exp))
    print('bunit : ' + str(expro.expro_bunit))

    # update the ions labels
    new_name = ['D', 'C', 'D']
    new_type = ['[therm]', '[therm]', '[fast]']
    expro.expro_name[:] = gapystr_set(new_name)
    expro.expro_type[:] = gapystr_set(new_type)
    assert gapystr_get(expro.expro_name)[:3] == new_name, error_message
    assert gapystr_get(expro.expro_type)[:3] == new_type, error_message

    # write a new file
    tmpdir = tempfile.gettempdir() + os.sep
    filename = tmpdir + 'input.gacode'
    print('WRITE: ' + filename)
    expro.expro_write(filename)

    print('')

print('SUCCESS: pygacode installation tests have passed')
