from numpy.distutils.core import setup,Extension

ext = Extension('gacode_ext',
                 sources=['expro/expro.f90',
                          'expro/expro_util.f90',
                          'expro/expro_pycomm.f90',
                          'geo/geo.f90'])

setup(py_modules=['pygacode.gacodefuncs',
                  'pygacode.gacodeinput'],
      packages=['pygacode',
                'pygacode.gyro',
                'pygacode.cgyro',
                'pygacode.tgyro',
                'pygacode.test',
                'pygacode.neo',
                'pygacode.profiles_gen'],
      package_data={'pygacode.test': ['input.gacode']},
      ext_modules=[ext]
      )

