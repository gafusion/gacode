from numpy.distutils.core import setup, Extension

wrapper = Extension('gacode',
                    sources=['expro/expro.f90',
                             'expro/expro_util.f90',
                             'expro/expro_pycomm.f90',
                             'geo/geo.f90'])

setup(name='gacode',
      version='0.1',
      description='Python interface to gacode profile and geometry tools.',
      url='https://gafusion.github.io/doc',
      author='General Atomics Theory Group',
      author_email='candy@fusion.gat.com',
      license='MIT',
      packages=['expro',
                'geo',
                'pygacode.cgyro',
                'pygacode.gyro',
                'pygacode.tgyro',
                'pygacode.neo',
                'pygacode.input_gacode'
      ],
      package_dir={'expro':'expro',
                   'geo':'geo',
                   'pygacode.cgyro':'pygacode/cgyro',
                   'pygacode.gyro':'pygacode/gyro',
                   'pygacode.tgyro':'pygacode/tgyro',
                   'pygacode.neo':'pygacode/neo',
                   'pygacode.input_gacode':'pygacode/input_gacode'
      },
      ext_modules=[wrapper]
)
