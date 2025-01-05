project('gacode', 'c', 'fortran',
  version: '1.0',
  license: 'MIT',
  meson_version: '>=0.64.0',
  default_options: ['warning_level=2']
)

add_languages('fortran')

py_mod = import('python')
py = py_mod.find_installation(pure: false)
py_dep = py.dependency()

incdir_numpy = run_command(py,
  ['-c', 'import os; os.chdir(".."); import numpy; print(numpy.get_include())'],
  check : true
).stdout().strip()

incdir_f2py = run_command(py,
    ['-c', 'import os; os.chdir(".."); import numpy.f2py; print(numpy.f2py.get_include())'],
    check : true
).stdout().strip()

inc_np = include_directories(incdir_numpy, incdir_f2py)

gacode_wrapper = custom_target('gacode_wrappers',
  input: [
    'expro/expro.f90',
    'expro/expro_util.f90',
    'expro/expro_pycomm.f90',
    'geo/geo.f90',
  ],
  output: ['gacode_extmodule.c', 'gacode_ext-f2pywrappers.f', 'gacode_ext-f2pywrappers2.f90'],
  command: [py, '-m', 'numpy.f2py', '@INPUT@', '-m', 'gacode_ext', '--lower']
)

py.extension_module('gacode_ext',
  [
    'expro/expro.f90',
    'expro/expro_util.f90',
    'expro/expro_pycomm.f90',
    'geo/geo.f90',
    gacode_wrapper
  ],
  incdir_f2py / 'fortranobject.c',
  include_directories: inc_np,
  dependencies: py_dep,
  install: true
)

subdir('pygacode')
