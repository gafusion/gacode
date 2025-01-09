from setuptools import Extension, setup
from setuptools.command.build_ext import build_ext
from importlib.machinery import EXTENSION_SUFFIXES
import os
import shutil
import sys
import subprocess
import sysconfig


class F2PyExtension(Extension):
    def __init__(self, name, fortran_sources, **kwargs):
        super().__init__(name, sources=[], **kwargs)
        self.fortran_sources = fortran_sources


class BuildF2PyExtension(build_ext):
    def build_extension(self, ext):
        output_module = ext.name
        env = os.environ.copy()
        if sys.version_info < (3, 12):
            env["SETUPTOOLS_USE_DISTUTILS"] = "1"
        subprocess.check_call([
            "f2py", "-c", "-m", output_module, *ext.fortran_sources
        ], env=env)

        ext_suffix = sysconfig.get_config_var("EXT_SUFFIX")
        lib_name = f"{output_module}{ext_suffix}"

        ext_path = self.get_ext_fullpath(ext.name)
        ext_dir = os.path.dirname(ext_path)
        if not os.path.exists(ext_dir):
            os.makedirs(ext_dir)

        shutil.move(lib_name, ext_path)


ext = F2PyExtension('gacode_ext',
                    fortran_sources=['expro/expro.f90',
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
      ext_modules=[ext],
      cmdclass={'build_ext': BuildF2PyExtension},
      )
