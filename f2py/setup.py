import os
import sys

with open(os.path.dirname(os.path.abspath(__file__)) + '/pygacode/version', 'r') as f:
    __version__ = f.read().strip()

# ================================
# utility function for updating conda version
# ================================
if 'conda' in sys.argv:
    import shutil
    import re
    import tempfile
    import subprocess

    tmpdir = tempfile._get_default_tempdir() + os.sep + 'pygacode_conda'
    if os.path.exists(tmpdir):
        shutil.rmtree(tmpdir)
    os.mkdir(tmpdir)
    os.chdir(tmpdir)
    subprocess.run('pip download pygacode==' + __version__, stdout=subprocess.PIPE, shell=True)
    sha = subprocess.run('pip hash pygacode-' + __version__ + '.tar.gz', stdout=subprocess.PIPE, shell=True)
    sha = str(sha.stdout, 'utf8').strip().split(':')[-1]

    subprocess.run('git clone git@github.com:conda-forge/pygacode-feedstock.git', stdout=subprocess.PIPE, shell=True)
    tmp = open('pygacode-feedstock/recipe/meta.yaml').read()
    mod = re.sub('{% set version.*', '{% set version = "' + __version__ + '" %}', tmp)
    mod = re.sub(".*sha256.*", "  sha256: %s" % sha, mod)
    open('pygacode-feedstock/recipe/meta.yaml', "w").write(mod)
    subprocess.run('cd pygacode-feedstock; git diff', shell=True)

    print('WORKING DIRECTORY IS: ' + tmpdir)
    sys.exit()

# ================================
# pip
# ================================
from numpy.distutils.core import setup, Extension

wrapper = Extension('gacode_ext',
                    sources=['expro/expro.f90',
                             'expro/expro_util.f90',
                             'expro/expro_pycomm.f90',
                             'geo/geo.f90',
                             'vis/vis.f90'])

setup(name='pygacode',
      version=__version__,
      description='Python interface to GACODE profile, geometry, and code tools',
      url='https://gacode.io',
      author='General Atomics Theory Group',
      author_email='candy@fusion.gat.com',
      license='MIT',
      py_modules=['pygacode.__init__', 'pygacode.gacodefuncs', 'pygacode.gacodeinput'],
      package_data={'pygacode.test': ['input.gacode'],
                    'pygacode': ['version']},
      packages=['pygacode', 'pygacode.test', 'pygacode.gyro', 'pygacode.cgyro', 'pygacode.tgyro', 'pygacode.neo', 'pygacode.profiles_gen'],
      ext_modules=[wrapper]
      )

# only run tests when installing
if 'install' in sys.argv:
    from pygacode.test import test_install
