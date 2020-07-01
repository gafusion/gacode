import os
import sys

from . import gyro
from . import cgyro
from . import tgyro
from . import neo
from . import profiles_gen

with open(os.path.dirname(__file__) + '/version', 'r') as f:
    __version__ = f.read().strip()

__all__ = ['gyro', 'cgyro', 'tgyro', 'neo', 'profiles_gen', '__version__']

try:
    from gacode_ext import *

    tmp = {}
    exec('from gacode_ext import *', tmp)
    __all__.extend(list(tmp.keys()))


    # The code above will import:
    # - bound_deriv
    # - bound_extrap
    # - bound_interp
    # - expro
    # - expro_acomm
    # - expro_compute_derived
    # - expro_icomm
    # - expro_lcomm
    # - expro_rcomm
    # - expro_scomm
    # - expro_skip_header
    # - expro_tcomm
    # - expro_vcomm
    # - expro_writea
    # - expro_writei
    # - expro_writes
    # - expro_writev
    # - geo
    # - vis
    # - volint

    def gapystr_get(s):
        try:
            n = len(s)
        except:
            return str(s, 'utf-8')

        u = []
        for i in range(len(s)):
            if sys.version_info[0] == 2:
                u.append(str(s[i]).strip())
            else:
                u.append(str(s[i], 'utf-8').strip())
        return u


    def gapystr_set(s, l=10, n=20):
        return [item.ljust(l) for item in s + [''] * n][:n]


    __all__.extend(['gapystr_get', 'gapystr_set'])

except ImportError:
    pass
