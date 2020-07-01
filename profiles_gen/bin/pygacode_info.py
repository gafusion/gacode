import sys
from pygacode import expro
from pygacode import gapystr_get as gapystr

filename = sys.argv[1]

expro.expro_read(filename)

ni = expro.expro_n_ion

print('Header information')
print(gapystr(expro.expro_head_original))
print(gapystr(expro.expro_head_statefile))
print(gapystr(expro.expro_head_gfile))
print(gapystr(expro.expro_head_cerfile))
print(' shot : {}'.format(expro.expro_shot))
print(' nexp : {}'.format(expro.expro_n_exp))
print(' ')
print('Ion types')
print(' ions : {}'.format(gapystr(expro.expro_name)[:ni]))
print(' type : {}'.format(gapystr(expro.expro_type)[:ni]))
print(' mass : {}'.format(expro.expro_mass[:ni]))
print(' ')
print('Shot figures of merit')
print(' Ip   [MA] : {:.4f}'.format(expro.expro_current))
print(' B     [T] : {:.4f}'.format(expro.expro_bcentr))
print(' betap [-] : {:.4f}'.format(expro.expro_betap))
print(' betan [%] : {:.4f}'.format(expro.expro_betan*100))
print(' tau   [s] : {:.4f}'.format(expro.expro_tau))
print(' tau98 [s] : {:.4f}'.format(expro.expro_tau98y2))


