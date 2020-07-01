# file processed by 2to3
from __future__ import print_function, absolute_import
from builtins import map, filter, range
#-------------------------------------------------------------------
# tglf_input_gyro.py
#
# PURPOSE
# read input.gyro.gen,  map it to TGLF inputs and write input.tglf
#
# AUTHORS:
# Gary Staebler
# Orso Meneghini
#--------------------------------------------------------------------
#!/usr/bin/env python
from gacodeinput import *
import sys
import tglf_defaults

# read the input.gyro file

g = {}
f = open('input.gyro.gen','r')
for line in f.readlines():
     
   # Remove leading and trailing whitespace from line
    line = line.strip()
#    print(line)
    # Skip blank lines
    if len(line) > 0 :
        tmp = line.split('  ')
        if len(tmp) == 2 :
            val = eval(tmp[0])
            arg = tmp[1]
            g[arg] = val

f.close
#print(g)

# map gyro to tglf
tg = {}

# ky-grid
tg['KYGRID_MODEL'] = 0
tg['NKY'] = g['TOROIDAL_GRID']-1
tg['KY'] = tg['NKY']*g['TOROIDAL_SEP']*g['RHO_STAR']*g['SAFETY_FACTOR']/g['RADIUS']
if tg['KY'] < 0.0 :
    tg['KY'] = -tg['KY']

#control switches
if g['RADIAL_PROFILE_METHOD'] == 1 :
    tg['GEOMETRY_FLAG']= 0
if g['RADIAL_PROFILE_METHOD'] == 5 :
    tg['GEOMETRY_FLAG']= 1
if g['RADIAL_PROFILE_METHOD' ]== 3 :
    print('error gyro RADIAL_PROFILE_METHOD = 3')
if g['N_FIELD'] >= 2 :
    tg['USE_BPER'] = 'T'
if g['N_FIELD'] == 3 :
    tg['USE_BPAR'] = 'T'
if g['ELECTRON_METHOD'] == 1 :
    tg['ADIABATIC_ELECT'] = 'T'
if g['GEO_GRADBCURV_FLAG'] == 0 :
    tg['USE_MHD_RULE'] = 'F'
tg['ALPHA_E'] = g['GAMMA_E_SCALE']
tg['ALPHA_P'] = g['GAMMA_P_SCALE']
tg['ALPHA_MACH'] = g['MACH_SCALE']
tg['DEBYE_FACTOR'] = g['LAMBDA_DEBYE_SCALE']
tg['XNU_FACTOR'] = g['NU_EI_SCALE']
if g['NU_II_SCALE'] > 0.0 :
    print('error GYRO ion collisions turned on')
nion = 1
if g['NI_OVER_NE_2'] > 0.0 :
    nion = nion + 1
if g['NI_OVER_NE_3'] > 0.0:
    nion = nion + 1
if g['NI_OVER_NE_4'] > 0.0:
    nion = nion + 1
if g['NI_OVER_NE_5'] > 0.0:
    nion = nion + 1
tg['NS'] = nion + 1

# plasma data
tg['SIGN_BT'] = g['BTCCW']
tg['SIGN_IT'] = g['IPCCW']
tg['ZS_1'] = -1.0
tg['ZS_2'] = g['Z']
tg['ZS_3'] = g['Z_2']
tg['ZS_4'] = g['Z_3']
tg['ZS_5'] = g['Z_4']
tg['ZS_6'] = g['Z_5']
tg['MASS_1'] = 1.0/(g['MU_ELECTRON']**2)
tg['MASS_2'] = 1.0/(g['MU'])**2
tg['MASS_3'] = (g['MU']/g['MU_2'])**2
tg['MASS_4'] = (g['MU']/g['MU_3'])**2
tg['MASS_5'] = (g['MU']/g['MU_4'])**2
tg['MASS_6'] = (g['MU']/g['MU_5'])**2
tg['RLNS_1'] = g['DLNNDR_ELECTRON']*(1.0-g['EPS_DLNNDR_ELECTRON'])
tg['RLNS_2'] = g['DLNNDR']*(1.0-g['EPS_DLNNDR'])
tg['RLNS_3'] = g['DLNNDR_2']*(1.0-g['EPS_DLNNDR_2'])
tg['RLNS_4'] = g['DLNNDR_3']*(1.0-g['EPS_DLNNDR_3'])
tg['RLNS_5'] = g['DLNNDR_4']*(1.0-g['EPS_DLNNDR_4'])
tg['RLNS_6'] = g['DLNNDR_5']*(1.0-g['EPS_DLNNDR_5'])
tg['RLTS_1'] = g['DLNTDR_ELECTRON']*(1.0-g['EPS_DLNTDR_ELECTRON'])
tg['RLTS_2'] = g['DLNTDR']*(1.0-g['EPS_DLNTDR'])
tg['RLTS_3'] = g['DLNTDR_2']*(1.0-g['EPS_DLNTDR_2'])
tg['RLTS_4'] = g['DLNTDR_3']*(1.0-g['EPS_DLNTDR_3'])
tg['RLTS_5'] = g['DLNTDR_4']*(1.0-g['EPS_DLNTDR_4'])
tg['RLTS_6'] = g['DLNTDR_5']*(1.0-g['EPS_DLNTDR_5'])
tg['TAUS_1'] = 1.0
tg['TAUS_2'] = g['TI_OVER_TE']
tg['TAUS_3'] = g['TI_OVER_TE_2']
tg['TAUS_4'] = g['TI_OVER_TE_3']
tg['TAUS_5'] = g['TI_OVER_TE_4']
tg['TAUS_6'] = g['TI_OVER_TE_5']
tg['AS_1'] = 1.0
tg['AS_2'] = g['NI_OVER_NE']
tg['AS_3'] = g['NI_OVER_NE_2']
tg['AS_4'] = g['NI_OVER_NE_3']
tg['AS_5'] = g['NI_OVER_NE_4']
tg['AS_6'] = g['NI_OVER_NE_5']
tg['VPAR_1'] = -tg['SIGN_IT']*g['MACH']
tg['VPAR_2'] = -tg['SIGN_IT']*g['MACH']
tg['VPAR_3'] = -tg['SIGN_IT']*g['MACH']
tg['VPAR_4'] = -tg['SIGN_IT']*g['MACH']
tg['VPAR_5'] = -tg['SIGN_IT']*g['MACH']
tg['VPAR_6'] = -tg['SIGN_IT']*g['MACH']
tg['VPAR_SHEAR_1'] = -tg['SIGN_IT']*g['GAMMA_P']
tg['VPAR_SHEAR_2'] = -tg['SIGN_IT']*g['GAMMA_P']
tg['VPAR_SHEAR_3'] = -tg['SIGN_IT']*g['GAMMA_P']
tg['VPAR_SHEAR_4'] = -tg['SIGN_IT']*g['GAMMA_P']
tg['VPAR_SHEAR_5'] = -tg['SIGN_IT']*g['GAMMA_P']
tg['VPAR_SHEAR_6'] = -tg['SIGN_IT']*g['GAMMA_P']
tg['VEXB_SHEAR'] = -tg['SIGN_BT']*g['GAMMA_E']
tg['BETAE'] = g['BETAE_UNIT']*g['AMPERE_SCALE']
tg['ZEFF'] = g['Z_EFF']
tg['DEBYE'] = g['LAMBDA_DEBYE']/g['RHO_STAR']
tg['XNUE'] = g['NU_EI']
#S-ALPHA GEOMETRY
tg['RMIN_SA'] = g['RADIUS']
tg['RMAJ_SA'] = g['ASPECT_RATIO']
tg['Q_SA'] = g['SAFETY_FACTOR']
if tg['Q_SA'] < 0.0 :
   tg['Q_SA'] = -tg['Q_SA']
tg['SHAT_SA'] = g['SHEAR']
tg['ALPHA_SA'] = g['SHIFT']
tg['XWELL_SA'] = 0.0
tg['THETA0_SA'] = 0.0
tg['B_MODEL_SA'] = 1
tg['FT_MODEL_SA'] = 1                   
# MILLER GEOMETRY
tg['RMIN_LOC'] = g['RADIUS']
tg['RMAJ_LOC'] = g['ASPECT_RATIO']
tg['ZMAJ_LOC'] = g['ZMAG']
tg['DRMINDX_LOC'] = 1.0
tg['DRMAJDX_LOC'] = g['SHIFT']
tg['DZMAJDX_LOC'] = g['DZMAG']
tg['KAPPA_LOC'] = g['KAPPA']
tg['S_KAPPA_LOC'] = g['S_KAPPA']
tg['DELTA_LOC'] = g['DELTA']
tg['S_DELTA_LOC'] = g['S_DELTA']
tg['ZETA_LOC'] = g['ZETA']
tg['S_ZETA_LOC'] = g['S_ZETA']
tg['Q_LOC'] = g['SAFETY_FACTOR']
if tg['Q_LOC'] < 0.0 :
   tg['Q_LOC'] = -tg['Q_LOC']
tg['Q_PRIME_LOC'] = g['SHEAR']*(tg['Q_LOC']/tg['RMIN_LOC'])**2
tg['KX0_LOC'] = 0.0
# compute P_PRIME_LOC
scale = g['GEO_BETAPRIME_SCALE']/(8.0*3.1415927)
betae = g['BETAE_UNIT']
rloc = g['RADIUS']
qloc = tg['Q_LOC']
dpdr = g['DLNNDR_ELECTRON'] + g['DLNTDR_ELECTRON'] +                              \
            g['NI_OVER_NE']*g['TI_OVER_TE']*(g['DLNTDR'] + g['DLNNDR']) +          \
            g['NI_OVER_NE_2']*g['TI_OVER_TE_2']*(g['DLNTDR_2'] + g['DLNNDR_2']) +  \
            g['NI_OVER_NE_3']*g['TI_OVER_TE_3']*(g['DLNTDR_3'] + g['DLNNDR_3']) +  \
            g['NI_OVER_NE_4']*g['TI_OVER_TE_4']*(g['DLNTDR_4'] + g['DLNNDR_4']) +  \
            g['TI_OVER_TE_5']*g['NI_OVER_NE_5']*(g['DLNTDR_5'] + g['DLNNDR_5'])  
tg['P_PRIME_LOC'] = -scale*betae*(qloc/rloc)*dpdr

#print(tg)                   
#transfer overwrites to input.tglf dictionary tin
t = SimpleInput()
t = tglf_defaults.set_defaults()

t.set_extension('')

tin = {}

for item in tg:
    if tg[item] != t.data_dict[item] :
        tin[item] = tg[item]

#print(tin)
# write input.tglf
f = open('input.tglf','w')
for item in tin:
    f.write(item+'='+str(tin[item]))
    f.write('\n')

f.close
print('input.gyro succesfully converted to input.tglf')
