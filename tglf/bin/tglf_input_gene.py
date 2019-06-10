#-------------------------------------------------------------------
# tglf_input_gene.py
#
# PURPOSE
# read GENE parameters_gene.txt file, map it to TGLF inputs and write input.tglf
# this conversion uses the GENE chioces for Tref, Lref, mref
# even if they are not the usual choices for TGLF. The conversion between B_unit and Bref
# will be done inside TGLF so that the output from TGLF stand alone will be in GENE units Bref
#
# AUTHORS:
# Gary Staebler
# Orso Meneghini
#--------------------------------------------------------------------
#!/usr/bin/env python
from gacodeinput import *
import sys
import tglf_defaults
import gene_defaults
import f90nml
from pprint import pprint

# read the file parameters_gene.txt file

genenamelists = f90nml.read('./parameters_gene.txt')
# pprint(genenamelists)
# unpack the namelist
geneinput={}
# box namelist
geneinput.update(genenamelists['box'])
# general namelist
geneinput.update(genenamelists['general'])
#  external_contr namelist
geneinput.update(genenamelists['external_contr'])
# geometry namelist
geneinput.update(genenamelists['geometry'])
# pprint(geneinput)
# species namelist
species = genenamelists['species']
#pprint(species)
# use the geneinputs and the defaults gend to make a composite genein
gend = SimpleInput()
gened = gene_defaults.set_defaults()
genein = {}
for item in gened.data_dict :
    genein[item] = gened.data_dict[item]
    if item in geneinput :
        genein[item] = geneinput[item]
if genein['exbrate'] == -1111 :
    print('warning: exbrate value not specified')
    genein['exbrate'] = 0.0
if genein['pfsrate'] == -1111 :
    print('warning: pfsrate value not specified')
    genein['pfsrate'] = 0.0
#pprint(genein)
# map genein to tglf
tg = {}
# set the units flag to GENE
tg['UNITS'] = 'GENE'
# ky-grid
tg['KYGRID_MODEL'] = 0
tg['NKY'] = genein['nky0']
tg['KY'] = genein['nky0']*genein['kymin']
if tg['KY'] <= 0.0 :
    print('error: KY invalid ',KY)
#control switches
if(genein['nonlinear']=='True') :
    tg['USE_TRANSPORT_MODEL'] = 'T'
if genein['magn_geometry'] =='s-alpha' :
    tg['GEOMETRY_FLAG']= 0
if genein['magn_geometry'] == 'miller' :
    tg['GEOMETRY_FLAG']= 1
if genein['beta'] > 0.0 :
    tg['USE_BPER'] = 'T'
    tg['USE_BPAR'] = 'T'
if genein['bpar'] == 'False' :
    tg['USE_BPAR'] = 'F'

tg['USE_MHD_RULE'] = 'F'
if genein['dpdx_term'] == 'gradB_eq_curv' :
    tg['USE_MHD_RULE'] = 'T'
if genein['collision_op'] != 'pitch-angle' :
    print('warning: TGLF only has pitch-angle collisions')
# pprint(species[2]['name'])
tg['ZEFF'] = genein['zeff']
tg['VEXB_SHEAR'] = genein['exbrate']
# MILLER GEOMETRY
tg['RMIN_LOC'] = genein['trpeps']*genein['major_r']
tg['RMAJ_LOC'] = genein['major_r']
tg['ZMAJ_LOC'] = genein['major_z']
tg['DRMINDX_LOC'] = 1.0
tg['DRMAJDX_LOC'] = genein['drr']
tg['DZMAJDX_LOC'] = genein['drz']
tg['KAPPA_LOC'] = genein['kappa']
tg['S_KAPPA_LOC'] = genein['s_kappa']
tg['DELTA_LOC'] = genein['delta']
tg['S_DELTA_LOC'] = genein['s_delta']*(1.0-genein['delta']**2.0)**0.5
tg['ZETA_LOC'] = genein['zeta']
tg['S_ZETA_LOC'] = genein['s_zeta']
tg['Q_LOC'] = genein['q0']
if tg['Q_LOC'] < 0.0 :
    tg['Q_LOC'] = -tg['Q_LOC']
tg['Q_PRIME_LOC'] = genein['shat']*(tg['Q_LOC']/tg['RMIN_LOC'])**2
tg['P_PRIME_LOC'] = -genein['amhd']/(8.0*3.14159265*tg['Q_LOC']*tg['RMAJ_LOC']*tg['RMIN_LOC'])
tg['KX0_LOC'] = 0.0
#S-ALPHA GEOMETRY
tg['RMIN_SA'] = genein['trpeps']*genein['major_r']
tg['RMAJ_SA'] = genein['major_r']
tg['Q_SA'] = genein['q0']
if tg['Q_SA'] < 0.0 :
    tg['Q_SA'] = -tg['Q_SA']
tg['SHAT_SA'] = genein['shat']
tg['ALPHA_SA'] = genein['amhd']
tg['XWELL_SA'] = 0.0
tg['THETA0_SA'] = 0.0
tg['B_MODEL_SA'] = 1
tg['FT_MODEL_SA'] = 1
# plasma species data
nspecies = genein['n_spec']
nions = 0
for i in range(0, nspecies) :
    passive = 'F'
    if 'passive' in species[i] :
        if species[i]['passive'] == 'True' :
          passive = 'T'
    if species[i]['charge'] == -1:
# electrons
       if passive == 'T' :
          tg['ADIABATIC_ELEC'] = 'T'
       tg['ZS_1'] = -1.0
       tg['MASS_1'] = species[i]['mass']
       tg['RLNS_1'] = species[i]['omn']
       tg['RLTS_1'] = species[i]['omt']
       tg['TAUS_1'] = species[i]['temp']
       tg['AS_1'] = species[i]['dens']
       tg['VPAR_1'] = '0.0'
       tg['VPAR_SHEAR_1'] = genein['pfsrate']
    if species[i]['charge'] > 0:
        # skip passive species
       if passive == 'F' :
           nions = nions + 1
           if nions == 1 :
             tg['ZS_2'] = species[i]['charge']
             tg['MASS_2'] = species[i]['mass']
             tg['RLNS_2'] = species[i]['omn']
             tg['RLTS_2'] = species[i]['omt']
             tg['TAUS_2'] = species[i]['temp']
             tg['AS_2'] = species[i]['dens']
             tg['VPAR_2'] = 0.0
             tg['VPAR_SHEAR_2'] = genein['pfsrate']
           if nions == 2 :
             tg['ZS_3'] = species[i]['charge']
             tg['MASS_3'] = species[i]['mass']
             tg['RLNS_3'] = species[i]['omn']
             tg['RLTS_3'] = species[i]['omt']
             tg['TAUS_3'] = species[i]['temp']
             tg['AS_3'] = species[i]['dens']
             tg['VPAR_3'] = 0.0
             tg['VPAR_SHEAR_3'] = genein['pfsrate']
           if nions == 3 :
             tg['ZS_4'] = species[i]['charge']
             tg['MASS_4'] = species[i]['mass']
             tg['RLNS_4'] = species[i]['omn']
             tg['RLTS_4'] = species[i]['omt']
             tg['TAUS_4'] = species[i]['temp']
             tg['AS_4'] = species[i]['dens']
             tg['VPAR_4'] = 0.0
             tg['VPAR_SHEAR_4'] = genein['pfsrate']
           if nions == 4 :
             tg['ZS_5'] = species[i]['charge']
             tg['MASS_5'] = species[i]['mass']
             tg['RLNS_5'] = species[i]['omn']
             tg['RLTS_5'] = species[i]['omt']
             tg['TAUS_5'] = species[i]['temp']
             tg['AS_5'] = species[i]['dens']
             tg['VPAR_5'] = '0.0'
             tg['VPAR_SHEAR_5'] = genein['pfsrate']
           if nions == 5 :
             tg['ZS_6'] = species[i]['charge']
             tg['MASS_6'] = species[i]['mass']
             tg['RLNS_6'] = species[i]['omn']
             tg['RLTS_6'] = species[i]['omt']
             tg['TAUS_6'] = species[i]['temp']
             tg['AS_6'] = species[i]['dens']
             tg['VPAR_6'] = '0.0'
             tg['VPAR_SHEAR_6'] = genein['pfsrate']
           if nions > 5 :
               print('warning: TGLF stand alone version can only have 5 kinetic ion species')
               nions = 5
tg['NS'] = nions + 1
print('number of kinetic ions = ',nions)
# set TGLF parameters that need electron units
# note beta,debye2 and amhd are wrt Bref so TGLF will convert to B_unit internally
tg['BETAE'] = genein['beta']*tg['TAUS_1']*tg['AS_1']
tg['DEBYE'] = genein['debye2']**0.5
# note Cref=(Tref/mref)**0.5 the XNUE = nue*Lref/Cref in GENE units
tg['XNUE'] = genein['coll']*(4.0/tg['MASS_1']**0.5)*tg['AS_1']/tg['TAUS_1']**1.5
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
print('parameters_gene.txt succesfully converted to input.tglf')
