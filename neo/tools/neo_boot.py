import os
import numpy as np
import sys
import string

workdir = 'bdir'
tools   = os.environ['GACODE_ROOT']+'/neo/tools/'
harvestdata={}
jneo_harvest=[]
jsauter_harvest=[]
Ipsirho_harvest=[]

if len(sys.argv) < 22:
   print "python neo_boot.py <rmin> <q> <nuee> <ni1/ne> <zi1> <mi1/mD> <ti1/te> <zi2> <mi2/mD> <ti2/te> <delta> <kappa> <sdelta> <skappa> <zeta> <szeta> <shift> <zmagovera> <szmag> <shear> <betastar> <index> "
   sys.exit()

# EXAMPLE:
# python $GACODE_ROOT/neo/tools/neo_boot.py 0.17 2.0 0.1 0.9 1 1.0 1.0 6 6.0 1.0 0.1 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 (1992)

# In the input.neo, there are 3 species:
# electrons are species 1, main ions are species 2,
# and impurity ions are species 3
#
# Normalizations in the input.neo assumed to be:
# a = rmaj (the major radius), i.e. RMAJ_OVER_A=1.0
# T_norm = T_e (electron temperature), i.e. TEMP_1=1.0
# n_norm = n_e (electron density), i.e. DENS_1=1.0
# m_norm = m_deuterium
# v_norm = sqrt(T_norm/m_norm) = c_s (sound speed)

rmin   = sys.argv[1]  # r/a (Minor radius divided by minor radius of LCFS)
q      = sys.argv[2]  # safety factor
nuee   = sys.argv[3]  # electron collision frequency/(c_s/a)
ni1    = sys.argv[4]  # main ion density: n_i1/n_e

                       # (note: n_i2/n_e computed from quasi-neutrality)
zi1    = sys.argv[5]   # main ion charge (integer)
mi1    = sys.argv[6]   # main ion mass: m_i/m_deuterium
ti1    = sys.argv[7]   # main ion temperature: t_i/t_e
zi2    = sys.argv[8]   # impurity ion charge (integer)
mi2    = sys.argv[9]   # impurity ion mass: m_i2/m_deuterium
ti2    = sys.argv[10]  # impurity ion temperature: t_i2/t_e
delta  = sys.argv[11]  # triangularity
kappa  = sys.argv[12]  # elongation
sdelta = sys.argv[13]  # triangularity radial derivative (dimensionless)
skappa = sys.argv[14]  # elongation radial derivative (dimensionless)
zeta   = sys.argv[15]  # squareness
szeta  = sys.argv[16]  # squareness shear
shift  = sys.argv[17]  # shafranov shift
zmagovera = sys.argv[18] # normalized elevation
szmag  = sys.argv[19]  # gradient of elevation
shear  = sys.argv[20]  # magnetic shear
betastar = sys.argv[21] # normalized total beta NOTE: This parameter is not used in the standard kinetic equation calculation! But it is used in the case of an anisotropic temperature species

if len(sys.argv)==23 and sys.argv[22]!='None':
   harvestdata['IndexRS']=int(sys.argv[22])
harvestdata['rmin']=float(rmin)
harvestdata['q']=float(q)
harvestdata['nuee']=float(nuee)
harvestdata['ni1/ne']=float(ni1)
harvestdata['zi1']=float(zi1)
harvestdata['mi1/mD']=float(mi1)
harvestdata['ti1/te']=float(ti1)
harvestdata['zi2']=float(zi2)
harvestdata['mi2/mD']=float(mi2)
harvestdata['ti2/te']=float(ti2)
harvestdata['delta']=float(delta)
harvestdata['kappa']=float(kappa)
harvestdata['sdelta']=float(sdelta)
harvestdata['skappa']=float(skappa)
harvestdata['zeta']=float(zeta)
harvestdata['szeta']=float(szeta)
harvestdata['shift']=float(shift)
harvestdata['zmagovera']=float(zmagovera)
harvestdata['szmag']=float(szmag)
harvestdata['shear']=float(shear)
harvestdata['betastar']=float(betastar)

# Prepare simulation directory
os.system('rm -rf '+workdir)
os.system('mkdir '+workdir)
os.system('cp '+tools+'input.neo.neo_boot '+workdir+'/input.neo')

list = ['DLNNDR_1',
        'DLNNDR_2',
        'DLNNDR_3',
        'DLNTDR_1',
        'DLNTDR_2',
        'DLNTDR_3']

# Open input.neo, append parameters, close

neoin = open(workdir+'/input.neo','a')

# Set input: r_min
neoin.write('RMIN_OVER_A='+rmin+'\n')

# Set input: q
neoin.write('Q='+q+'\n')

# Set input: nu_ee/(cs/a)
neoin.write('NU_1='+nuee+'\n')

# Set input: delta
neoin.write('DELTA='+delta+'\n')

# Set input: kappa
neoin.write('KAPPA='+kappa+'\n')

# Set input: sdelta
neoin.write('S_DELTA='+sdelta+'\n')

# Set input: skappa
neoin.write('S_KAPPA='+skappa+'\n')

# Set input: zeta
neoin.write('ZETA='+zeta+'\n')

# Set input: szeta
neoin.write('S_ZETA='+szeta+'\n')

# Set input: shift
neoin.write('SHIFT='+shift+'\n')

# Set input: zmagovera
neoin.write('ZMAG_OVER_A='+zmagovera+'\n')

# Set input: szmag
neoin.write('S_ZMAG='+szmag+'\n')

# Set input: shear
neoin.write('SHEAR='+shear+'\n')

# Set input: betastar
neoin.write('BETA_STAR='+betastar+'\n')

# Set input: main ion charge, mass, temperature, density
neoin.write('Z_2='+zi1+'\n')
neoin.write('MASS_2='+mi1+'\n')
neoin.write('TEMP_2='+ti1+'\n')
neoin.write('DENS_2='+ni1+'\n')

# Set input: impurity ion charge, mass, temperature, density
neoin.write('Z_3='+zi2+'\n')
neoin.write('MASS_3='+mi2+'\n')
neoin.write('TEMP_3='+ti2+'\n')
# compute ni2/ne from quasi-neutrality
ni2 = (1.0-float(zi1)*float(ni1))/(1.0*float(zi2))
ni2 = str(ni2)
neoin.write('DENS_3='+ni2+'\n')

neoin.close()

z_all = [-1,float(zi1),float(zi2),-1,float(zi1),float(zi2)]
n_all = [1.0,float(ni1),float(ni2),1.0,float(ni1),float(ni2)]

cneo = []
csauter = []

# Run NEO
for i in range(6):
   neoin = open(workdir+'/input.neo','a')
   # Overlay gradients
   neoin.write('# '+str(i)+'\n')
   for j in range(6):
      if j == i:
         neoin.write(list[j]+'=1\n')
      else:
         neoin.write(list[j]+'=0\n')
   neoin.close()
   os.system('neo -e '+workdir)

   # Harvest output
   neoout = np.loadtxt(workdir+'/out.neo.transport')
   jneo=neoout[2]

   neoout = open(workdir+'/out.neo.diagnostic_geo','r')
   line = neoout.readlines()[0]
   ipsi = float(string.splitfields(line,'=')[1].rstrip())

   neoout = np.loadtxt(workdir+'/out.neo.equil')
   rhostar=neoout[3]

   neoout = np.loadtxt(workdir+'/out.neo.theory')
   jsauter=neoout[10]

   print 'jneo        '+str(jneo)
   print 'jsauter     '+str(jsauter)
   print 'I*Psi*rho_* '+str(ipsi*rhostar)

   jneo_harvest.append(jneo)
   jsauter_harvest.append(jsauter)
   Ipsirho_harvest.append(ipsi*rhostar)

   cneo.append(jneo/(ipsi*rhostar*abs(z_all[i])*n_all[i]))
   csauter.append(jsauter/(ipsi*rhostar*abs(z_all[i])*n_all[i]))
   
harvestdata['jneo']=jneo_harvest
harvestdata['jsauter']=jsauter_harvest
harvestdata['IPsirho']=Ipsirho_harvest
harvestdata['NEO_Coef']=cneo
harvestdata['SAUTER_Coef']=csauter

print list
print cneo
print csauter

sys.path.append(os.environ['GACODE_ROOT']+'/shared/harvest_client/')
from harvest_lib import harvest_send
import cPickle
with open('neo_boot.pkl','w') as f:
   cPickle.dump(harvestdata,f)

harvest_send(harvestdata,'Neo_boot',verbose=True,protocol='TCP',port=31000)
