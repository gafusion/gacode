import os
import numpy as np
import sys
import string

workdir = 'bdir'
tools   = os.environ['GACODE_ROOT']+'/neo/tools/'

if len(sys.argv) < 11:
   print "python neo_boot.py <rmin> <q> <nuee> <ni1/ne> <ti1/te> <shift> <kappa> <skappa> <delta> <sdelta>"
   sys.exit()

# EXAMPLE:
# python $GACODE_ROOT/neo/tools/neo_boot.py 0.17 2.0 0.1 0.9 1.0 0.0 1.0 0.0 0.0 0.0

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

rmin   = sys.argv[1]   # r/a (Minor radius divided by minor radius of LCFS) 
q      = sys.argv[2]   # safety factor
nuee   = sys.argv[3]   # electron collision frequency/(c_s/a)
ni1    = sys.argv[4]   # main ion density: n_i1/n_e
ti1    = sys.argv[5]   # main ion temperature: t_i/t_e
shift  = sys.argv[6]   # Shafranov shift: dR/dr 
kappa  = sys.argv[7]   # elongation (dimensionless)
skappa = sys.argv[8]   # elongation radial derivative(dimensionless)
delta  = sys.argv[9]   # triangularity dimensionless)
sdelta = sys.argv[10]  # triangularity radial derivative (dimensionless)

zi1    =  str(1.0)   # main ion charge (integer)
mi1    =  str(1.0)   # main ion mass: m_i/m_deuterium
zi2    =  str(6.0)   # impurity ion charge (integer)
mi2    =  str(6.0)   # impurity ion mass: m_i2/m_deuterium
ti2    =  ti1        # impurity ion temperature: t_i2/t_e

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

# Set input: shift
neoin.write('SHIFT='+shift+'\n')

# Set input: kappa
neoin.write('KAPPA='+kappa+'\n')

# Set input: skappa
neoin.write('S_KAPPA='+skappa+'\n')

# Set input: delta
neoin.write('DELTA='+delta+'\n')

# Set input: sdelta
neoin.write('S_DELTA='+sdelta+'\n')

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

   #neoout = open(workdir+'/out.neo.diagnostic_geo2','r') 
   #line = neoout.readlines()[0]
   #ipsi = float(string.splitfields(line,'=')[1].rstrip())
   neoout = np.loadtxt(workdir+'/out.neo.diagnostic_geo2') 
   ipsi   = neoout[0]
   geo1   = neoout[1]
   geo2   = neoout[2]
   geo3   = neoout[3]
   geo4   = neoout[4]
   geo5   = neoout[5]
   geo6   = neoout[6]
   geo7   = neoout[7]
   geo8   = neoout[8]
   geo9   = neoout[9]
   geo10  = neoout[10]

   neoout = np.loadtxt(workdir+'/out.neo.equil') 
   rhostar=neoout[3]
   
   #print 'jneo        '+str(jneo)

   cneo.append(jneo/(ipsi*rhostar*abs(z_all[i])*n_all[i]))

print list
print cneo
print geo1,geo2,geo3
