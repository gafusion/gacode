import numpy as np
import os
import sys


# ---- BEGIN input parsing ----

x = sys.argv[1]
var1 = x.split(',')[0]
var2 = x.split(',')[1] 

x = sys.argv[2]
min1 = float(x.split(',')[0])
min2 = float(x.split(',')[1])

x = sys.argv[3]
max1 = float(x.split(',')[0])
max2 = float(x.split(',')[1])

x = sys.argv[4]
n1 = int(x.split(',')[0])
n2 = int(x.split(',')[1])

run = int(sys.argv[5])

dir = sys.argv[6]

# ---- END input parsing ----

vec1 = min1+(max1-min1)*np.arange(n1)/float(n1-1)
vec2 = min2+(max2-min2)*np.arange(n2)/float(n2-1)

if run == 0:
    print(var1,vec1)
    print(var2,vec2)
    sys.exit()

# IMPORTANT: Set working directory
os.chdir(dir)

os.system('cp out.tglf.localdump input.tglf.template')

outfile = open('out.tglf.mrun_z','w')

root=os.environ['GACODE_ROOT']
for i2 in vec2:
    for i1 in vec1:
        a = var1+'='+str(i1)
        b = var2+'='+str(i2)
        os.system('cp input.tglf.template input.tglf')
        os.system('echo '+a+' >> input.tglf')
        os.system('echo '+b+' >> input.tglf')
        os.system('python '+root+'/tglf/bin/tglf_parse.py')
        os.system('tglf -p '+dir+' -e . > out')

        print('INFO: (tglf_mrun) '+a+' ; '+b)

        data = np.loadtxt('out.tglf.gbflux')
        ns   = np.loadtxt('out.tglf.grid')
        # Ge,Qe,Gi,Qi
        ni = int(ns[0])
        ovec = [0,1,ni,ni+1]

        for i in ovec:
            outfile.write(str(data[i])+' ')
        outfile.write('\n')

outfile.close()

# X:
np.savetxt('out.tglf.mrun_x',vec1)
# Y:
np.savetxt('out.tglf.mrun_y',vec2)

os.system('rm out')
os.system('mv input.tglf.template input.tglf')
