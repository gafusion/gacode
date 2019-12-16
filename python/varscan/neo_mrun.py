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

os.system('cp out.neo.localdump input.neo.template')

outfile = open('out.neo.mrun_z','w')

root=os.environ['GACODE_ROOT']
for i2 in vec2:
    for i1 in vec1:
        a = var1+'='+str(i1)
        b = var2+'='+str(i2)
        os.system('cp input.neo.template input.neo')
        os.system('echo '+a+' >> input.neo')
        os.system('echo '+b+' >> input.neo')
        os.system('python '+root+'/neo/bin/neo_parse.py')
        os.system('neo -p '+dir+' -e . > out')

        print('INFO: (neo_mrun) '+a+' ; '+b)

        data = np.loadtxt('out.neo.transport')

        equil = np.loadtxt('out.neo.equil')
        g_gb = equil[12]*equil[13]**1.5*equil[3]**2
        q_gb = equil[12]*equil[13]**2.5*equil[3]**2

        # Ge,Qe,Gi,Qi
        ovec = [13,14,5,6]
        data[13] = data[13]/g_gb
        data[14] = data[14]/q_gb
        data[5]  = data[5]/g_gb
        data[6]  = data[6]/q_gb
        for i in ovec:
            outfile.write(str(data[i])+' ')
        outfile.write('\n')

outfile.close()

# X:
np.savetxt('out.neo.mrun_x',vec1)
# Y:
np.savetxt('out.neo.mrun_y',vec2)

os.system('rm out')
os.system('mv input.neo.template input.neo')

