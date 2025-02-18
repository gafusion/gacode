import numpy as np
import numpy.polynomial.chebyshev as pc

# list should start with W and end with H
ion_vec = ['W','He','Be','C','O','N','F','Ne','Al','Si','Ar','Ca','Fe','Ni','Kr','Mo','Xe','Li','H']

# Temperature limits in keV
te_max = 50.0
te_min = 0.05

def radfmt(i,c,lang):
    
    if lang == 'f90':
        u = """'"""+i+"""'"""
        v = """'D'"""
        w = """'T'"""
        if i == 'W':
            rstr = '  if (name == '+u+') then\n'
        elif i == 'H':
            rstr = '  else if (name == '+u+') or (name == '+v+') or (name == '+w+') then\n'
        else:
            rstr = '  else if (name == '+u+') then\n'
        
        rstr += '     c(:) = (/{:+.12e},{:+.12e},{:+.12e},{:+.12e},&\n'.format(c[0],c[1],c[2],c[3])
        rstr += '              {:+.12e},{:+.12e},{:+.12e},{:+.12e},&\n'.format(c[4],c[5],c[6],c[7])
        rstr += '              {:+.12e},{:+.12e},{:+.12e},{:+.12e}/)\n'.format(c[8],c[9],c[10],c[11])
    elif lang == 'julia':
        u = '''"'''+i+'''"'''
        rstr  = '    elseif name == '+u+'\n'
        rstr += '        coefficients = [\n'
        for j in range(11):
            rstr += '            {:+.12e},\n'.format(c[j])
        rstr += '            {:+.12e}\n'.format(c[j])
        rstr += '        ]\n'

    return rstr

# log-temperature bounds (convert to eV)
lte1 = np.log(te_max*1e3)
lte0 = np.log(te_min*1e3)
dt   = lte1-lte0
nc   = 12

x = np.loadtxt('aurora/x.txt')
t = np.loadtxt('aurora/te.txt')

fstr = ''
jstr = ''
for i in ion_vec:
    y = np.loadtxt('aurora/'+i+'.txt')
    # 1MW/cm^3 = 10^13 erg/cm^3/s 
    c = pc.chebfit(x,np.log(y*10**13),nc)
    rstr = radfmt(i,c,'julia')
    jstr = jstr+rstr
    rstr = radfmt(i,c,'f90')
    fstr = fstr+rstr
   
myfile = open("julia.txt","w")
myfile.write(jstr)
myfile.close()

myfile = open("f90.txt","w")
myfile.write(fstr)
myfile.close()
