import sys
import numpy as np
import matplotlib.pyplot as plt
import math 
from gacodefuncs import *
from cgyro.data import cgyrodata

# This is the fraction of the simulation length at the end 
# that is considered for the average
frac=0.5

# Get data directory from command line
data_dir=sys.argv[1]

# Read data using cgyro data class
data = cgyrodata(data_dir+'/')

# print "Time vector:"
# print data.t

# print
# print "Theta vector:"
# print data.theta

data.getgeo()

data.getdata()

ntheta=len(data.theta)
ntime=len(data.t)

tmin=math.ceil(ntime*(1-frac)) 
tlength=len(data.t[tmin:ntime])

rephiaver=math.fsum(map(math.fsum, data.phib[0,:,tmin:ntime-1]))/(ntheta*tlength)
reapaver=math.fsum(map(math.fsum, data.aparb[0,:,tmin:ntime-1]))/(ntheta*tlength)
rebpaver=math.fsum(map(math.fsum, data.bparb[0,:,tmin:ntime-1]))/(ntheta*tlength)

imphiaver=math.fsum(map(math.fsum, data.phib[1,:,tmin:ntime-1]))/(ntheta*tlength)
imapaver=math.fsum(map(math.fsum, data.aparb[1,:,tmin:ntime-1]))/(ntheta*tlength)
imbpaver=math.fsum(map(math.fsum, data.bparb[1,:,tmin:ntime-1]))/(ntheta*tlength)


print "Average Re[delta phi] and Im[delta phi] during time interval"
print rephiaver
print imphiaver
print "Average Re[delta A||] and Im[delta A||] during time interval"
print reapaver
print imapaver
print "Average Re[delta B||] and Im[delta B||] during time interval"
print rebpaver
print imbpaver

# writing time traces of poloidal averages into files 
# Time 

fphi= open("PhiAv.txt","w+")
fapa= open("ApaAv.txt","w+")
fbpa= open("BpaAv.txt","w+")
cosbpa= open("BpaCos.txt","w+")
rebpcosaver=0.0
for i in range(ntime):
     rephifsa=math.fsum(data.phib[0,:,i])/ntheta
     imphifsa=math.fsum(data.phib[1,:,i])/ntheta

     fphi.write("%11.3e"% (data.t[i]))
     fphi.write("%11.3e"% (rephifsa))
     fphi.write("%11.3e\n"% (imphifsa))

     reapafsa=math.fsum(data.aparb[0,:,i])/ntheta 
     imapafsa=math.fsum(data.aparb[1,:,i])/ntheta

     fapa.write("%11.3e"% (data.t[i]))
     fapa.write("%11.3e"% (reapafsa))
     fapa.write("%11.3e\n"% (imapafsa))

     rebpafsa=math.fsum(data.bparb[0,:,i])/ntheta 
     imbpafsa=math.fsum(data.bparb[1,:,i])/ntheta

     fbpa.write("%11.3e"% (data.t[i]))
     fbpa.write("%11.3e"% (rebpafsa))
     fbpa.write("%11.3e\n"% (imbpafsa))
     
     rebpacos=0.0
     imbpacos=0.0
     for j in range(ntheta):
        rebpacos=rebpacos+2.0*math.cos(data.theta[j])*data.bparb[0,j,i]/ntheta 
        imbpacos=imbpacos+2.0*math.cos(data.theta[j])*data.bparb[1,j,i]/ntheta

     cosbpa.write("%11.3e"% (data.t[i]))
     cosbpa.write("%11.3e"% (rebpacos))
     cosbpa.write("%11.3e\n"% (imbpacos))
     
     if i >= tmin:
         rebpcosaver=rebpcosaver+rebpacos/tlength

print "Average cos component of Re[delta B||] during time interval"
print rebpcosaver
      
fphi.close() 
fapa.close() 
fbpa.close() 
cosbpa.close()

# print "theta dependence of Re[delta B||] at last time point"
# print data.bparb[0,:,ntime-1]


