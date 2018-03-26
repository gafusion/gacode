import os
import numpy as np
import sys
from scipy.interpolate import Rbf
import tarfile

d = {}

# Inputs
#  indata(1) = neo_rmin_over_a_in
#  indata(2) = neo_q_in
#  indata(3) = log10(neo_nu_1_in)
#  indata(4) = neo_dens_in(2)
#  indata(5) = neo_temp_in(2)
#----
#  indata(6) = neo_delta_in
#  indata(7) = neo_s_delta_in
#  indata(8) = neo_s_kappa_in

fulltag_in = ['in1_eps','in2_q','in3_nu','in4_ni1','in5_ti1',
              'in7_delta','in8_sdelta','in9_kappa','in10_skappa']

tar = tarfile.open("circle3.tar.gz")
root='circle3/'

indata = np.fromfile(root+'indata.dat',dtype='float',sep=" ")
n = len(indata)/9
indata = np.reshape(indata,(9,n),'F')

tag_in = {}
# Basic 5 inputs
for i in range(5):
   d[fulltag_in[i]] = indata[i,:]
   tag_in[i] = fulltag_in[i]

# Geometry inputs
#  ingeodata(0)  = I/psiprime
#  ingeodata(1)  = ftrap
#  ingeodata(2)  = <B^2>
#  ingeodata(3)  = <1/B^2>
#  ingeodata(4)  = <(b dot grad B)^2>
#  ingeodata(5)  = <(b dot grad B)^2/B^2>
#  ingeodata(6)  = <q R0 / Jpsi B>
#  ingeodata(7)  = <Btor^2>
#  ingeodata(8)  = <Bpol^2>
#  ingeodata(9)  = <|grad r|^2>
#  ingeodata(10) = <|grad r|^2/B^2>
#  ingeodata(11) = <-grad_r * gsin/B>

ingeodata = np.fromfile(root+'ingeodata.dat',dtype='float',sep=" ")
ingeodata = np.reshape(ingeodata,(12,n),'F')

# Magic 6th parameter
d['in6_geo'] = ingeodata[1,:]
tag_in[5] = 'in6_geo'

# Outputs
#  outdata(1) = Cne
#  outdata(2) = CTe
#  outdata(3) = Cni1
#  outdata(4) = CTi1
#  outdata(5) = Cni2
#  outdata(6) = CTi2

tag_out = ['cne','cte','cni1','cti1','cni2','cti2']

outdata = np.fromfile(root+'outdata.dat',dtype='float',sep=" ")
outdata = np.reshape(outdata,(6,n),'F')

# 6 outputs
for i in range(6):
   d[tag_out[i]] = outdata[i,:]

for key in sorted(d.keys()):
#   print key,d[key][0]
   print key
   
in_mins = np.zeros(6)
in_maxs = np.zeros(6)
dscale = {}
for i in range(6):
   in_mins[i] = min(d[tag_in[i]])
   in_maxs[i] = max(d[tag_in[i]])
   if (in_mins[i] != in_maxs[i]):
      dscale[tag_in[i]] = (d[tag_in[i]]-in_mins[i])/(in_maxs[i]-in_mins[i])
   else:
      dscale[tag_in[i]] =  d[tag_in[i]]   
      
print 'in_mins: ', in_mins
print 'in_maxs: ', in_maxs

# Test point
eps=0.3391239
q=4.4351
lnu=np.log10(0.798065)
ni1=0.741015
ti1=1.564186
geo=0.757284526

n1 = 6
n2 = 4
n3 = 4
n4 = 3
n5 = 3
n6 = 6

myin = [eps,q,lnu,ni1,ti1,geo]
myinscale = np.zeros(6)
for i in range(6):
   if (in_mins[i] != in_maxs[i]):
      myinscale[i] = (myin[i]-in_mins[i])/(in_maxs[i]-in_mins[i])
   else:
      myinscale[i] = myin[i]

print 'myin: ', myin      
print 'myinscale: ', myinscale

for x in tag_out:
   rbf = Rbf(dscale['in1_eps']*n1,dscale['in2_q']*n2,dscale['in3_nu']*n3,
             dscale['in4_ni1']*n4,dscale['in5_ti1']*n5,dscale['in6_geo']*n6,
             d[x],function='cubic',epsilon=1.0)

   np.savetxt('node_'+x,rbf.nodes,fmt='%.12e',newline=os.linesep)
 
   print x,rbf(myinscale[0]*n1,myinscale[1]*n2,myinscale[2]*n3,
               myinscale[3]*n4,myinscale[4]*n5,myinscale[5]*n6)

