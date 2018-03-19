#!/usr/bin/python

import sys

# first convert the file into internal structure
# tests[test] - list of results, each a list of columns
tests={}

with open(sys.argv[1],'r') as fd:
    lines=fd.readlines()

for line in lines:
    els = line.split()
    if (len(els)!=12):
      print "WARNING Invalid line found: %s"%line
      continue
    test=els[0]
    #we will ignore the index value
    if not tests.has_key(test):
        tests[test]=[]
    tests[test].append(els[1:])

# extract relevant values
raws={}
restartios={}

for test in tests.keys():
      raws[test]={}
      raws[test]={'sums':[],'mid':[]}
      cel = raws[test]
      for i in range(11):
          cel['sums'].append(0.0)
          cel['mid'].append(0.0)
          
      for idx in (0,1,3,4):
          els = tests[test][idx]
          for i in range(11):
              el=float(els[i])
              elsum=cel['sums'][i]
              elsum += el
              cel['sums'][i]=elsum
      for idx in [2]:
          els = tests[test][idx]
          for i in range(11):
              el=float(els[i])
              cel['mid'][i] = el

testlist=raws.keys()
testlist.sort()
for test in testlist:
        cel = raws[test]
        meantotal = cel['sums'][10]/4
        meanio = cel['sums'][9]/4
        restartio = cel['mid'][9] - meanio
        if (restartio<0.0):
           restartio = 0.0

        restartios[test] = (meanio,restartio,meantotal)


# compute means
mtests={}
computed={}

for rawtest in restartios.keys():
    test=rawtest.rsplit("_",2)[0]
    nmpi=int(rawtest.rsplit("_",2)[1])
    #we will ignore the index value
    if not mtests.has_key(test):
        mtests[test]={}
    if not mtests[test].has_key(nmpi):
        mtests[test][nmpi] = []
    mtests[test][nmpi].append(restartios[rawtest])

for test in mtests.keys():
    computed[test]={}
    for nmpi in mtests[test].keys():
      computed[test][nmpi]={'count':len(mtests[test][nmpi]),
                            'mins': [],'maxs':[],'sums':[]}
      cel = computed[test][nmpi]
      for i in range(3):
          cel['mins'].append(1.e99)
          cel['maxs'].append(0.0)
          cel['sums'].append(0.0)

      for els in mtests[test][nmpi]:
          for i in range(3):
              el=float(els[i])
              elmin=cel['mins'][i]
              elmax=cel['maxs'][i]
              elsum=cel['sums'][i]
              if (el<elmin): elmin=el
              if (el>elmax): elmax = el
              elsum += el
              cel['mins'][i]=elmin
              cel['maxs'][i]=elmax
              cel['sums'][i]=elsum

# now print out results
titles=['base_io','restart_io','TOTAL']

fdout = open("%s.rstrt"%sys.argv[1],'w')

fdout.write("#Absolute walclock values\n")

outstr="#test nmpi nnodes "

for i in range(3):
    outstr += " mean_%s"%titles[i]
for i in range(3):
    outstr += " mindiff_%s"%titles[i]
for i in range(3):
    outstr += " maxdiff_%s"%titles[i]
fdout.write("%s\n"%outstr)

testlist=computed.keys()
testlist.sort()
for test in testlist:
    ndiv = 64
    ompstrarr= test.rsplit("_",1)
    if (len(ompstrarr)>1):
      ompstr = ompstrarr[1]
      if (ompstr[0]=='o'): # _oXX
        ndiv /= int(ompstr[1:])/2
    mpilist=computed[test].keys()
    mpilist.sort()
    for nmpi in mpilist:
        nnodes = nmpi/ndiv
        cel = computed[test][nmpi]
        outstr=""
        means={}
        for i in range(3):
            means[i] = cel['sums'][i]/cel["count"]
            outstr += " %.4e"%means[i]
        for i in range(3):
            dperc = 0.0
            if (means[i]>=1.e-2):
               diff = cel['mins'][i] - means[i]
               dperc = diff/means[i]
            outstr += " %i%%"%int(dperc*100)
        for i in range(3):
            dperc = 0.0
            if (means[i]>=1.e-2):
               diff = cel['maxs'][i] - means[i]
               dperc = diff/means[i]
            outstr += " %i%%"%int(dperc*100)

        fdout.write("%s %i %i %s\n"%(test,nmpi,nnodes, outstr))

fdout.close()

