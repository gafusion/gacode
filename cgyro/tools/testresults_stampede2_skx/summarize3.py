#!/usr/bin/python

import sys

# first convert the file into internal structure
# tests[test][nmpi] - list of results, each a list of columns
tests={}

with open(sys.argv[1],'r') as fd:
    lines=fd.readlines()

for line in lines:
    els = line.split()
    if (len(els)!=12):
      print "WARNING Invalid line found: %s"%line
      continue
    test=els[0].rsplit("_",2)[0]
    nmpi=int(els[0].rsplit("_",2)[1])
    #we will ignore the index value
    if not tests.has_key(test):
        tests[test]={}
    if not tests[test].has_key(nmpi):
        tests[test][nmpi] = []
    tests[test][nmpi].append(els[1:])

# now calc values
computed={}

for test in tests.keys():
    computed[test]={}
    for nmpi in tests[test].keys():
      computed[test][nmpi]={'count':len(tests[test][nmpi]),
                            'mins': [],'maxs':[],'sums':[]}
      cel = computed[test][nmpi]
      for i in range(11):
          cel['mins'].append(1.e99)
          cel['maxs'].append(0.0)
          cel['sums'].append(0.0)
          
      for els in tests[test][nmpi]:
          for i in range(11):
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
titles=['str','str_comm','nl','nl_comm','field','field_comm','shear','coll','coll_comm','io','TOTAL','prop']

fdout = open("%s.abs"%sys.argv[1],'w')
fdout2 = open("%s.scaled"%sys.argv[1],'w')
fdout3 = open("%s.summary"%sys.argv[1],'w')

fdout.write("#Absolute walclock values\n")
fdout2.write("#Scaled walclock values (by nnodes)\n")
fdout3.write("#Scaled walclock values (by nnodes)\n")


outstr="#test variant nmpi nnodes nprocsxnode nomp mpiorder"

for i in range(11):
    outstr += " mean_%s"%titles[i]
for i in range(11):
    outstr += " mindiff_%s"%titles[i]
for i in range(11):
    outstr += " maxdiff_%s"%titles[i]
fdout.write("%s\n"%outstr)
fdout2.write("%s\n"%outstr)
fdout3.write("%s\n"%outstr)

coresxnode=96

testlist=computed.keys()
testlist.sort()
for test in testlist:
    tarr = test.split("_",1)
    btest = tarr[0]
    vtest = "noHT"
    if (len(tarr)>1):
      vtest = tarr[1]
    ndiv = coresxnode
    ompstrarr= test.rsplit("_",1)
    if (len(ompstrarr)>1):
      ompstr = ompstrarr[1]
      if (ompstr[0]=='o'): # _oXX
        ndiv /= int(ompstr[1:])
      else: # try _oxx_r2
        ompstrarr= ompstrarr[0].rsplit("_",1)
        if (len(ompstrarr)>1):
           ompstr = ompstrarr[1]
           if (ompstr[0]=='o'): # _oXX
             ndiv /= int(ompstr[1:])

    mpiorder=1
    if vtest.find("_r2")!=-1:
      mpiorder=2
    mpilist=computed[test].keys()
    mpilist.sort()
    for nmpi in mpilist:
        nnodes = nmpi/ndiv
        cel = computed[test][nmpi]
        outstr=""
        outstr2=""
        means={}
        for i in range(11):
            means[i] = cel['sums'][i]/cel["count"]
            outstr += " %.4e"%means[i]
            outstr2 += " %.4e"%(means[i]*nnodes)
        for i in range(11):
            dperc = 0.0
            if (means[i]>=1.e-2):
               diff = cel['mins'][i] - means[i]
               dperc = diff/means[i]
            outstr += " %i%%"%int(dperc*100)
            outstr2 += " %i%%"%int(dperc*100)
        for i in range(11):
            dperc = 0.0
            if (means[i]>=1.e-2):
               diff = cel['maxs'][i] - means[i]
               dperc = diff/means[i]
            outstr += " %i%%"%int(dperc*100)
            outstr2 += " %i%%"%int(dperc*100)

        nprocs = nmpi/nnodes
        nomp = coresxnode/nprocs
        fdout.write( "%s %s %i %i %i %i %i %s\n"%(btest,vtest,nmpi,nnodes, nprocs, nomp, mpiorder, outstr))
        fdout2.write("%s %s %i %i %i %i %i %s\n"%(btest,vtest,nmpi,nnodes, nprocs, nomp, mpiorder, outstr2))

# here we order by base test name only
btestdict={}
for test in testlist:
    btest = test.split("_",1)[0]
    if btest not in btestdict.keys():
        btestdict[btest] = {}
    ndiv = coresxnode
    ompstrarr= test.rsplit("_",1)
    if (len(ompstrarr)>1):
      ompstr = ompstrarr[1]
      if (ompstr[0]=='o'): # _oXX
        ndiv /= int(ompstr[1:])
      else: # try _oxx_r2
        ompstrarr= ompstrarr[0].rsplit("_",1)
        if (len(ompstrarr)>1):
           ompstr = ompstrarr[1]
           if (ompstr[0]=='o'): # _oXX
             ndiv /= int(ompstr[1:])
    mpilist=computed[test].keys()
    for nmpi in mpilist:
        nnodes = nmpi/ndiv
        if nnodes not in btestdict[btest].keys():
           btestdict[btest][nnodes] = []
        btestdict[btest][nnodes].append((test,nmpi))

btests = btestdict.keys()
btests.sort()
for btest in btests:
  nnodeslist = btestdict[btest].keys()
  nnodeslist.sort()
  for nnodes in nnodeslist:
    mpilist=btestdict[btest][nnodes]
    mpilist.sort()
    for (test,nmpi) in mpilist:
        tarr = test.split("_",1)
        vtest = "noHT"
        if (len(tarr)>1):
           vtest = tarr[1]
        mpiorder=1
        if vtest.find("_r2")!=-1:
           mpiorder=2
        cel = computed[test][nmpi]
        outstr=""
        outstr2=""
        means={}
        for i in range(11):
            means[i] = cel['sums'][i]/cel["count"]
            outstr += " %.4e"%means[i]
            outstr2 += " %.4e"%(means[i]*nnodes)
        for i in range(11):
            dperc = 0.0
            if (means[i]>=1.e-2):
               diff = cel['mins'][i] - means[i]
               dperc = diff/means[i]
            outstr += " %i%%"%int(dperc*100)
            outstr2 += " %i%%"%int(dperc*100)
        for i in range(11):
            dperc = 0.0
            if (means[i]>=1.e-2):
               diff = cel['maxs'][i] - means[i]
               dperc = diff/means[i]
            outstr += " %i%%"%int(dperc*100)
            outstr2 += " %i%%"%int(dperc*100)

        nprocs = nmpi/nnodes
        nomp = coresxnode/nprocs
        fdout3.write("%s %s %i %i %i %i %i %s\n"%(btest,vtest,nmpi,nnodes, nprocs, nomp, mpiorder, outstr2))


fdout.close()
fdout2.close()
fdout3.close()

