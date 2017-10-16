#!/usr/bin/python
import sys


nels = (len(sys.argv)-1)/2
if (nels<1):
  sys.stderr.write("Usage scaling.py [fname title]+\n")
  sys.exit(1)

titles=[]
fnames=[]

for i in range(nels):
  fnames.append(sys.argv[1+i*2])
  titles.append(sys.argv[1+i*2+1])


#################

data={}

for i in range(nels):
  fname = fnames[i]
  title = titles[i]


  with open(fname,"r") as fd:
    lines = fd.readlines()

  for line in lines:
    if line[0] == '#':
      continue
    larr = line.split()
    tname = larr[0]
    nodes = int(larr[3])
    totals = float(larr[17])

    if not data.has_key(tname):
      data[tname]={}
    if not data[tname].has_key(nodes):
      data[tname][nodes]={}
    if not data[tname][nodes].has_key(i):
      data[tname][nodes][i]=totals
    elif (totals<data[tname][nodes][i]):
      data[tname][nodes][i] = totals


tlist = data.keys()
tlist.sort()
for tname in tlist:
  sys.stdout.write("#%s,nnodes"%tname)
  for title in titles:
    sys.stdout.write(",%s"%title)
  sys.stdout.write("\n")

  nlist = data[tname].keys()
  nlist.sort()
  for nodes in nlist:
    sys.stdout.write("%s,%s"%(tname,nodes))
    for i in range(nels):
      if data[tname][nodes].has_key(i):
        sys.stdout.write(",%6.1f"%data[tname][nodes][i])
      else:
        sys.stdout.write(",")
    sys.stdout.write("\n")

  sys.stdout.write("\n")



