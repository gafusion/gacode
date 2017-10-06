#!/usr/bin/python
import sys

f1=sys.argv[1]
f2=sys.argv[2]

with open(f1,"r") as fd:
  lines1 = fd.readlines()

with open(f2,"r") as fd:
  lines2 = fd.readlines()

data2={}
data={}
data1={} # only rejected ones

for line in lines2:
  if line[0] == '#':
    sys.stdout.write(line)
    continue
  larr = line.split()
  tname = larr[0]
  nodes = int(larr[3])

  if not data.has_key(tname):
    data[tname]={}
    data2[tname]={}
  if not data[tname].has_key(nodes):
    data[tname][nodes]=[]
    data2[tname][nodes]=[]
  data[tname][nodes].append(line)
  data2[tname][nodes].append(line)

for line in lines1:
  if line[0] == '#':
    continue  # comment already dealt with above
  larr = line.split()
  tname = larr[0]
  ltype = larr[1]
  nodes = int(larr[3])

  if (ltype!="noHT") and data2.has_key(tname) and data2[tname].has_key(nodes):
    if not data1.has_key(tname):
      data1[tname]={}
    if not data1[tname].has_key(nodes):
      data1[tname][nodes]=[]
    data1[tname][nodes].append(line)
    continue # this already dealt with above

  if not data.has_key(tname):
    data[tname]={}
  if not data[tname].has_key(nodes):
    data[tname][nodes]=[]
  data[tname][nodes].append(line)

tlist = data.keys()
tlist.sort()
for tname in tlist:
  nlist = data[tname].keys()
  nlist.sort()
  for nodes in nlist:
    for line in data[tname][nodes]:
      sys.stdout.write(line)

sys.stdout.write("#Rejected data from old run\n")
tlist = data1.keys()
tlist.sort()
for tname in tlist:
  nlist = data1[tname].keys()
  nlist.sort()
  for nodes in nlist:
    for line in data1[tname][nodes]:
      sys.stdout.write("#%s"%line)



