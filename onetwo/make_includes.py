#!python
print 'Starting make_includes.py'
import os
import sys
import glob
path = []
objs = []
stray = []
l = stray
for i in sys.argv[1:]:
  if i=='-obj':
    l = objs
    continue
  if i=='-path':
    l = path
    continue
  l.append(i)
if len(stray):
  raise Exception('Files were given before -obj or -path were given')
exp_path = []
for v in path:
  exp_path.extend(v.split(':'))
if '.' or './' not in exp_path:
  exp_path.insert(0,'./')
srcs = []
for p in exp_path:
  for ext in ['.f','.f90','.F','.F90']:
    srcs.extend([os.path.split(f)[1] for f in glob.glob(p+'/*%s'%ext)])
src_base = map(lambda x: x.split('.')[:-1],srcs)
obj_base = map(lambda x: x.split('.')[:-1],objs)
make_include = []
for oi,o in enumerate(objs):
  if obj_base[oi] in src_base:
    src = srcs[src_base.index(obj_base[oi])]
    make_include.append('include .dep_incs/%s.inc'%src)
make_include_s = '\n'.join(make_include)
if os.path.exists('make.includes'):
  with open('make.includes','r') as f:
    orig=f.read()
  if orig!=make_include_s:
    print 'writing new make.includes'
    with open('make.includes','w') as f:
      f.write(make_include_s)
else:
  with open('make.includes','w') as f:
    f.write(make_include_s)
  
print 'finished make_includes.py'
