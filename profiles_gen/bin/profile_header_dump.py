import sys
import string

inputfile = sys.argv[1]
outfile = open('profile_header','w')

index = 0
for line in open(inputfile,'r').readlines():
   if index > 0:
      index = index+1
      line = string.strip(line)
      if len(line) > 1:
         x = line.split()
         outfile.write(x[2]+' '+x[3]+' '+x[4]+'\n')
      else:
         sys.exit()

   if 'IONS' in line:
      index = 1


