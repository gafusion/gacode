import os

mydir=[]
mydir.append("1_iterdb")
mydir.append("../2_iterdbnc")
mydir.append("../3_plasmastate")
mydir.append("../4_genf")
mydir.append("../6_pfile")
mydir.append("../8_null")

mystr=[]
mystr.append("-i iterdb_kinefit_128913.01500 -g g128913.01500 -cer cer141459.03890")
mystr.append("-i iterdb.nc.04125 -g g141397.04125_kinetic")
mystr.append("-i p133444.cdf -g p133444.geq")
mystr.append("-i pi3-801-0010-prof11.csv -g pi3-801-0010-g-eqdsk.txt")
mystr.append("-i p132010.003058.peq -g g132010.003058")
mystr.append("-i null -g g128913.01500")

for i in range(len(mydir)):
   print(mydir[i])
   os.chdir(mydir[i])
   os.system('profiles_gen '+mystr[i]+' > out')
   os.system('sed -i 1,2d input.gacode')
   os.system('diff -s input.gacode input.gacode.check')
   print()

