import sys
from profiles_gen.data import profiles_gen_format

x = profiles_gen_format()

rows = len(x.tags)
cols = len(x.tags[0])

print("------------------------------------------------------------")
print("input.profiles format")
print("------------------------------------------------------------")


head = "Col\ttag\tUnits"
fmt = "{idx:d}\t{v1:s}\t{v2:s}"

for i in range(rows):
    print('\nBlock '+str(i+1)+' -------------------------------------------------------------\n')
    print("Col".ljust(5)+'Tag'.ljust(20)+'Units'.ljust(15),'LaTeX')
    for j in range(cols):
        u = x.tags[i][j]
        print(str(j+1).ljust(5)+x.desc[u][2].ljust(20)+x.desc[u][1].ljust(15)+x.desc[u][0])

for i in range(rows):
    for j in range(cols):
        u = x.tags[i][j]
        if u in list(x.defs.keys()):
            print('-------------------------------------------------------------------------------')
            print('\nBlock '+ str(i+1)+'  Column '+str(j+1)+'\n\nTag: '+x.desc[u][2])
            print()
            try:
                print(x.defs[u])
                print('-------------------------------------------------------------------------------')
            except:
                pass
