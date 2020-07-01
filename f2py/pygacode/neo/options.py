from .data import NEOData
import sys
import numpy as np

sim1 = NEOData(sys.argv[1])
keys = []

for k1, v1 in sim1.transport.items():
    s = 0
    s = np.sum(v1.data)
    if s != 0:
#Unless every entry in v2 is zero, it gets added to the list.
        keys.append(k1)
keys.sort()
for k in keys:
    print(k)
