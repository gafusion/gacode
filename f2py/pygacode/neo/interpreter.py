from .data import NEOData
import sys
import matplotlib.pyplot as plt
import math

args = ['GAMMA', 'Q', '1upar', 'jboot']

sim1 = NEOData(sys.argv[1])
legend = bool(int(sys.argv[4]))
verbose = bool(int(sys.argv[5]))

try:
    n1 = int(sys.argv[2])
except ValueError:
    print(sys.argv[2] + " is not a vaild dimension.")
    sys.exit()
try:
    n2 = int(sys.argv[3])
except ValueError:
    print(sys.argv[3] + " is not a vaild dimension.")
    sys.exit()

if len(sys.argv) < 7:
    for arg in args:
        sim1.plot(arg, n1=n1, n2=n2, legend=legend, verbose=verbose)
    plt.show()
else:
    args = sys.argv[6:]
    for arg in args:
        sim1.plot(arg, n1=n1, n2=n2, legend=legend, verbose=verbose)
    plt.show()
