# Example usage:
# python template.py <simdir>
#
# see gacode/f2py/pygacode/cgyro/data.py for available data objects

import sys
import numpy as np
import matplotlib.pyplot as plt
from ..gacodefuncs import *
from .data import cgyrodata

# Get data directory from command line
data_dir=sys.argv[1]

# Read data using cgyro data class
data = cgyrodata(data_dir+'/')

print("Time vector:")
print(data.t)

print()
print("Theta vector:")
print(data.theta)

data.getgeo()
print()
print("B(theta):")
print(data.geo[:,2])

print()
print("Re[phi(theta)] at last time:")
print(data.phib[0,:,-1])
print()
print("Im[phi(theta)] at last time:")
print(data.phib[1,:,-1])
