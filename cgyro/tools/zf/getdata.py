# file processed by 2to3
from __future__ import print_function, absolute_import
from builtins import map, filter, range
import sys
import numpy as np
import matplotlib.pyplot as plt
from gacodefuncs import *
from cgyro.data import cgyrodata

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

data.getdata()
print()
print("B_unit:")
print(data.b_unit)

print()
print("Re[phi(theta)] at last time:")
print(data.phib[0,:,-1])
print()
print("Im[phi(theta)] at last time:")
print(data.phib[1,:,-1])
