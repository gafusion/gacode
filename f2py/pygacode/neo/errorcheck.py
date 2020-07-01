import sys
from .data import NEOData
import matplotlib.pyplot as plt
import math

#First are several functions that perform error checking that comes up frequently.
def opendir(directory):
    try:
        sim = NEOData(directory)
    except IOError:
        print
