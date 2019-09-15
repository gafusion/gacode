"""profiles_generrorcheck.py contains routines that are useful to other python scripts to catch
errors in input.""" 

import sys
from profiles_gen.data import profiles_genData
import matplotlib.pyplot as plt
import math

#First are several functions that perform error checking that comes up frequently.
def opendir(directory):
    try:
        prof = profiles_genData(directory)
    except IOError:
        print(directory + " does not contain file input.profiles.  Type profiles_gen for help.")
        sys.exit()
    except IndexError:
        print("ERROR: Too few arguments.  Type profiles_gen for help.")
        sys.exit()
    return prof

def setdim(n):
    try:
        rval = int(n)
    except ValueError:
        print("ERROR: " + n + " is not a valid dimension.  Type profiles_gen for help.")
        sys.exit()
    except IndexError:
        print("ERROR: Please give second dimension.  Type profiles_gen for help.")
        sys.exit()
    return rval
