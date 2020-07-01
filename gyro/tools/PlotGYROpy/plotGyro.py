#!/usr/bin/env python
# file processed by 2to3
from __future__ import print_function, absolute_import
from builtins import map, filter, range

import numpy
import optparse
import readGyroH5
from matplotlib import pylab

def makeContour(R,Z,data):
     """ Make a contour plot """ 
     pylab.figure()
     CS2=pylab.contourf(R,Z,data,cmap=pylab.cm.get_cmap('cool'))
     #pylab.title('High resolution contour of density')
     pylab.show()
     return 

def main():
    """ SEK: Once the pretty print can be figured out, then we can look
    at completing the standalone executable
    """
    parser = optparse.OptionParser(usage="%prog [options] inputFile")
    parser.add_option('-i', '--input', dest='input',
                      help='Name of input file if not specified as argument.',
                      default='')
    parser.add_option('-v', '--variable', dest='variable',
                      help='Variable to plot',
                      default='')

    options, args = parser.parse_args()

    # Process arguments
    if len(args) > 1:
      parser.print_usage()
      return
    elif len(args)==1:
      inputFile=args[0]
    else:
      if options.input == '':
        print("Must specify an input file")
        return
      else:
        inputFile=options.input
    
    return

if __name__ == "__main__":
        main()
