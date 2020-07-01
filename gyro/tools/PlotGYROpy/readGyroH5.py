#!/usr/bin/env python
# file processed by 2to3
from __future__ import print_function, absolute_import
from builtins import map, filter, range

import numpy
import optparse
import tables
import sys

def readProfile(h5file,nodeName=""):
  """ given an gyro HDF5 file that contains the profile data
      return a dictionary containing the data.
      Assume that the profile data does not have groups
  """

  hf=tables.openFile(h5file,'r')
  dataNames=list(hf.root._v_children.keys())
  # Idiot check done by seeing if n_kinetic is there
  if dataNames.count("n_kinetic") == 0:
     print(h5file+" does not seem to contain the gyro profile data")
     return
  data={}
  for dName in dataNames:
    getNode="hf.root."+dName+".read()"
    try:
      data[dName]=eval(getNode)
    except:
      print("Problems getting data: "+dName)
  hf.close()
  return data

def printData(dataDict,outfile=""):
  """ Given a dictionary, pretty print to either a file if specified 
      or stdout if not
  """
  if outfile:
    of=open(outfile,"w")
  else:
    of=sys.stdout

  for dName in list(dataDict.keys()):
    of.write(dName+'\n')
    of.write(dataDict[dName])
    sys.stdout.write(dataDict[dName])
    of.write('\n')
    of.write('-------------------------------------------------------\n')
  of.close()
  return 

def readToroidal(h5file,nodeName=""):
  """ For the fine files, the toroidal slices are placed into groups
      to facilitate how it is plotted in visit.  This reads those in
  """

  hf=tables.openFile(h5file,'r')
  torlist=[]
  for  hfnode in list(hf.root._v_children.keys()):
    if hfnode.endswith("_toroidal"):
      torlist.append(hfnode)
  if len(torlist) == 0:
    print(h5file+" has no mode data.")
    return
  torData={}
  if nodeName:
    if torlist.count(nodeName) > 0:
      torlist=[nodeName]
    else:
      print(nodeName +" not found in file.")
      return
  for hfnode in torlist:
    childNodes="hf.root."+hfnode+"._v_children"
    for chNode in eval(childNodes):
      getNode="hf.root."+hfnode+"."+chNode+".read()"
      torData[chNode]=eval(getNode)
  hf.close()
  return torData

def readModes(h5file,nodeName=""):
  """ given an hdf5 file, see if any mode data is in there 
      If so, return a dictionary containing the data.
      The mode data is in a group with the real and imaginary data 
      specified independently -- this returns the data as
      numpy complex arrays
      Can specify a single variable to save memory.
  """

  hf=tables.openFile(h5file,'r')
  modelist=[]
  for  hfnode in list(hf.root._v_children.keys()):
    if hfnode.endswith("_modes"):
      modelist.append(hfnode)
  if len(modelist) == 0:
    print(h5file+" has no mode data.")
    return
  modeData={}
  if nodeName:
    if modelist.count(nodeName) > 0:
      modelist=[nodeName]
    else:
      print(nodeName +" not found in file.")
      return
  for hfnode in modelist:
    childNodes="hf.root."+hfnode+"._v_children"
    real_imag=list(eval(childNodes).keys())
    if not real_imag[0].endswith("real") or not real_imag[1].endswith("imag"):
        print("Problems finding real and imaginary components of "+hfnode)
        continue
    getNode="hf.root."+hfnode+"."+real_imag[0]+".read()"
    realData=eval(getNode)
    getNode="hf.root."+hfnode+"."+real_imag[1]+".read()"
    imagData=eval(getNode)
    modeData[hfnode]=realData + 1j * imagData
  hf.close()
  return modeData

def readAllButModesAndTor(h5file,nodeName=""):
  """ given an hdf5 file that has modes in it, read everything
      but those modes.
      Return a dictionary containing the data.
      Can specify a single variable to save memory.
  """
  hf=tables.openFile(h5file,'r')
  datalist=[]
  for  hfnode in list(hf.root._v_children.keys()):
    if not hfnode.endswith("_modes") or not hfnode.endswith("_toroidal"):
      datalist.append(hfnode)
  if len(hfnode) == 0:
    print(h5file+" has no non-mode data.")
    return
  if nodeName:
    if datalist.count(nodeName) > 0:
      datalist=[nodeName]
    else:
      print(nodeName +" not found in file.")
      return
  allData={}
  for dName in datalist:
    getNode="hf.root."+dName+".read()"
    try:
      allData[dName]=eval(getNode)
    except:
      print("Problems getting data: "+dName)
  hf.close()
  return allData

def readAllButModes(h5file,nodeName=""):
  """ given an hdf5 file that has modes in it, read everything
      but those modes.
      Return a dictionary containing the data.
      Can specify a single variable to save memory.
  """
  hf=tables.openFile(h5file,'r')
  datalist=[]
  for  hfnode in list(hf.root._v_children.keys()):
    if not hfnode.endswith("_modes"):
      datalist.append(hfnode)
  if len(hfnode) == 0:
    print(h5file+" has no non-mode data.")
    return
  if nodeName:
    if datalist.count(nodeName) > 0:
      datalist=[nodeName]
    else:
      print(nodeName +" not found in file.")
      return
  allData={}
  for dName in datalist:
    getNode="hf.root."+dName+".read()"
    try:
      allData[dName]=eval(getNode)
    except:
      print("Problems getting data: "+dName)
  hf.close()
  return allData

def main():
    """ SEK: Once the pretty print can be figured out, then we can look
    at completing the standalone executable
    """
    parser = optparse.OptionParser(usage="%prog [options] inputFile")
    parser.add_option('-i', '--input', dest='input',
                      help='Name of input file if not specified as argument.',
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
