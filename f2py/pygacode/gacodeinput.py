#----------------------------------------------------------------------
# gacodeinput.py
#
# PURPOSE:
#  Collection of classes for parsing of GACODE free-format input files.
#
# NOTES:
#  SimpleInput  : input.cgyro, input.neo, input.tglf
#  ManagerInput : input.tgyro [see tgyro/bin/tgyro_parse.py]
#----------------------------------------------------------------------

import os

#--------------------------------------------------------------------
# PARSER FOR input.cgyro, etc.
#--------------------------------------------------------------------
class SimpleInput:
    """Input parser for simple input.* files"""
    def __init__(self):
        self.data_dict = {}
        self.data_orderlist = []
        self.dep_dict = {}
        self.dep_orderlist = []
        self.user_dict = {}
        self.error = 0
        self.error_msg = ""
        self.extension = ".gen"

    def add(self,param,default,n=1):
        if n == 1:
            self.data_dict[param]=default
            self.data_orderlist.append(param)
        else:
            for i in range(n):
                p = param+'_'+str(i+1)
                self.data_dict[p]=default
                self.data_orderlist.append(p)
              
    def dep(self,param,default):
        self.dep_dict[param]=default
        self.dep_orderlist.append(param)

    def printdebug(self):
        print(self.data_dict)
        print(self.data_orderlist)
        print(self.dep_dict)
        print(self.dep_orderlist)
        print(self.user_dict)
        if self.error == 1:
            print(self.error_msg)

    def printmsg(self):
        if self.error == 1:
            print(self.error_msg)

    def set_extension(self,text):
        self.extension = text

    def read_input(self,inputfile):
        # 1. read user input file
        with open(inputfile,'r') as fin:
           for line in fin.readlines():

               # Remove leading and trailing whitespace from line
               line = line.strip()

               # Skip blank lines
               if len(line) > 0 and line[0] != '#':
                   x = line.split('=')
                   y = x[1].split('#')
                   arg = x[0].strip()
                   val = y[0].strip()

                   self.user_dict[arg] = val

        # 2. build complete input file, looking for errors
        for x in list(self.user_dict.keys()):
            if x in self.data_dict:
                self.data_dict[x] = self.user_dict[x]
            elif x in self.dep_dict:
                 if self.dep_dict[x] == 'ignore':
                     print('WARNING: (gacodeinput.py) Ignoring parameter '+x)
                 else:
                     self.error=1
                     self.error_msg=self.error_msg+'ERROR: (gacodeinput) Deprecated parameter '+x+'\n'
                     self.error_msg=self.error_msg+'       '+self.dep_dict[x]+'\n'
            else: 
                self.error=1
                self.error_msg=self.error_msg+"ERROR: (gacodeinput) Bogus parameter "+x+'\n'

        if self.error == 0:
            with open(inputfile+self.extension,'w') as f:
               for x in self.data_orderlist:
                   f.write(self.data_dict[x]+'  '+x+'\n')
        
#--------------------------------------------------------------------
# PARSER FOR input.tgyro
#--------------------------------------------------------------------
class ManagerInput:
    """Input parser for input.tgyro"""
    def __init__(self):
        self.data_dict = {}
        self.data_orderlist = []
        self.dep_dict = {}
        self.dep_orderlist = []
        self.user_dict = {}
        self.error = 0
        self.error_msg = ""
        self.extension = ".gen"
        self.slavepath = []
        self.slaveproc = []
        self.slaveradius = []
        self.overlayfile = []
        self.sum_proc = 0

    def add(self,param,default):
        self.data_dict[param]=default
        self.data_orderlist.append(param)

    def dep(self,param,default):
        self.dep_dict[param]=default
        self.dep_orderlist.append(param)

    def printmsg(self):
        if self.error == 1:
            print(self.error_msg)

    def write_proc(self,datafile):
        with open(datafile,'w') as f:
           f.write(str(len(self.slavepath))+'\n')
           f.write(str(self.sum_proc))

    def set_extension(self,text):
        self.extension = text

    def strip_tag(self,datafile):

        with open(datafile+'.input','w') as file_input:

           n = 0
           with open(datafile,'r') as fdata:
              datalines = fdata.readlines()

           adata = []
           for line in datalines:
               line_s = line.strip()
               
               # Look for occurence of tag and put item in list.
               if (line_s[0:3] == 'ALL'):
                   # ALL overlay is applied to all DIR entries 
                   adata = line_s.split(' ')
               elif (line_s[0:3] == 'DIR'):
                   n = n+1
                   data = line_s.split(' ')

                   # data[0] -> DIR
                   # data[1] -> directory1, etc
                   # data[2] -> <n_cores>
                   # data[3,...] -> OVERLAY_VARIABLE [special option X=<xmin>]

                   # slavepath stores directory
                   self.slavepath.append(data[1])
                   # slaveproc stores number of cores
                   self.slaveproc.append(data[2])

                   # Now manage overlays and optional radii
                   if (len(data) > 3):
                       # Overlay or optional radius
                       if data[3][0:1] == 'X':
                           # This is the special option X=<xmin> for min(r/a) or min(rho)
                           self.slaveradius.append(data[3].split('=')[1])
                           # Need to subtract 4 because X is not an overlay
                           nover = len(data)-4
                           nj    = 4
                       else:
                           nover = len(data)-3
                           nj    = 3
                           self.slaveradius.append("-1")
                   else:
                       # No overlay
                       self.slaveradius.append("-1")
                       nover = 0

                   # Overlay parameters reside in data[3], ... 
                   self.overlayfile.append('overlay.'+str(n))
                   with open('overlay.'+str(n),'w') as file_overlay:

                      for j in range(len(adata)-1):
                          file_overlay.write(adata[j+1]+'\n')
                      #----------------------------------------------------------
                      # This loop writes each overlay parameter list to overlay.*
                      for j in range(nover):
                          file_overlay.write(data[j+nj]+'\n')

                   #----------------------------------------------------------

               else:

                   file_input.write(line)


    def read_input(self,datafile):

        # Eventual output file
        # For example: datafile = input.cgyro
        #             extension = .gen
        outfile = datafile+self.extension

        # Split datafile into datafile.input (pure input) 
        # and overlay.n, where n is the number of tags.
        self.strip_tag(datafile)

        # Parse simple input part
        z = SimpleInput()
        z.data_dict = self.data_dict
        z.data_orderlist = self.data_orderlist
        z.dep_dict = self.dep_dict

        # NOTE: 'datafile'.input (input.tgyro.input) contains the TGYRO inputs 
        # (with DIR lines stripped)
        z.read_input(datafile+'.input')
        self.error = z.error
        self.error_msg = z.error_msg

        # Something went wrong in SimpleInput (bogus parameter),
        # so let's just return
        if self.error == 1:
            return
        
        # File generated above
        genfile = datafile+'.input'+z.extension 

        n_path = len(self.slavepath)

        os.system('mv '+genfile+' '+outfile)
        os.system('rm '+datafile+'.input')
        with open(outfile,'a') as file_outfile:
           file_outfile.write(str(n_path)+'\n')

           # Logging
           print('INFO: (gacodeinput) Number of code instances: '+str(n_path))

           for p in range(len(self.slavepath)):
               self.sum_proc = self.sum_proc + int(self.slaveproc[p])
               basedir = self.slavepath[p]

               # Detect the code to be run in each directory            
               if os.path.isfile(basedir+'/input.gyro'):
                   code='gyro'
               elif os.path.isfile(basedir+'/input.qlgyro'):
                   code='qlgyro'
               elif os.path.isfile(basedir+'/input.cgyro'):
                   code='cgyro'
               elif os.path.isfile(basedir+'/input.tglf'):
                   code='tglf'
               elif os.path.isfile(basedir+'/input.ifs'):
                   code='ifs'
               elif os.path.isfile(basedir+'/input.glf23'):
                   code='glf23'
               elif os.path.isfile(basedir+'/input.etg'):
                   code='etg'
               elif os.path.isfile(basedir+'/input.mmm'):
                   code='mmm'
               else:
                   code='unknown'
                   self.error=1
                   self.error_message='Could not identify code'

               file_outfile.write(basedir+' '+self.slaveproc[p]+' '+self.slaveradius[p]+' '+code+'\n') 

               if code == 'unknown':
                   print('ERROR: (gacodeinput.py) No code found in '+basedir)
                   continue
               if code == 'ifs':
                   print('INFO: (gacodeinput.py) Found ifs input in '+basedir)
                   continue
               elif code == 'mmm':
                   print('INFO: (gacodeinput.py) Found mmm input in '+basedir)
                   continue
               else:
                   print('INFO: (gacodeinput.py) Found '+code+' input in '+basedir)

               basefile = basedir+'/input.'+code
               tempfile = basefile+'.temp'

               with open(basefile,'r') as file_base:
                   with open(tempfile,'w') as file_temp:

                      for line in file_base.readlines():
                          if line[0:18] != "# -- Begin overlay":
                              file_temp.write(line)
                          else:
                              break

               # Overlay parameters
               os.system('echo "# -- Begin overlay" >> '+tempfile)
               os.system('cat '+self.overlayfile[p]+' >> '+tempfile)
               os.system('mv '+tempfile+' '+basefile)

               # Run code in test mode
               os.system(code+' -i '+basedir+' -n 1 -p $PWD > out.log')

               os.system('rm '+self.overlayfile[p])

        print('INFO: (gacodeinput) Required MPI tasks in TGYRO: '+str(self.sum_proc))

