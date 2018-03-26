#----------------------------------------------------------------------
# gacodeinput.py
#
# PURPOSE:
#  Collection of classes for parsing of GACODE free-format input files.
#
# NOTES: 
#  SimpleInput  : input.gyro, input.neo, input.tglf input.glf23
#  ProfileInput : input.profiles
#  ManagerInput : input.tgyro [see tgyro/bin/tgyro_parse.py]
#----------------------------------------------------------------------

import string
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

    def add(self,param,default):
        self.data_dict[param]=default
        self.data_orderlist.append(param)

    def dep(self,param,default):
        self.dep_dict[param]=default
        self.dep_orderlist.append(param)

    def printdebug(self):
        print self.data_dict
        print self.data_orderlist
        print self.dep_dict
        print self.dep_orderlist
        print self.user_dict
        if self.error == 1:
            print self.error_msg

    def printmsg(self):
        if self.error == 1:
            print self.error_msg
        
    def set_extension(self,text):
        self.extension = text

    def read_input(self,inputfile):
        # 1. read user input file
        for line in open(inputfile,'r').readlines():

            # Remove leading and trailing whitespace from line
            line = string.strip(line)

            # Skip blank lines
            if len(line) > 0 and line[0] != '#':
                x = string.splitfields(line,'=')
                y = string.splitfields(x[1],'#')
                #y = string.splitfields(y[0],' ')
                arg = string.strip(x[0])
                val = string.strip(y[0])

                self.user_dict[arg] = val

        # 2. build complete input file, looking for errors
        for x in self.user_dict.keys():
            if self.data_dict.has_key(x) == 1:
                self.data_dict[x] = self.user_dict[x]
            elif self.dep_dict.has_key(x) == 1:
                self.error=1
                self.error_msg=self.error_msg+'ERROR: (gacodeinput) Deprecated parameter '+x+'\n'
                self.error_msg=self.error_msg+'       '+self.dep_dict[x]+'\n'
            else:
                self.error=1
                self.error_msg=self.error_msg+"ERROR: (gacodeinput) Bogus parameter "+x+'\n'

        if self.error == 0:
            f=open(inputfile+self.extension,'w')
            for x in self.data_orderlist:
                f.write(self.data_dict[x]+'  '+x+'\n')
        

#--------------------------------------------------------------------
# PARSER FOR input.profiles
#--------------------------------------------------------------------
class ProfileInput:
    """Input parser for input.profiles"""
    def __init__(self):
        self.data_dict = {}
        self.data_orderlist = []
        self.user_dict = {}
        self.error = 0
        self.error_msg = ""
        self.extension = ".gen"

    def add(self,param,default):
        self.data_dict[param]=default
        self.data_orderlist.append(param)

    def printmsg(self):
        if self.error == 1:
            print self.error_msg

    def set_extension(self,text):
        self.extension = text

    def read_input(self,inputfile):
        # 1. read user input file
        for line in open(inputfile,'r').readlines():

            # Remove leading and trailing whitespace from line
            line = string.strip(line)
                
            # Skip blank lines
            if len(line) > 0 and line[0] != '#':
                x = string.splitfields(line,'=')
                y = string.splitfields(x[1],'#')
                arg = string.strip(x[0])
                val = string.strip(y[0])

                self.user_dict[arg] = val

        # 2. build complete input file, looking for errors
        for x in self.user_dict.keys():
            if self.data_dict.has_key(x) == 1:
                self.data_dict[x] = self.user_dict[x]
            else:
                self.error=1
                self.error_msg=self.error_msg+"ERROR: Bogus parameter "+x+'\n'


    def read_profile(self,inputfile):

        outputfile = inputfile+self.extension

        file_out = open(outputfile,'w')
        file_temp = open('parse_temp','w')

        profile_data = []

        for line in open(inputfile,'r').readlines():

            line = string.strip(line)

            # Split inputfile (input.profiles) into scalar and 
            # vector data:
            if 'SHOT' in line:
                x = string.splitfields(line,':')[1]
                x = string.strip(x) 
                if len(x) == 0:
                    x = '0'
                file_temp.write('SHOT='+x+'\n')

            if (len(line) > 0) and (line[0] != '#'):                
                if string.find(line,'=') > -1:
                    # Write scalar data into temp file
                    file_temp.write(line+'\n')
                else:
                    # Originally '  ' was the pattern:       
                    #  data = string.splitfields(line,'   ')
                    data = string.split(line)
                    ncol = len(data)
                    for j in range (0,ncol):
                        # Save vector data in variable v.
                        profile_data.append(data[j])

        file_temp.close()

        # Parse the scalar data in file_temp
        self.read_input('parse_temp')

        # Compute number of rows for profile data
        ncol = 5
        nrow = int(self.data_dict['N_EXP'])
        nblock = len(profile_data)/(nrow*ncol)

        for x in self.data_orderlist:
            file_out.write(self.data_dict[x]+'  '+x+'\n')

        # Write vector data
        for k in range(nblock):
            for j in range(ncol):
                for i in range(nrow):
                    indx = ncol*i+j+k*nrow*ncol
                    file_out.write(profile_data[indx]+'\n')
    
        # Clean up temporary file
        os.system('rm parse_temp')


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
            print self.error_msg

    def write_proc(self,datafile):
        f = open(datafile,'w')
        f.write(str(len(self.slavepath))+'\n')
        f.write(str(self.sum_proc))
        f.close()

    def set_extension(self,text):
        self.extension = text

    def strip_tag(self,datafile):

        file_input = open(datafile+'.input','w')

        n = 0
        for line in open(datafile,'r').readlines():
            line_s = string.strip(line)

            # Look for occurence of tag and put item in list.
            if (line_s[0:3] == 'DIR'):   
                n = n+1
                data = string.splitfields(line_s,' ')

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
                        self.slaveradius.append(string.splitfields(data[3],'=')[1])
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
                file_overlay = open('overlay.'+str(n),'w')

                #----------------------------------------------------------
                # This loop writes each overlay parameter list to overlay.*
                for j in range(nover):
                    file_overlay.write(data[j+nj]+'\n')
                file_overlay.close()
                #----------------------------------------------------------
 
            else:

                file_input.write(line)

        file_input.close()


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
        file_outfile = open(outfile,'a')
        file_outfile.write(str(n_path)+'\n')

        # Logging
        print 'INFO: (gacodeinput) Number of code instances: '+str(n_path)

        for p in range(len(self.slavepath)):
            self.sum_proc = self.sum_proc + int(self.slaveproc[p])
            basedir = self.slavepath[p]

            # Detect the code to be run in each directory            
            if os.path.isfile(basedir+'/input.gyro'):
                code='gyro'
            elif os.path.isfile(basedir+'/input.cgyro'):
                code='cgyro'
            elif os.path.isfile(basedir+'/input.tglf'):
                code='tglf'
            elif os.path.isfile(basedir+'/input.ifs'):
                code='ifs'
            elif os.path.isfile(basedir+'/input.glf23'):
                code='glf23'
            else:
                code='unknown'
                self.error=1
                self.error_message='Could not identify code'

            file_outfile.write(basedir+' '+self.slaveproc[p]+' '+self.slaveradius[p]+' '+code+'\n') 

            if code == 'unknown':
                print 'ERROR: (gacodeinput.py) No code found in '+basedir
                continue
            if code == 'ifs':
                print 'INFO: (gacodeinput.py) Found ifs input in '+basedir
                continue
            else:
                print 'INFO: (gacodeinput.py) Found '+code+' input in '+basedir

            if os.path.isfile(basedir+'/input.profiles'):
               os.system('python $GACODE_ROOT/shared/bin/profile_parse.py '+basedir+'/input.profiles')
 
            basefile = basedir+'/input.'+code 
            tempfile = basefile+'.temp' 

            file_base = open(basefile,'r')
            file_temp = open(tempfile,'w')

            for line in file_base.readlines():
                if line[0:18] <> "# -- Begin overlay":
                    file_temp.write(line)
                else:
                    break

            file_base.close()
            file_temp.close()

            # Overlay parameters
            os.system('echo "# -- Begin overlay" >> '+tempfile)
            os.system('cat '+self.overlayfile[p]+' >> '+tempfile)
            os.system('mv '+tempfile+' '+basefile)

            # Run code in test mode
            os.system(code+' -i '+basedir+' -n 1 -p $PWD > out.log')
                  
            os.system('rm '+self.overlayfile[p])

        print 'INFO: (gacodeinput) Required MPI tasks in TGYRO: '+str(self.sum_proc)

