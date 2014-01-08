#USE THIS MAKEFILE TO BUILD ONETWO On LOCAL ARCHITECTURES
# FORT=pgf90,FORT=ifort  are  valid,also OPT=optimize (default) or OPT=debug
#or OPT=debug  may be set on all machines
#gmake FORT=pgf95 OPT=optimize

#nbi_dimensions_mod,nbi_com_mod.f90,nbi_types.f90,ezspline_mod.f90,
#ezcdf.f90,ezcdf_genput.f90,ezcdf_opncls.f90,....,netcdf.inc,others ?? 
# are  picked up from nubeam area . 


 

SHELL     = /bin/sh
VERSION   = v5.4
GRID      = 129_201
NAME      = onetwo_$(GRID)
HOST      = $(shell uname -s)
ONETWO    = $(shell pwd)
SHARED_MODULES = $(ONETWO)/shared_modules
VPATH = $(SHARED_MODULES)
vpath nfreya_module
INSTALLDIR = ./ 
NUBEAM_12_BUILD = FALSE
UBUILD = FALSE
INCPATH = ./
NETCDF_INC =


ifeq ($(HOST),Linux)
        PORTLIB =	portlib
	IMSL = imsl_12
	CPPLF95LIB = 
	PGFLIBS   =
	ifeq ($(shell uname -n),isis1)
           include $(INCPATH)make.isis
	endif

	ifeq ($(shell uname -n),delphi2.gat.com)
           include $(INCPATH)make.delphi2
	endif
	ifeq ($(shell uname -n),lohan1.gat.com)
           #include $(INCPATH)make.lohans
           include $(INCPATH)make.lohan1
        endif
	ifeq ($(shell uname -n),lohan3.gat.com)
           #include $(INCPATH)make.lohans
           include $(INCPATH)make.lohan3
	endif
	ifeq ($(shell uname -n),benten.gat.com)
           #include $(INCPATH)make.lohans
           include $(INCPATH)make.lohan3
	endif
        ifeq ($(shell uname -n),hestia.gat.com)
           include $(INCPATH)make.zeus
        endif
	ifeq ($(shell uname -n),lohan4.gat.com)
           include $(INCPATH)make.lohan4
	endif
	ifeq ($(shell uname -n),lohan5.gat.com)
           include $(INCPATH)make.lohan6          
	endif

	ifeq ($(shell uname -n),head)
           include $(INCPATH)make.lohan6          
	endif

        
        ifeq ($(shell uname -n),venusa)
           include $(INCPATH)make.venus
        endif

	ifeq ($(shell uname -n),node01.cluster)
           include $(INCPATH)make.lohan6          
	endif

	ifeq ($(shell uname -n),node06.cluster)
           include $(INCPATH)make.lohan6          
	endif

	ifeq ($(shell uname -n),zeus.gat.com)
           include $(INCPATH)make.zeus
	endif
	#ifeq ($(shell uname -n),hydra)
         #  include $(INCPATH)make.nfrc
	#endif
	ifeq ($(shell uname -n),jaguar1)
           include $(INCPATH)make.jaguar
	endif
	ifeq ($(shell uname -n),jaguar2)
           include $(INCPATH)make.jaguar
	endif
	ifeq ($(shell uname -n),jaguar3)
           include $(INCPATH)make.jaguar
	endif
	ifeq ($(shell uname -n),jaguar4)
           include $(INCPATH)make.jaguar
	endif
	ifeq ($(shell uname -n),jaguar5)
           include $(INCPATH)make.jaguar
	endif
	ifeq ($(shell uname -n),jaguar6)
           include $(INCPATH)make.jaguar
	endif
	ifeq ($(shell uname -n),aetna.gat.com)
           include $(INCPATH)make.aetna
	endif
	ifdef MACHINEOK
		#OPT can be changed on comand line:
		OPT=optimize
		#MPI = 1        # comment out this line if #MPI is not used
		#OPEN_MP = 1    #comment out this line if Open Mp is not used
                   #note NUBEAM modules are not compiled with openMp
                   #at this time (see compile_nubeam below,lines
                   #that start with C$ are activated with open mp directive)
		#PROFILE = 1     #comment out this line if 
				 #profileing  is not used
		#AUTO_PARA = 1   #comment out this line if auto 
				 #parallelization is not desired
		INCPATH	  = ./
	 	GCCLIB =  /usr/lib/libstdc++-3-libc6.2-2-2.10.0.so
                LISTING    = 
		ifeq ($(OPT),debug)
			OPTION = $(DEBUG) 
		else
			OPTION = $(OPTIMIZE) 
		endif
		ifndef OPEN_MP
			OMP =
		else
			#   do not include /usr/local/lib/libpgc.a for some reason or other
			#   it only allows us to use 1 cpu 
			OMP = -mp
			FTN90LIB =  \
				$(PGFLIBS)/lib/libpgftnrtl.a  \
				$(PGFLIBS)/lib/libpgthread.a \
				$(PGFLIBS)/lib/libpgmp.a
		endif

		ifndef PROFILE
			PROF = 
		else
			PROF = -Mprof=func
			PROF = -Mprof=lines
			FTN90LIB =  $(PGFLIBS)/lib/libpgc.a \
					$(PGFLIBS)/lib/libpgftnrtl.a \
					$(PGFLIBS)/lib/libpgthread.a \
					$(PGFLIBS)/lib/libpgmp.a \
					$(PGFLIBS)/lib/libpgnod_prof.a
              
		endif

		ifndef AUTO_PARA
			AUTP =
		else
			AUTP = -Mconcur -Minfo
			FTN90LIB = /usr/local/lib/libpgftnrtl.a \
				 /usr/local/lib/libpgthread.a \
				  /usr/local/lib/libpgmp.a
		endif


		LOAD      = $(FC) $(LINK)  -o $(NAME) # -m
		ARCHIVE   = ar rv
	endif

endif   # over HOSTS 


PROTECT   = chmod 755
DELETE    = rm -f
RENAME    = mv
#include $(INCPATH)/$(MAKE1) #make.includes
include $(INCPATH)/make.sources
include $(INCPATH)/$(MAKE2) #make.objects

#note to me portlib was changed !!!!!!!!!!!!!!!!
#NOTE***********************************************************************
#
#Now, ONETWO does not call the NUBEAM related libraries anymore.
#
#Note carefully that f90 modules, say a.f90, when compiled produce a.mod
#and a.o . You must link with a.o, but a.mod is also required. 
#failure to link with a.o does NOT produce a link error!!
#but the module info will not be there as expected - Linux with pgf90 
#the search order of the libraries given below does not agree
#with the order posted in nubeam readme but is necessary to get linkup!!!!!
#****************************************************************************


ifeq ($(strip $(FORT)),ifort)
##LIBRARIES = $(LIBS)/$(PORTLIB).a $(LIBS)/cortlib.a  \
##			$(LIBS)/$(IMSL).a \
##		        $(LIBS)/method_lines_library.a  $(LAPACK) \
##			$(UFILE_LIBRARIES) $(NETCDF_LIB)  $(TSPLINE_LIB)/$(SPLINE_LIB)
LIBRARIES = $(LIBS)/$(PORTLIB).a $(LIBS)/cortlib.a \
		        $(LIBS)/method_lines_library.a  $(LAPACK) \
			$(UFILE_LIBRARIES) $(NETCDF_LIB)  $(TSPLINE_LIB)/$(SPLINE_LIB)
endif
ifeq ($(strip $(FORT)),pgf95)
##LIBRARIES = $(LIBS)/$(PORTLIB).a $(LIBS)/cortlib.a   $(NETCDF_LIB) \
##			$(LIBS)/$(IMSL).a \
##		        $(LIBS)/method_lines_library.a  $(LAPACK) \
##			$(UFILE_LIBRARIES)  $(TSPLINE_LIB)/$(SPLINE_LIB)
LIBRARIES = $(LIBS)/$(PORTLIB).a $(LIBS)/cortlib.a  \
		        $(LIBS)/method_lines_library.a  $(LAPACK) \
			$(UFILE_LIBRARIES)  $(TSPLINE_LIB)/$(SPLINE_LIB)
                        # $(NETCDF_LIB) 
endif
.SUFFIXES:


all: $(NAME) data_files
	@echo
	@echo ">>  `date '+%a %d-%h-%y %r'`  `pwd`  `uname -mns`  $(LOGNAME)"
	@echo
	@echo make of $(NAME) in `pwd` on `hostname` completed
#	@echo $(UFILE_LIBRARIES)
#	@echo $(NETCDF_LIB)

$(NAME):           $(OBJECTS) $(OBJECTS_P) $(LIBRARIES)
	@echo "loading .... "
	$(LOAD)    $(OBJECTS) $(NETCDF_LIB) $(LIBRARIES)    $(FTN90LIB)
	$(PROTECT) $(NAME)
#	cp $(NAME) $(INSTALLDIR)

#JMP START

ifdef PARALLEL_12 #------------------------------------------------------------

mpi12.o: mpi12.F90 
	$(COMPILE) mpi12.F90

onetwo_parallel.o:   onetwo_parallel.F90 mpi12.o param.o
	$(COMPILE) onetwo_parallel.F90

onetwo_slave.o:    onetwo_slave.F90 mpi12.o param.o
	$(COMPILE) onetwo_slave.F90

glf_slave.o:    glf_slave.F90 mpi12.o
	$(COMPILE) glf_slave.F90

endif #-------------------------------------------------------------------------



# 
# Start Smith optimization
#
SRCS = $(shell ls *.[Ff]{,90} */*.[Ff]90)
make_depends: make_depends.sh make.objects
	@echo "Running make_depends.sh."
	@./make_depends.sh "$(OBJECTS)"
include make_depends

#
#END Smith optimization
#

##escZ.o:	escZ.c
##	gcc -c escZ.c
##
##2esc.o:	2esc.c
##	gcc -c 2esc.c

#----------------------------------------------------------------
#NUBEAM RELATED ITEMS
# note that the *.o versions of these must be in the libraries that
# are searched at link time.
#----------------------------------------------------------

#nubeam_driver.o:    nubeam_driver.f90  $(INCLUDES_nubeam_driver)
#	$(COMPILE) nubeam_driver.f90  
#
#nubeam_sub.o:    nubeam_sub.f90 
#	$(COMPILE) nubeam_sub.f90
#
#dnubeam_driver.o:    dnubeam_driver.f90  $(INCLUDES_nubeam_driver)
#	$(COMPILE) dnubeam_driver.f90  
#
#
#Nfreya_output.o:	Nfreya_output.f90 transp.o nrtype.o nub2.o nub.o \
#			P_Nfreya_12_interface.o  io.o mesh.o geom.o common_constants.mod \
#		        plasma_prop.o numbrs.o param.o
#	$(COMPILE) Nfreya_output.f90
#
#
#P_Nfreya_12_interface.o: $(INCLUDES_P_Nfreya_12_interface)
#
#P_Nfreya_nubeam_namelist.o:	P_Nfreya_nubeam_namelist.F90  nrtype.o \
#			        transp.o io_gcnmp.o solcon.o param.o nf_param.o \
#				ions.o fusion.o nub3.o	dnubeam_mod.o
#	$(COMPILE)  $< -o $@
#
#curden_terms.mod:         curden_terms.o  # used only by mmm_in.o, not used in Onetwo
#curden_terms.o:         curden_terms.F90   nrtype.o  common_constants.mod grid_class.o
#	$(COMPILE) $< -o $@ 
#
#nubeam_12.o:	nubeam_12.f $(INCLUDES_nubeam_12)
#ifeq ($(NUBEAM_12_BUILD),FALSE)
#	$(COMPILE) $(NUBEAM12)
#	mv $(basename $(NUBEAM12)).o  nubeam_12.o
#else 
#	$(COMPILE) nubeam_12.f 
#endif
#
#line_av_den.o:	line_av_den.F90 nrtype.o bicubic_spline.o plasma_prop.o \
#		common_constants.mod grid_class.o cubic_spline_class.o vector_class.o
#	$(COMPILE) $(GFLAGS) $<
#
##forcebal gcnmp interface routines:
#set_forcebal_input.o:	set_forcebal_input.F90 nrtype.o solcon_gcnmp.o \
#			common_constants.mod neutral_beams.o plasma_prop.o         \
#			grid_class.o MPI_data.o solcon_gcnmp.o
#	$(COMPILE) -DONETWO   $<
#mmm_input_defaults.o:	mmm_input_defaults.F90 nrtype.o mmm_in.o common_constants.mod \
#			modmmm7_1.o
#	$(COMPILE) $<
#
#load_input.o:         load_input.F90 iterdbmd_gcnmp.o io_gcnmp.o ions_gcnmp.o \
#		      solcon_gcnmp.o grid_class.o bc_values_gcnmp.o nrtype.o \
#	              error_handler.o MPI_data.o dep_var_gcnmp.o              \
#		      gcnmp_namelist.mod neutral_beams.o scale_profiles.o common_constants.mod \
#		      plasma_shape.o
#	$(COMPILE) -DONETWO $<
#
#get_tot_pressure.o:	get_tot_pressure.f90 nrtype.o grid_class.o ions_gcnmp.o fast_ion_data_gcnmp.o \
#			common_constants.mod dep_var_gcnmp.o neutral_beams.o error_handler.o io_gcnmp.o \
#			vector_class.o
#	$(COMPILE) $< -o $@
#
#plasma_shape.o:	plasma_shape.F90 nrtype.o plasma_prop.o common_constants.mod error_handler.o \
#		io_gcnmp.o MPI_data.o
#	$(COMPILE)  $(GFLAGS) $<
#
##setup_spawn_P_Nfreya.o:		setup_spawn_P_Nfreya.f90 nrtype.o  nfreya_namelist.o
##	$(COMPILE) $< -o $@
#
#
#ufiles_12.o:    ufiles_12.f90 $(INCLUDES_ufiles_12)
#ifeq ($(UBUILD),FALSE)
#	$(COMPILE) $(UFILES12)
#	@echo ufiles12 a  = $(UFILES12)
#	mv $(basename $(UFILES12)).o  ufiles_12.o
#else
#ifeq ($(FORT),lf95)
#	$(COMPILELF95) ufiles_12.f90 
#else
#	$(COMPILE) ufiles_12.f90 
#	@echo ufiles12 b = $(UFILES12)
#endif
#endif
#
#nbnamelist.mod:	 nbnamelist.o 
#nbnamelist.o:    nbnamelist.F90  nrtype.o param.o transp.o numbrs.o
#	$(COMPILE)  -DONETWO  $< -o $@  
#
#nbi_restart.mod:	 nbi_restart.o
#nbi_restart.o:    nbi_restart.f90 $(INCLUDES_nbi_restart)
#ifeq ($(MACHINE),taurus)
#	$(COMPILENODB) nbi_restart.f90
#else
#	$(COMPILE) nbi_restart.f90  
#	#mv $(basename $(nbi_restart)).o  nbi_restart.o
#endif
#
#get_nboutput.o:    get_nboutput.f90 $(INCLUDES_get_nboutput)
#	$(COMPILE) get_nboutput.f90 
#
#transp.mod:  transp.o 
#transp.o:  transp.f90  beam_structure.o fast_ion_diffusion.o
#	$(COMPILE_nubeam) transp.f90
#
#
#ifeq ($(HOST),HP-UX)
#	cp EVENTS.mod events.mod
#endif
#
##modules from transp(NOT from nubeam)
#xpload_module.mod:  xpload_module.o
#xpload_module.o:  xpload_module.f90  
#	$(COMPILE_nubeam) xpload_module.f90
#
##end nubeam modules --------------------------------------------------------




$(LIBS)/method_lines_library.a: $(LIBS)/method_lines_library.o
	$(DELETE)  $(LIBS)/method_lines_library.a
	$(ARCHIVE) $(LIBS)/method_lines_library.a   $(LIBS)/method_lines_library.o
$(LIBS)/method_lines_library.o: $(LIBS)/method_lines_library.f
	cd         $(LIBS) ; \
	$(COMPILE) $(LIBS)/method_lines_library.f




$(LIBS)/cortlib.a: $(LIBS)/cortlib.o
	$(DELETE)  $(LIBS)/cortlib.a
	$(ARCHIVE) $(LIBS)/cortlib.a   $(LIBS)/cortlib.o
$(LIBS)/cortlib.o: $(LIBS)/cortlib.c
	echo $(PGCC)
	cd         $(LIBS) ; \
	$(PGCC)  -c    $(LIBS)/cortlib.c




$(LIBS)/$(PORTLIB).a: $(LIBS)/$(PORTLIB).o
	$(DELETE)  $(LIBS)/$(PORTLIB).a
	$(ARCHIVE) $(LIBS)/$(PORTLIB).a   $(LIBS)/$(PORTLIB).o
$(LIBS)/$(PORTLIB).o: $(LIBS)/$(PORTLIB).F
	cd         $(LIBS) ; \
	$(COMPILE) -I../ $(LIBS)/$(PORTLIB).F



##$(LIBS)/$(IMSL).a:    $(LIBS)/$(IMSL).o
##	$(DELETE)  $(LIBS)/$(IMSL).a
##	$(ARCHIVE) $(LIBS)/$(IMSL).a      $(LIBS)/$(IMSL).o
##$(LIBS)/$(IMSL).o:    $(LIBS)/$(IMSL).f
##	cd         $(LIBS) ; \
##	$(COMPILE) $(LIBS)/$(IMSL).f


.PHONY:	clean
clean:
	@echo cleaning up files on $(MACHINE)
	 $(DELETE) $(NAME) *.a *.o *.mod
	 $(DELETE) $(LIBS)/*.o $(LIBS)/*.a 
	 $(DELETE) $(SHARED_MODULES)/*.o
	 #@echo $(LIBS)

data_files: 
	cd $(ONETWO)

test_piksr2: test_piksr2.o 
	$(LOAD) $^ -o test_piksr2
