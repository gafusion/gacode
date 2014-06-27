#USE THIS MAKEFILE TO BUILD stand alone P_Nfreya module On LOCAL ARCHITECTURES
# FORT=pgf90,FORT=ifort  are  valid,also OPT=optimize (default) or OPT=debug
#may be set on all machines 
# example: gmake -f Make_PNfreya FORT=pgf95 OPT=optimize
#cpp preprocessor is used to preprocess *.F90 files
#NOTE THAT ONETWO LOOKS FOR P_NFREYA_PGF_L1 OR P_NFREYA_IFORT_L6(IN THIS DIRECTORY)
# SO CHANGE THE NAME AFTER COMPILATION


SHELL     = /bin/sh
VERSION   = v0.9
NAME      = P_Nfreya
HOST = $(shell uname -s)
INSTALLDIR = ./ 
FFLAGS =
MPI_LIB =
NETCDF_LIB = 
#USEMPI =1
MPI2_EXT = 1


ifeq ($(HOST),Linux)
	ifeq ($(shell uname -n),lohan1.gat.com)
		include ./make.lohan1
		ONETWO          =  /p/linux/onetwo/source/onetwo
		SHARED_MODULES  = $(ONETWO)/shared_modules
		VPATH           = $(SHARED_MODULES)
                ONETWO          = /usc-data/p2/linux/onetwo/source/onetwo
                IMSL            = imsl_12
                LIBS            = $(ONETWO)/onetwo_libs/$(FORT)/
		LIBSM           = $(ONETWO)/onetwo_libs
               # MISC_LIBS       = $(LIBS)/imsl_12.a
               OPTIONS = -r8
	endif  #lohan1
	ifeq ($(shell uname -n),lohan4.gat.com)
		include ./make.lohan4
		ONETWO          =  /p/linux/onetwo/source/onetwo
		SHARED_MODULES  = $(ONETWO)/shared_modules
		VPATH           = $(SHARED_MODULES)
		LIBSM           = $(ONETWO)/onetwo_libs
               OPTIONS = -r8
	endif  #lohan4

	ifeq ($(shell uname -n),lohan5.gat.com)
		include ./make.lohan5
		ONETWO          =  /p/linux/onetwo/source/onetwo
		SHARED_MODULES  = $(ONETWO)/shared_modules
		VPATH           = $(SHARED_MODULES)
		LIBSM           = $(ONETWO)/onetwo_libs
               OPTIONS = -r8
	endif  #lohan5

	ifeq ($(shell uname -n),head)
		include ./make.lohan6
		ONETWO          =  /p/linux/onetwo/source/onetwo
		SHARED_MODULES  = $(ONETWO)/shared_modules
		VPATH           = $(SHARED_MODULES)
                ONETWO          = /usc-data/p2/linux/onetwo/source/onetwo
                IMSL            = imsl_12
                LIBS            = $(ONETWO)/onetwo_libs/$(FORT)/
              #  MISC_LIBS       = $(LIBS)/imsl_12.a
		LIBSM           = $(ONETWO)/onetwo_libs
               OPTIONS = -r8
	endif  #lohan6
	ifeq ($(shell uname -n),node06.cluster)
		include ./make.lohan6
		ONETWO          =  /p/linux/onetwo/source/onetwo
		SHARED_MODULES  = $(ONETWO)/shared_modules
		VPATH           = $(SHARED_MODULES)
                ONETWO          = /usc-data/p2/linux/onetwo/source/onetwo
                IMSL            = imsl_12
                LIBS            = $(ONETWO)/onetwo_libs/$(FORT)/
            #    MISC_LIBS       = $(LIBS)/imsl_12.a
		LIBSM           = $(ONETWO)/onetwo_libs
               OPTIONS = -r8
	endif  #node06

	ifeq ($(shell uname -n),venusa)
		include ./make.venusa
		ONETWO          =  /p/linux/onetwo/source/onetwo
    ONETWO = ..
		SHARED_MODULES  = $(ONETWO)/shared_modules
		VPATH           = $(SHARED_MODULES):$(ONETWO)
                ONETWO          = /usc-data/p2/linux/onetwo/source/onetwo
                ONETWO = ..
                IMSL            = imsl_12
                LIBS            = $(ONETWO)/onetwo_libs/$(FORT)/
            #    MISC_LIBS       = $(LIBS)/imsl_12.a
		LIBSM           = $(ONETWO)/onetwo_libs
               OPTIONS = -r8
	endif  

	ifeq ($(shell uname -n),vidar)
		include ./make.vidar
		ONETWO          =  /p/linux/onetwo/source/onetwo
    ONETWO = ..
		SHARED_MODULES  = $(ONETWO)/shared_modules
		VPATH           = $(SHARED_MODULES):$(ONETWO)
                ONETWO          = /usc-data/p2/linux/onetwo/source/onetwo
                ONETWO = ..
                IMSL            = imsl_12
                LIBS            = $(ONETWO)/onetwo_libs/$(FORT)/
            #    MISC_LIBS       = $(LIBS)/imsl_12.a
		LIBSM           = $(ONETWO)/onetwo_libs
		OPTIONS = -fdefault-real-8
	endif  

	ifeq ($(shell uname -n),star7.gat.com)
		include ./make.star7
		ONETWO          =  /p/linux/onetwo/source/onetwo
		SHARED_MODULES  = $(ONETWO)/shared_modules
		VPATH           = $(SHARED_MODULES)
                ONETWO          = /usc-data/p2/linux/onetwo/source/onetwo
            #    IMSL            = imsl_12
                LIBS            = $(ONETWO)/onetwo_libs/$(FORT)/
            #    MISC_LIBS       = $(LIBS)/imsl_12.a
		LIBSM           = $(ONETWO)/onetwo_libs
               OPTIONS = -r8
	endif  

	ifdef MACHINEOK
		#OPT can be changed on comand line:
		OPT=optimize
		INCPATH	  = ./
                LISTING    = 
		ifeq ($(OPT),debug)
			OPTION = $(DEBUG) 
		else
			OPTION = $(OPTIMIZE)
		endif
                OMP =
	 	PROF = 
		AUTP =

		LOAD      = $(FC) $(LINK) $(PROF)   -o $(NAME) # -m
		ARCHIVE   = ar rv
	endif   # machineok


endif   # over HOSTS 

GFLAGS    = $(FFLAGS) -DNFREYA    # defines NFREYA  for use in conditional compilation
PROTECT   = chmod 755
DELETE    = rm -f
RENAME    = mv

#Nfreya specific modules:
SOURCES_Nf    = nfreya_driver.F90                  \
                nfreya_version.f90                 \
		load_nfreya_input.F90              \
		distribute_run_directives.F90      \
		nfreya_namelist_defaults.f90       \
	        prep_zonal_input.f90               \
		nfreya_init.F90                    \
		contour_module.f90                 \
		zonal_data.f90                     \
	        set_ion_prop2.f90                  \
		2d_interp.f90                      \
                tport_mhd_grid_data.f90            \
		xsct_class.F90                     \
                cpub_dat.f90                       \
                hexnb_data.f90                     \
		d01bff_dat.f90                     \
		nfreya_routines.f90                  \
		1d_spline.f90                      \
		nfreya_run.F90                     \
		ran_num_gen.f90                    \
		copy_injector_attributes.f90       \
		nfreya_write_statefile.f90         \
		prep_transport_output.F90          \
		nblist_read.f90                    \
		nfreya_load.f90                    \
		distribute_ufile.F90               \
#                distribute_nubeam_namelist.f90




# shared modules:
SOURCES_Sh   =  nrtype.F90			   \
		getargs.f90                        \
		vector_class.f90                   \
		consts.f90                         \
		MPI_data.F90                       \
		error_handler.F90                  \
		io_gcnmp.f90                       \
		ions_gcnmp.f90                     \
		solcon_gcnmp.f90                   \
		glf23_gcnmp.f90                    \
		dep_var_gcnmp.f90                  \
		iterdbmd_gcnmp.f90                 \
		bc_values_gcnmp.F90                \
		curden_terms.F90                   \
		scale_profiles.f90                 \
		utils.F90                          \
		rw_iterdb_netcdf.F90               \
		rw_iterdb.F90                      \
		grid_class.F90                     \
		fast_ion_data_gcnmp.f90            \
		fast_ion_diffusion.f90             \
		fdyimp.f90 gcnmp_version.f90       \
		tglfin.f90 tglf_modules.f90        \
		distribute_statefile.F90           \
		tspline.f90                        \
		strip_comments.F90                 \
		neutral_beams.F90                  \
		plasma_prop.f90                    \
		source_terms_gcnmp.f90             \
		shot_info.f90                      \
		neutral_data.f90                   \
		nfreya_namelist.F90                \
		beam_structure.F90                 \
		allocate_err.f90                   \
		param.f90                          \
		numbrs.f90                         \
		string_util.F90                    \
		periodic_table_mod.f90             \
		nubeam_sub.f90                     \
                nf_param.f90                       \
		P_Nfreya_interface.F90             \
		fusion_gcnmp.f90                   \
		tcoef.f90                          \
		nub.F90                            \
		get_tot_pressure.f90               \
                dnubeam_mod.f90                    \
                replace_imsl.f90                   \
		m_refsor.f90                       \
		bicube.f90                         \
		cubic_spline_class.f90             \
		mhdpar.f90                         \
		bicubic_spline.f90                 \
		tglfv1.94_modules.f90              \
		allocate_plasma_freq.f90
##		$(LIBS)/imsl_12.f  
#		nbnamelist.F90                      no longer used


SOURCES   = $(SOURCES_Nf) $(SOURCES_Sh)



OBJECT1=$(SOURCES:.f90=.o)
OBJECT2=$(OBJECT1:.f=.o)
OBJECTS=$(OBJECT2:.F90=.o)

#-------------------------------------------------------------------------------
# the following magic was very hard to come by and took me 2 full days to resolve HSJ:
#UFILE_LIBRARIES  = $(UFILE_LIBS)/libuflib.a         \
#                   $(UFILE_LIBS)/libufhdf.a        \
#                   $(UFILE_LIBS)/libmdstransp.a    \
#                   $(UFILE_LIBS)/libvaxonly.a      \
#                   $(UFILE_LIBS)/libportlib.a      \
#                   $(UFILE_LIBS)/libmds_sub.a      \
#		   /f/mdsplus/linux64/mdsplus/lib/libMdsLib_client.a \
#		   /f/mdsplus/linux64/mdsplus/lib/libMdsIpShr.a \
#		   /f/mdsplus/linux64/mdsplus/lib/libTdiShr.a \
#		   $(UFILE_LIBS)/libhdf_dummy.a  

#/f/mdsplus/linux64/mdsplus/lib libraries *.a  can be repalced with *.so versions
#-------------------------------------------------------------------------------

LIBRARIES       =	$(LIBSM)/method_lines_library.a \
			$(LAPACK) $(MPELIB) $(MPI_LIB)  $(SPLINELIB) \
			  $(MISC_LIBS)   $(UFILE_LIBRARIES) #$(NETCDF_LIB) \



all: $(NAME)
	@echo 
	@echo ">>  `date '+%a %d-%h-%y %r'`  `pwd`  `uname -mns`  $(LOGNAME)"
	@echo
	@echo make of $(NAME) in `pwd` on `hostname` completed

$(NAME):           $(OBJECTS)  $(LIBRARIES)
	@echo in nme uflibs = $(UFILE_LIBRARIES)
	@echo in nme  = $(UFILE_LIBS)
	@echo in nme NETCD_LIB = $(NETCDF_LIB)
	$(LOAD)   $(GFLAGS) $(OMP)   $(OBJECTS)  $(LIBRARIES)  $(NETCDF_LIB)
	$(PROTECT) $(NAME)
	#echo uflibs = $(UFILE_LIBRARIES)

test_random: test_random.o
	$(FC) $< -o test_random
test_random.o: test_random.f90 
	$(COMPILE) $<

allocate_err.o:	allocate_err.f90 nrtype.o
	$(COMPILE) $<

allocate_plasma_freq.o:	allocate_plasma_freq.f90  $(INCLUDES_pl_freq)
	$(COMPILE)  -DONETWO -I$(NETCDF_INC)  $< -o  $@ 

INCLUDES_bc_values_gcnmp = nrtype.o error_handler.o grid_class.o solcon_gcnmp.o \
		     ions_gcnmp.o plasma_prop.o vector_class.o        \
		     io_gcnmp.o curden_terms.o curden_terms.o
bc_values_gcnmp.o:         bc_values_gcnmp.F90   $(INCLUDES_bc_values_gcnmp)
	$(COMPILE)  $(GFLAGS)  $< 

beam_structure.o:	beam_structure.F90  nrtype.o allocate_err.o
	$(COMPILE)  $(GFLAGS)  $<  
        


consts.o:       consts.f90   nrtype.o 
	$(COMPILE) $< 

copy_injector_attributes.o :	copy_injector_attributes.f90  nrtype.o
	$(COMPILE) $< 	

cpub_dat.o:       cpub_dat.f90 consts.o   nrtype.o 
	$(COMPILE) $< 


INCLUDES_curden_terms =		nrtype.o  consts.o dep_var_gcnmp.o grid_class.o \
				source_terms_gcnmp.o
curden_terms.o:         curden_terms.F90 $(INCLUDES_curden_terms)
	$(COMPILE) $(GFLAGS) $< 


d01bff_dat.o:	d01bff_dat.f90
	$(COMPILE)  $< 


fusion_gcnmp.o:	fusion_gcnmp.f90 nrtype.o consts.o grid_class.o ions_gcnmp.o \
		dep_var_gcnmp.o
	$(COMPILE)  $< 

hexnb_data.o:	hexnb_data.f90 nrtype.o
	$(COMPILE)  $< 

#nbnamelist.o:	nbnamelist.F90 nrtype.o param.o  numbrs.o P_Nfreya_interface.o
#	$(COMPILE)  $<  not used HSJ 2/25/13

nrtype.o:         nrtype.F90   
	@echo  FORT=$(FORT)	
	@echo COMPILE   = $(COMPILE)
	@echo USEMPI    = $(USEMPI)
	@echo MPI2_EXT  = $(MPI2_EXT)
	@echo MPI_INC   = $(MPI_INC)
	@echo MSG = $(MSG)
	$(COMPILE) $< -o $@ 

nub.o:	nub.F90  nrtype.o
	$(COMPILE)  $(GFLAGS)  $< -o $@

nfreya_namelist.o:	nfreya_namelist.F90 neutral_beams.o nrtype.o io_gcnmp.o \
			grid_class.o curden_terms.o error_handler.o zonal_data.o \
			xsct_class.o nf_param.o P_Nfreya_interface.o nub.o
	$(COMPILE) $< -o $@

nfreya_namelist_defaults.o:	nfreya_namelist_defaults.f90 nrtype.o bc_values_gcnmp.o  \
				consts.o   source_terms_gcnmp.o nf_param.o zonal_data.o     \
				curden_terms.o error_handler.o tglfin.o neutral_beams.o  \
				xsct_class.o P_Nfreya_interface.o nub.o
	$(COMPILE) $< -o $@ 


nfreya_init.o:	nfreya_init.F90  xsct_class.o prep_zonal_input.o consts.o \
		neutral_beams.o tport_mhd_grid_data.o nfreya_routines.o     \
		error_handler.o MPI_data.o plasma_prop.o shot_info.o
	$(COMPILE) $<

nfreya_load.o:	nfreya_load.f90 nrtype.o P_Nfreya_interface.o  \
		 neutral_beams.o error_handler.o io_gcnmp.o    \
		 zonal_data.o periodic_table_mod.o consts.o shot_info.o
	$(COMPILE) $<


nfreya_run.o:	nfreya_run.F90 neutral_beams.o ions_gcnmp.o nrtype.o       \
		nfreya_routines.o tport_mhd_grid_data.o consts.o plasma_prop.o
	$(COMPILE)    $(GFLAGS)   $<	nfreya_run.F90

nfreya_version.o:	nfreya_version.f90
	$(COMPILE) nfreya_version.f90

nfreya_routines.o:	nfreya_routines.f90 nrtype.o neutral_beams.o zonal_data.o\
		        xsct_class.o 1d_spline.o plasma_prop.o ran_num_gen.o   \
			tport_mhd_grid_data.o grid_class.o P_Nfreya_interface.o \
	                nub.o
	$(COMPILE) nfreya_routines.f90

nfreya_write_statefile.o:	nfreya_write_statefile.f90 nrtype.o \
				neutral_beams.o zonal_data.o \
				plasma_prop.o io_gcnmp.o
	$(COMPILE) nfreya_write_statefile.f90


nblist_read.o:	nblist_read.f90 nrtype.o  beam_structure.o\
		io_gcnmp.o error_handler.o consts.o dnubeam_mod.o param.o
	$(COMPILE)  $<

io_gcnmp.o:                  io_gcnmp.f90  nrtype.o
	$(COMPILE)  $< 
  

ions_gcnmp.o:                ions_gcnmp.f90  nrtype.o
	$(COMPILE) $<  

iterdbmd_gcnmp.o:         iterdbmd_gcnmp.f90  nrtype.o 
	$(COMPILE)  $<

error_handler.o:         error_handler.F90 nrtype.o io_gcnmp.o MPI_data.o
	$(COMPILE)    $(GFLAGS)   $< 

gcnmp_version.o:	gcnmp_version.f90
	$(COMPILE) $<

grid_class.o:		grid_class.F90 plasma_prop.o vector_class.o \
			nrtype.o tspline.o error_handler.o 
	$(COMPILE)   $(GFLAGS)   $<


MPI_data.o:		     MPI_data.F90  nrtype.o consts.o solcon_gcnmp.o  
	$(COMPILE)    $(GFLAGS) $< 


nubeam_sub.o:	nubeam_sub.f90 nrtype.o
	$(COMPILE)  $<


neutral_beams.o:         neutral_beams.F90  nrtype.o ions_gcnmp.o nf_param.o \
			plasma_prop.o nrtype.o consts.o grid_class.o nf_param.o
	$(COMPILE) $(GFLAGS)  $<


neutral_data.o:	neutral_data.f90 nrtype.o           
	$(COMPILE)  $<

periodic_table_mod.o:	periodic_table_mod.f90 
	$(COMPILE)  $<

P_Nfreya_interface.o:	P_Nfreya_interface.F90 nrtype.o beam_structure.o \
			io_gcnmp.o  string_util.o consts.o           \
			solcon_gcnmp.o fast_ion_diffusion.o consts.o
	@echo ufile_lib = $(UFILE_LIBS)    
	@echo $(UFILE_LIBRARIES)    
	$(COMPILE) $(GFLAGS)  $<

prep_zonal_input.o:	prep_zonal_input.f90 contour_module.o ions_gcnmp.o \
			tport_mhd_grid_data.o zonal_data.o replace_imsl.o
	$(COMPILE)  $<

prep_transport_output.o:	prep_transport_output.F90 tspline.o \
			        nrtype.o nf_param.o zonal_data.o prep_zonal_input.o \
				grid_class.o consts.o tport_mhd_grid_data.o   \
				vector_class.o
	$(COMPILE) $(GFLAGS)  $<

xsct_class.o:	xsct_class.F90 nrtype.o zonal_data.o ions_gcnmp.o cpub_dat.o \
		hexnb_data.o d01bff_dat.o tspline.o replace_imsl.o
	$(COMPILE)  $(GFLAGS)  $<

tport_mhd_grid_data.o:	   tport_mhd_grid_data.f90 nrtype.o
	$(COMPILE)  $<

2d_interp.o:	2d_interp.f90
	$(COMPILE)  $<

1d_spline.o:	1d_spline.f90 nrtype.o
	$(COMPILE)  $<
       
nf_param.o:        nf_param.f90
	$(COMPILE)   $<

plasma_prop.o:         plasma_prop.f90 vector_class.o   nrtype.o 
	$(COMPILE)  $<

allocate_plasma_prop.o:         allocate_plasma_prop.f90 grid_class.o   nrtype.o \
		       	        solcon_gcnmp.o plasma_prop.o
	$(COMPILE) allocate_plasma_prop.f90


set_ion_prop2.o:		set_ion_prop2.f90 nrtype.o ions_gcnmp.o neutral_beams.o \
			error_handler.o consts.o
	$(COMPILE) $< 

contour_module.o:	contour_module.f90 nrtype.o consts.o io_gcnmp.o  \
			error_handler.o MPI_data.o  tport_mhd_grid_data.o \
			replace_imsl.o
	$(COMPILE) $< 



dep_var_gcnmp.o:         dep_var_gcnmp.f90  nrtype.o  
	$(COMPILE) $< 

fast_ion_data_gcnmp.o:         fast_ion_data_gcnmp.f90   
	$(COMPILE) $<

fast_ion_diffusion.o:	fast_ion_diffusion.f90 nrtype.o dnubeam_mod.o solcon_gcnmp.o
	$(COMPILE)  $<

fdyimp.o:          fdyimp.f90  nrtype.o grid_class.o
	$(COMPILE)  $<

glf23_gcnmp.o:         glf23_gcnmp.f90   nrtype.o
	$(COMPILE) $<

param.mod: param.o         
param.o:	param.f90
	echo $(VPATH)
	$(COMPILE) $< -o $@ 

mhdpar.mod:         mhdpar.o
mhdpar.o:	mhdpar.f90 param.o
	$(COMPILE) $< -o $@ 

bicube.mod:          bicube.o
bicube.o:          bicube.f90   cubic_spline.mod mhdpar.o
	$(COMPILE)  $< -o  $@

bicubic_spline.o:         bicubic_spline.f90  nrtype.o io_gcnmp.o error_handler.o
	$(COMPILE)  $(GFLAGS)  $< 

cubic_spline.mod: cubic_spline_class.o
cubic_spline_class.o:	cubic_spline_class.f90 nrtype.o common_constants.mod
	$(COMPILE) $<

m_refsor.mod:	m_refsor.o
m_refsor.o:	m_refsor.f90 nrtype.o 
	$(COMPILE)   $< -o  $@ 

replace_imsl.mod:	replace_imsl.o 
replace_imsl.o:		replace_imsl.F90 nrtype.o m_refsor.o \
			bicube.o cubic_spline_class.o  \
			bicubic_spline.o error_handler.o io_gcnmp.o
	$(COMPILE)   $< -o  $@ 


shot_info.o:         shot_info.f90   nrtype.o ions_gcnmp.o 
	$(COMPILE)  $<

solcon_gcnmp.o:		solcon_gcnmp.f90   ions_gcnmp.o io_gcnmp.o \
			glf23_gcnmp.o  nrtype.o   consts.o dep_var_gcnmp.o
	$(COMPILE) $<

source_terms_gcnmp.o:         source_terms_gcnmp.f90  nrtype.o grid_class.o \
			plasma_prop.o ions_gcnmp.o solcon_gcnmp.o vector_class.o \
			neutral_beams.o
	$(COMPILE)   $<

scale_profiles.o:	scale_profiles.f90 nrtype.o plasma_prop.o \
			ions_gcnmp.o consts.o grid_class.o  
	$(COMPILE) $<


string_util.o:	string_util.F90 nrtype.o grid_class.o
	$(COMPILE) $<

vector_class.o:		vector_class.f90   nrtype.o
	$(COMPILE) $<    


INCLUDES_dist_run  = nrtype.o solcon_gcnmp.o plasma_prop.o \
			   grid_class.o ions_gcnmp.o curden_terms.o bc_values_gcnmp.o    \
			   glf23_gcnmp.o   MPI_data.o consts.o zonal_data.o \
			   P_Nfreya_interface.o nub.o
distribute_run_directives.o:	distribute_run_directives.F90  $(INCLUDES_dis_run)
	$(COMPILE)    $(GFLAGS)   $<


INCLUDES_distribute_statefile  = nrtype.o  error_handler.o solcon_gcnmp.o     \
			     MPI_data.o consts.o
distribute_statefile.o:         distribute_statefile.F90  $(INCLUDES_distribute_statefile)
	$(COMPILE)    $(GFLAGS)   $<

distribute_ufile.o:	distribute_ufile.F90 nrtype.o MPI_data.o neutral_beams.o \
			P_Nfreya_interface.o ions_gcnmp.o
	$(COMPILE) $<    
#distribute_nubeam_namelist.o:	distribute_nubeam_namelist.f90 nbnamelist.o \
#			        nrtype.o  MPI_data.o  fast_ion_diffusion.o \
#				consts.o
# remove nbnamelist.o (not used any longer) HSJ 2/25/13
distribute_nubeam_namelist.o:	distribute_nubeam_namelist.f90 \
			        nrtype.o  MPI_data.o  fast_ion_diffusion.o \
				consts.o
	$(COMPILE) $<

dnubeam_mod.o:	dnubeam_mod.F90 nrtype.o param.o io_gcnmp.o
	$(COMPILE) $(GFLAGS)  $< 

INCLUDES_nfreya_driver = nrtype.o  error_handler.o  ions_gcnmp.o \
		         MPI_data.o nfreya_run.o ran_num_gen.o   \
			 neutral_beams.o nfreya_version.o        \
			 prep_transport_output.o nfreya_routines.o

nfreya_driver.o:         nfreya_driver.F90   $(INCLUDES_nfreya_driver)
	$(COMPILE)    $(GFLAGS)  nfreya_driver.F90

numbrs.o:	numbrs.f90 nrtype.o param.o
	$(COMPILE)   $<

g_nextiou.o:         g_nextiou.f90  
	$(COMPILE) g_nextiou.f90

getargs.o:         getargs.f90 io_gcnmp.o nrtype.o error_handler.o nfreya_version.o
	$(COMPILE) getargs.f90

get_ene_zeff.o:		get_ene_zeff.f90 nrtype.o grid_class.o ions_gcnmp.o fast_ion_data_gcnmp.o \
			consts.o dep_var_gcnmp.o neutral_beams.o
	$(COMPILE) get_ene_zeff.f90

get_tot_pressure.o:	get_tot_pressure.f90 nrtype.o grid_class.o ions_gcnmp.o fast_ion_data_gcnmp.o \
			consts.o dep_var_gcnmp.o neutral_beams.o
	$(COMPILE) $< -o $@


get_vars.o:         get_vars.f90 nrtype.o solcon_gcnmp.o grid_class.o \
			error_handler.o bc_values_gcnmp.o 
	$(COMPILE) get_vars.f90

load_nfreya_input.o:	load_nfreya_input.F90 iterdbmd_gcnmp.o io_gcnmp.o ions_gcnmp.o \
		        solcon_gcnmp.o grid_class.o bc_values_gcnmp.o nrtype.o \
	                error_handler.o MPI_data.o dep_var_gcnmp.o             \
			nfreya_namelist.o neutral_beams.o scale_profiles.o     \
			P_Nfreya_interface.o 
	$(COMPILE) $<


refactor_profiles.o:         refactor_profiles.f90 nrtype.o bc_values_gcnmp.o curden_terms.o \
			dep_var_gcnmp.o fast_ion_data_gcnmp.o fdyimp.o fluxx.o glf23_gcnmp.o \
			grid_class.o ions_gcnmp.o  plasma_prop.o solcon_gcnmp.o \
			source_terms_gcnmp.o tspline.o neo_tport.o
	$(COMPILE) refactor_profiles.f90


INCLUDES_rw_iterdb	=  error_handler.o  consts.o          \
			nrtype.o io_gcnmp.o  ions_gcnmp.o dep_var_gcnmp.o        \
                        shot_info.o plasma_prop.o solcon_gcnmp.o      \
                        source_terms_gcnmp.o fast_ion_data_gcnmp.o neutral_data.o \
                        curden_terms.o  grid_class.o bc_values_gcnmp.o      \
		        fdyimp.o MPI_data.o gcnmp_version.o iterdbmd_gcnmp.o tglfin.o \
			vector_class.o neutral_beams.o fusion_gcnmp.o tcoef.o \
			allocate_plasma_freq.o
rw_iterdb.o:         rw_iterdb.F90   $(INCLUDES_rw_iterdb)
	$(COMPILE)  $(GFLAGS)  $<

rw_iterdb_netcdf.o:         rw_iterdb_netcdf.F90 $(INCLUDES_rw_iterdb) 
	$(COMPILE) $(GFLAGS)  -I$(NETCDF_INC)   $<

ran_num_gen.o:	ran_num_gen.F90 nrtype.o consts.o
	$(COMPILE)  $(GFLAGS)  $<

reallocate.o:		reallocate.f90 nrtype.o
	$(COMPILE) reallocate.f90


r8tomsqz.o:        r8tomsqz.F f77_dcomplx.h
	$(COMPILE)  r8tomsqz.F

INCLUDES_savesol	= nrtype.o solcon_gcnmp.o grid_class.o               \
			  plasma_prop.o
savesol.o:         savesol.f90   $(INCLUDES_savesol)
	$(COMPILE) savesol.f90



set_ion_prop.o:	set_ion_prop.f90 nrtype.o ions_gcnmp.o grid_class.o         \
		error_handler.o plasma_prop.o                     \
		consts.o io_gcnmp.o neutral_beams.o
	$(COMPILE) set_ion_prop.f90 

INCLUDES_strip_comments = nrtype.o io_gcnmp.o error_handler.o \
			  beam_structure.o
strip_comments.o:         strip_comments.F90   $(INCLUDES_strip_comments)
	$(COMPILE) $<


tcoef.o:	tcoef.F90 nf_param.o nrtype.o
	$(COMPILE) $<

tglfv1.94_modules.o:	tglfv1.94_modules.F90
	$(COMPILE)  -DONETWO    $<  -o $@

tglfin.o:   tglfin.f90	nrtype.o tglfv1.94_modules.o  tglf_modules.o
	$(COMPILE) $(OPTIONS) $<  
	echo $(COMPILE)

tglf_modules.o:		tglf_modules.f90	
	$(COMPILE) $(OPTIONS) $< 
 
tspline.o:	tspline.f90 nrtype.o tspline.f90 error_handler.o nrtype.o
	$(COMPILE) $<

update.o:         update.f90  nrtype.o  solcon_gcnmp.o plasma_prop.o \
		  curden_terms.o vector_class.o tspline.o
	$(COMPILE) update.f90



utils.o:         utils.F90  nrtype.o  MPI_data.o neutral_beams.o \
	         ions_gcnmp.o plasma_prop.o solcon_gcnmp.o \
                 io_gcnmp.o error_handler.o 
	$(COMPILE)    $(GFLAGS)  $<


write_output.o:         write_output.f90  nrtype.o error_handler.o \
		   	 gcnmp_version.o dep_var_gcnmp.o
	$(COMPILE) write_output.f90



write_plotfile.o:	write_plotfile.f90 nrtype.o solcon_gcnmp.o    \
			dep_var_gcnmp.o plasma_prop.o grid_class.o    \
			io_gcnmp.o  gcnmp_version.o
	$(COMPILE) write_plotfile.f90

zonal_data.o:	zonal_data.f90 nrtype.o 
	$(COMPILE)   $<



$(LIBSM)/method_lines_library.a: $(LIBSM)/method_lines_library.o
	$(DELETE)  $(LIBSM)/method_lines_library.a
	$(ARCHIVE) $(LIBSM)/method_lines_library.a   $(LIBSM)/method_lines_library.o

$(LIBSM)/method_lines_library.o: $(LIBSM)/method_lines_library.f
	cd         $(LIBSM) ; \
	$(COMPILE) $(LIBSM)/method_lines_library.f


#$(LIBS)/$(IMSL).o:	$(LIBS)/$(IMSL).f
#	cd         $(LIBS) ; \
#	$(COMPILE) $(OPTIONS)    $<


#$(LIBS)/imsl_12.a:	$(LIBS)/$(IMSL).o
#	$(DELETE)  $(LIBS)/$(IMSL).a
#	$(ARCHIVE) $(LIBS)/$(IMSL).a      $(LIBS)/$(IMSL).o

.PHONY:	clean

clean:
	@echo cleaning up files on $(MACHINE)
	$(DELETE) $(NAME) *.a *.o *.mod *.lst
	#cd $(SHARED_MODULES)
	#$(DELETE) *.a *.o *.mod *.lst
	$(DELETE) $(LIBS)/*.o $(LIBS)/*.a 
	$(DELETE) $(LIBSM)/*.o $(LIBSM)/*.a 
	@echo $(LIBS)
	echo $(COMPILE)

print:
	echo $(VPATH)

