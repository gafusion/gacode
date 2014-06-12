#USE THIS MAKEFILE TO BUILD ONETWO On LOCAL ARCHITECTURES
# FORT=pgf90,FORT=ifort  are  valid,also OPT=optimize (default) or OPT=debug
#or OPT=debug  may be set on all machines
#gmake FORT=pgf95 OPT=optimize


VERSION   = v5.4
GRID      = 129_201
NAME      = onetwo_$(GRID)
HOST      = $(shell uname -s)
ONETWO    = $(shell pwd)
SHARED_MODULES = $(ONETWO)/shared_modules
VPATH = $(SHARED_MODULES)
vpath nfreya_module

include ../shared/install/make.inc.$(GACODE_PLATFORM)

PROTECT   = chmod 755
DELETE    = rm -f
RENAME    = mv

include ../shared/install/make.inc.$(GACODE_PLATFORM)
#These should be in the preceding platform dependent included makefile
TSPLINE_LIB   = /p/linux/tools/lib
SPLINE_LIB    = splinelib_pgf_64.a
LAPACK = -L/c/source/PGI/pgi/linux86-64/13.2/lib -llapack
UFILE_LIBS =  /p/linux/nubeam/201405/ntcc/lib
UFILE_LIBRARIES  = $(UFILE_LIBS)/libuflib.a         \
$(UFILE_LIBS)/libufhdf.a        \
$(UFILE_LIBS)/libmdstransp.a    \
$(UFILE_LIBS)/libvaxonly.a      \
$(UFILE_LIBS)/libportlib.a      \
$(UFILE_LIBS)/libmds_sub.a      \
$(UFILE_LIBS)/libmdstransp.a /f/mdsplus/linux64/lib64/libMdsLib.so
FC = /opt/pgi/linux86-64/13.2/bin/pgf90
#FC =  /opt/pgi/linux86-64/11.10/bin/pgf90
COMPILE   = $(FC) $(OMP) $(PROF)\
  $(AUTP) -c $(OPTION) $(LISTING) $(SPECIAL)

#End platform dependent
LOAD      = $(FC)  -o $(NAME) 
include ./make.objects # Defines OBJECTS

LIBS = onetwo_libs
PORTLIB = portlib

LIBRARIES = onetwo_libs/portlib.a onetwo_libs/cortlib.a   \
	    onetwo_libs/method_lines_library.a  $(LAPACK) \
	    $(UFILE_LIBRARIES)  $(TSPLINE_LIB)/$(SPLINE_LIB)


.SUFFIXES:


all: $(NAME) 
	@echo
	@echo ">>  `date '+%a %d-%h-%y %r'`  `pwd`  `uname -mns`  $(LOGNAME)"
	@echo
	@echo make of $(NAME) in `pwd` on `hostname` completed
#	@echo $(UFILE_LIBRARIES)
#	@echo $(NETCDF_LIB)

$(NAME):           $(OBJECTS) $(OBJECTS_P) $(LIBRARIES)
	@echo "loading .... "
	$(LOAD)  -L /u/meneghini/fann/src/  $(OBJECTS) $(NETCDF_LIB) $(LIBRARIES) $(FTN90LIB) -ldoublefann
	$(PROTECT) $(NAME)
#	cp $(NAME) $(INSTALLDIR)

brainfuse.o: brainfuse.c
	gcc -c -g -I/u/meneghini/fann/src/include brainfuse.c

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

.dep_incs/%.inc: %
	@./make_depends.sh $< $@
make.includes: make.objects
	@python ./make_includes.py -path $(VPATH) -obj $(OBJECTS) $(OBJECTS_P)
include make.includes

#make_depends: make_depends.sh make.objects
#	@echo "Running make_depends.sh."
#	@./make_depends.sh "$(OBJECTS)"
#include make_depends

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
	cd         $(LIBS) ; \
	$(PGCC)  -c    $(LIBS)/cortlib.c

$(LIBS)/$(PORTLIB).a: $(LIBS)/$(PORTLIB).o
	$(DELETE)  $(LIBS)/$(PORTLIB).a
	$(ARCHIVE) $(LIBS)/$(PORTLIB).a   $(LIBS)/$(PORTLIB).o
$(LIBS)/$(PORTLIB).o: $(LIBS)/$(PORTLIB).F
	cd         $(LIBS) ; \
	$(COMPILE) -I../ $(LIBS)/$(PORTLIB).F


.PHONY:	clean
clean:
	@echo cleaning up files on $(MACHINE)
	 $(DELETE) $(NAME) *.a *.o *.mod
	 $(DELETE) $(LIBS)/*.o $(LIBS)/*.a 
	 $(DELETE) $(SHARED_MODULES)/*.o
	 $(DELETE) .dep_incs/*.inc
	 #@echo $(LIBS)

#data_files: 
#	cd $(ONETWO)

#test_piksr2: test_piksr2.o 
#	$(LOAD) $^ -o test_piksr2
