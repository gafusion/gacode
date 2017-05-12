# Either option 1 or 2
# 1
# The GACODE_ROOT environmental variable should point to the root directory where GACODE is installed.
# Environment variable GACODE_PLATFORM should be set to one of the following systems:
#   ALCF_BGP ALCF_CETUS BABBAGE BANACH CAOS CARVER CMODWS DELPHI
#   DROP EDISON_CRAY EDISON_INTEL FT GFORTRAN_CORE2 GFORTRAN_OSX
#   GFORTRAN_OSX_64 GFORTRAN_OSX_BELLI GFORTRAN_OSX_MACPORTS
#   GFORTRAN_OSX_TECHX GFORTRAN_PENRYN GFORTRAN_TECHX HOPPER_CRAY
#   HOPPER_INTEL HPC_ITER IFORT_CORE2 IFORT_PENRYN JAC LINDGREN
#   LOHAN LOKI LOKI_SCRATCH METIUS NEWT OSX_MOUNTAINLION PACER
#   PGFORTRAN_OSX PG_OPT64 PG_OPT64_FFTW PG_OPT64_MUMPS PPPL
#   PPPL_PATHSCALE RANGER SATURN TITAN_CRAY VENUS
# 2
# Define CC, CFLAGS, FC, FFLAGS, ARCH here or at the command line (make FC=FC, for instance)
#
ifdef GACODE_ROOT
	include $(GACODE_ROOT)/platform/build/make.inc.$(GACODE_PLATFORM)
else
	CC=cc
	CFLAGS=
	FC=gfortran
	FFLAGS=
	ARCH=ar cr
endif

LLIB = libharvest.a

EXEC = clientC clientF

OBJECTS = harvest_lib.o harvest_clientC.o harvest_clientF.o

$(LLIB): harvest_lib.o Makefile
	$(ARCH) $(LLIB) $<

clientC : harvest_client.c $(LLIB) 
	$(CC) $(CFLAGS) -o $@ -L./ $< -lharvest

clientF : harvest_client.F90 $(LLIB) 
	$(FC) $(FFLAGS) -o $@ -L./ $< -lharvest

all: $(LLIB) $(EXEC)

clean:
	rm -f *.o  *.a $(EXEC)

gacode_install:
	@[ -d $(GACODE_ROOT)/shared/harvest_client ] || mkdir $(GACODE_ROOT)/shared/harvest_client
	cp -f Makefile harvest_lib.c harvest_lib.h harvest_lib.inc harvest_lib.py $(GACODE_ROOT)/shared/harvest_client/

omfit_install:
	cp -f harvest_lib.py $(OMFIT_ROOT)/src/classes/
