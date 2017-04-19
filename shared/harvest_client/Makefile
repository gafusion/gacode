# Either option 1 or 2
# 1. The GACODE_ROOT environmental variable should point to the root 
#    directory where GACODE is installed.
#    The GACODE_PLATFORM variable should be set to one of the systems listed
#    in platform/install
#
# 2. Define CC, CFLAGS, FC, FFLAGS, ARCH here or at the command line
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

clientF : harvest_client.f90 $(LLIB) 
	$(FC) $(FFLAGS) -o $@ -L./ $< -lharvest

all: $(LLIB) $(EXEC)

clean:
	rm -f *.o  *.a $(EXEC)

gacode_install:
	@[ -d $(GACODE_ROOT)/shared/harvest_client ] || mkdir $(GACODE_ROOT)/shared/harvest_client
	cp -f Makefile harvest_lib.c harvest_lib.h harvest_lib.inc $(GACODE_ROOT)/shared/harvest_client/

omfit_install:
	cp -f harvest_lib.py $(OMFIT_ROOT)/src/classes/
