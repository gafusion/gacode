include ${GACODE_ROOT}/platform/build/make.inc.${GACODE_PLATFORM}

ifeq ($(OPT),debug)
   FFLAGS=${FDEBUG}
else
   FFLAGS=${FOPT}
endif

vgen_libs = -L$(GACODE_ROOT)/shared/harvest_client/ \
              -I$(GACODE_ROOT)/shared/harvest_client/ \
              -lharvest \
              
ifdef FANN_ROOT              
        JBSNN_LIBS=-I$(GACODE_ROOT)/shared/neural/ -L$(GACODE_ROOT)/shared/neural/ -lbrainfuse -L${FANN_ROOT}/lib -lfann
endif

EXEC = vgen

OBJECTS = vgen_globals.o \
	  vgen_init.o \
	  vgen_harvest.o \
	  vgen_compute_neo.o \
	  vgen_getgeo.o \
	  vgen_reduce.o \
	  $(EXEC).o

.SUFFIXES : .o .f90

all: $(OBJECTS) $(EXTRA_LIBS)
	$(FC) $(FFLAGS) -o $(EXEC) $(OBJECTS) $(EXTRA_LIBS) $(LMATH) $(JBSNN_LIBS) $(vgen_libs)

.f90.o :
	$(FC) $(FMATH) $(FFLAGS) $(vgen_libs) -c $<

clean:
	rm -f *.o $(EXEC) ${GACODE_ROOT}/modules/vgen_globals.mod
