include ${GACODE_ROOT}/shared/install/make.inc.${GACODE_PLATFORM}

ifeq ($(OPT),debug)
   FFLAGS=${FDEBUG}
else
   FFLAGS=${FOPT}
endif

vgen_libs = -L$(GACODE_ROOT)/shared/harvest_client/ \
              -I$(GACODE_ROOT)/shared/harvest_client/ \
              -lharvest     
EXEC = vgen

OBJECTS = vgen_globals.o \
	  vgen_init.o \
	  vgen_compute_neo.o \
	  vgen_getgeo.o \
	  vgen_reduce.o \
	  $(EXEC).o

.SUFFIXES : .o .f90

all: $(OBJECTS) $(EXTRA_LIBS)
	$(FC) $(FFLAGS) $(vgen_libs) -o $(EXEC) $(OBJECTS) $(EXTRA_LIBS) $(LMATH)

.f90.o :
	$(FC) $(FMATH) $(FFLAGS) $(vgen_libs) -c $<

clean:
	rm -f *.o $(EXEC) ${GACODE_ROOT}/modules/vgen_globals.mod
