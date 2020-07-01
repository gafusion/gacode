include ${GACODE_ROOT}/platform/build/make.inc.${GACODE_PLATFORM}

ifeq ($(OPT),debug)
   FFLAGS=${FDEBUG}
else
   FFLAGS=${FOPT}
endif

FC += ${FOMP}

.SUFFIXES : .o .f90

all: geo_lib.a

geo_lib.a: geo.o
	$(ARCH) geo_lib.a geo.o

.f90.o :
	$(FC) $(FFLAGS) $(FMATH) -c $<

clean:
	rm -f *.o *.a *.mod
	rm -f ${GACODE_ROOT}/modules/geo.mod
