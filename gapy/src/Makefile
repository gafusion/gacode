include ${GACODE_ROOT}/platform/build/make.inc.${GACODE_PLATFORM}

ifeq ($(OPT),debug)
   FFLAGS=${FDEBUG}
else
   FFLAGS=${FOPT}
endif

FC += ${FOMP}

.SUFFIXES : .o .f90

all: geo_lib.a util_lib.a

geo_lib.a: gapy_geo.o
	$(ARCH) geo_lib.a gapy_geo.o

#vpro_lib.a: gapy_vpro.o
#	$(ARCH) vpro_lib.a gapy_vpro.o

util_lib.a: gapy_util.o
	$(ARCH) util_lib.a gapy_util.o

.f90.o :
	$(FC) $(FFLAGS) $(FMATH) -c $<

clean:
	rm -f *.o *.a *.mod
