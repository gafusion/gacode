include ${GACODE_ROOT}/platform/build/make.inc.${GACODE_PLATFORM}

ifeq ($(OPT),debug)
   FFLAGS=${FDEBUG}
else
   FFLAGS=${FOPT}
endif

EXEC = locpargen

EXTRA_LIBS = ${GACODE_ROOT}/f2py/expro/expro_lib.a \
             ${GACODE_ROOT}/f2py/geo/geo_lib.a 

OBJECTS = locpargen_globals.o \
          locpargen_geo.o \
          locpargen_neo.o \
          locpargen_tglf.o \
	  locpargen_tglf_stack.o \
          locpargen_cgyro.o \
          $(EXEC).o 

.SUFFIXES : .o .f90

exec: $(OBJECTS) $(EXTRA_LIBS)
	$(FC) $(FFLAGS) -o $(EXEC) $(OBJECTS) $(EXTRA_LIBS) $(LMATH)

.f90.o :
	$(FC) $(FMATH) $(FFLAGS) -c $<

clean:
	rm -f *.o $(EXEC) *~
	rm -f ${GACODE_ROOT}/modules/locpargen*.mod
