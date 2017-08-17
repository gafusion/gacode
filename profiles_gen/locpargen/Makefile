include ${GACODE_ROOT}/platform/build/make.inc.${GACODE_PLATFORM}

ifeq ($(OPT),debug)
   FFLAGS=${FDEBUG}
else
   FFLAGS=${FOPT}
endif

EXEC = locpargen

EXTRA_LIBS = ${GACODE_ROOT}/shared/EXPRO/EXPRO_lib.a \
             ${GACODE_ROOT}/shared/GEO/GEO_lib.a  \
             ${GACODE_ROOT}/shared/math/math_lib.a 

OBJECTS = $(EXEC).o

.SUFFIXES : .o .f90

exec: $(OBJECTS) $(EXTRA_LIBS)
	$(FC) $(FFLAGS) -o $(EXEC) $(OBJECTS) $(EXTRA_LIBS) $(LMATH)

.f90.o :
	$(FC) $(FMATH) $(FFLAGS) -c $<

clean:
	rm -f *.o $(EXEC)
	rm -f ${GACODE_ROOT}/modules/locpargen*.mod
