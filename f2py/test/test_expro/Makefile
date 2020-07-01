include ${GACODE_ROOT}/platform/build/make.inc.${GACODE_PLATFORM}

ifeq ($(OPT),debug)
   FFLAGS=${FDEBUG}
else
   FFLAGS=${FOPT}
endif

FC += ${FOMP}

.SUFFIXES : .o .f90

EXTRA_LIBS = \
	${GACODE_ROOT}/f2py/expro/expro_lib.a \
	${GACODE_ROOT}/f2py/geo/geo_lib.a

test: test.o
	$(FC) -o test test.o ${EXTRA_LIBS} ${LMATH}

.f90.o :
	$(FC) $(FFLAGS) $(FMATH) -c $<

clean:
	rm -f *.o test input.gacode
