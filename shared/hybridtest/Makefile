include ${GACODE_ROOT}/platform/build/make.inc.${GACODE_PLATFORM}

ifeq ($(OPT),debug)
   FFLAGS=${FDEBUG}
else
   FFLAGS=${FOPT}
endif

FC += ${FOMP} ${FACC}

all:
	$(FC) $(FMATH) $(FFLAGS) -c hybrid.f90
	$(FC) $(FMATH) $(FFLAGS) -o hybrid hybrid.o $(LMATH)

clean:
	rm -f *.o *~ hybrid
