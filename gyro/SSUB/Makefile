include ${GACODE_ROOT}/platform/build/make.inc.${GACODE_PLATFORM}

ifeq ($(OPT),debug)
   FFLAGS=${FDEBUG}
else
   FFLAGS=${FOPT}
endif

FC += ${FOMP}

LLIB=SSUB_lib

OBJECTS = SSUB_private.o \
          fSSUB.o \
          rSSUB.o \
          SSUB_init.o \
          SSUB_cleanup.o

.SUFFIXES : .o .f90 

lib: $(OBJECTS)
	$(ARCH) $(LLIB).a $(OBJECTS)

.f90.o :
	$(FC) $(FMATH) $(FFLAGS) -c $<

clean:
	rm -rf *.o $(LLIB).a
	rm -rf ${GACODE_ROOT}/modules/ssub*.mod
                
