include ${GACODE_ROOT}/platform/build/make.inc.${GACODE_PLATFORM}

ifeq ($(OPT),debug)
   FFLAGS=${FDEBUG}
else
   FFLAGS=${FOPT}
endif

FC += ${FOMP}

OBJECTS = expro.o \
          expro_util.o \
          expro_comm.o \
          expro_locsim.o

.SUFFIXES : .o .f90

all: expro_lib.a

expro_lib.a: $(OBJECTS)
	$(ARCH) expro_lib.a $(OBJECTS)

.f90.o .F90.o :
	$(FC) $(FFLAGS) $(FMATH) -c $<

clean:
	rm -f *.o *.a *.mod
	rm -f ${GACODE_ROOT}/modules/expr*.mod
