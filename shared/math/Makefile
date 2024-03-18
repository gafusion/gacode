include ${GACODE_ROOT}/platform/build/make.inc.${GACODE_PLATFORM}

ifeq ($(OPT),debug)
   FFLAGS=${FDEBUG}
else
   FFLAGS=${FOPT}
endif

LLIB=math_lib

OBJECTS = math_constants.o \
          gauss_legendre.o \
          half_hermite.o \
          i0.o \
          parallel_dim.o \
          pseudospec_lib.o

.SUFFIXES : .o .f90 .F90 .f

all: lib

lib: $(OBJECTS)
	$(ARCH) $(LLIB).a $(OBJECTS)

.f.o :
	$(F77) $(FFLAGS) -c $<

.f90.o :
	$(FC) $(FMATH) $(FFLAGS) -c $<

clean:
	rm -f *.o $(LLIB).a *.mod *.lst
	rm -f ${GACODE_ROOT}/modules/half_hermite.mod
	rm -f ${GACODE_ROOT}/modules/math_constants.mod
	rm -f ${GACODE_ROOT}/modules/parallel_lib.mod

-include ${GACODE_ROOT}/platform/build/make.inc.${GACODE_PLATFORM}
